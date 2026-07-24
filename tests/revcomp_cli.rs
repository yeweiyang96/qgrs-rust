use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use arrow_array::{Array, Int32Array, StringArray, UInt64Array};
use arrow_schema::DataType;
use flate2::Compression;
use flate2::write::GzEncoder;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

const CSV_HEADER: &str = "start,end,length,tetrads,y1,y2,y3,score,sequence";
const ASYMMETRIC_G4: &str = "GGGGAGGGGTTGGGGCCCGGGG";
const ASYMMETRIC_G4_SOURCE: &str = "TTCCCCGGGCCCCAACCCCTCCCCAAA";
const ASYMMETRIC_I_MOTIF: &str = "CCCCTCCCCTTCCCCTTTCCCC";
const ASYMMETRIC_I_MOTIF_SOURCE: &str = "AAGGGGAAAGGGGAAGGGGAGGGGTTT";
const OVERLAPPING_G4_SOURCE: &str = "CCCCCGGGCCCCCAACCCCCTCCCCC";
const REVERSE_COMPLEMENT_BLOCK_SIZE: usize = 64 * 1024;

static WORKSPACE_COUNTER: AtomicU64 = AtomicU64::new(0);

struct TestWorkspace {
    root: PathBuf,
}

#[derive(Debug, PartialEq, Eq)]
struct CsvTable {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

#[derive(Debug)]
struct ParquetTable {
    header: Vec<String>,
    data_types: Vec<DataType>,
    rows: Vec<Vec<String>>,
}

impl TestWorkspace {
    fn new(label: &str) -> Self {
        let counter: u64 = WORKSPACE_COUNTER.fetch_add(1, Ordering::Relaxed);
        let timestamp: u128 = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock must be after the Unix epoch")
            .as_nanos();
        let root: PathBuf = std::env::temp_dir().join(format!(
            "qgrs-revcomp-cli-{}-{timestamp}-{counter}-{label}",
            std::process::id()
        ));
        fs::create_dir(&root).expect("create isolated revcomp CLI test workspace");
        Self { root }
    }

    fn path(&self, name: &str) -> PathBuf {
        self.root.join(name)
    }
}

impl Drop for TestWorkspace {
    fn drop(&mut self) {
        if let Err(error) = fs::remove_dir_all(&self.root) {
            eprintln!(
                "failed to clean revcomp CLI test workspace: path={:?} error={error}",
                self.root
            );
        }
    }
}

fn qgrs_command() -> Command {
    Command::new(env!("CARGO_BIN_EXE_qgrs"))
}

fn run(command: &mut Command) -> Output {
    command.output().expect("execute qgrs CLI")
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "qgrs failed with status {:?}\nstdout:\n{}\nstderr:\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_failure(output: &Output) {
    assert!(
        !output.status.success(),
        "qgrs unexpectedly succeeded\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn csv_columns(csv: &str) -> Vec<&str> {
    assert_eq!(csv.lines().next(), Some(CSV_HEADER));
    let rows: Vec<&str> = csv.lines().skip(1).collect();
    assert_eq!(rows.len(), 1, "expected exactly one CSV result:\n{csv}");
    rows[0].split(',').collect()
}

fn write_gzip(path: &Path, bytes: &[u8]) {
    let file: fs::File = fs::File::create(path).expect("create gzip FASTA");
    let mut encoder: GzEncoder<fs::File> = GzEncoder::new(file, Compression::default());
    encoder.write_all(bytes).expect("write gzip FASTA");
    encoder.finish().expect("finish gzip FASTA");
}

fn directory_snapshot(path: &Path) -> Vec<(String, Vec<u8>)> {
    let mut files: Vec<(String, Vec<u8>)> = fs::read_dir(path)
        .expect("read output directory")
        .map(|entry_result| {
            let entry: fs::DirEntry = entry_result.expect("read output directory entry");
            let file_type: fs::FileType = entry.file_type().expect("read output entry type");
            assert!(file_type.is_file(), "unexpected non-file output: {entry:?}");
            let name: String = entry
                .file_name()
                .into_string()
                .expect("output filename must be valid UTF-8");
            let bytes: Vec<u8> = fs::read(entry.path()).expect("read output file");
            (name, bytes)
        })
        .collect();
    files.sort_by(|left, right| left.0.cmp(&right.0));
    files
}

fn assert_directory_empty(path: &Path) {
    let mut entries: fs::ReadDir = fs::read_dir(path).expect("read temporary directory");
    assert!(
        entries.next().is_none(),
        "temporary directory was not cleaned: {path:?}"
    );
}

fn csv_row_count(bytes: &[u8]) -> usize {
    std::str::from_utf8(bytes)
        .expect("CSV output must be valid UTF-8")
        .lines()
        .skip(1)
        .count()
}

fn read_csv_table(path: &Path) -> CsvTable {
    let mut reader: csv::Reader<fs::File> = csv::Reader::from_path(path).expect("open CSV table");
    let header: Vec<String> = reader
        .headers()
        .expect("read CSV header")
        .iter()
        .map(str::to_string)
        .collect();
    let rows: Vec<Vec<String>> = reader
        .records()
        .map(|record_result| {
            record_result
                .expect("read CSV row")
                .iter()
                .map(str::to_string)
                .collect()
        })
        .collect();
    CsvTable { header, rows }
}

fn read_parquet_table(path: &Path) -> ParquetTable {
    let file: fs::File = fs::File::open(path).expect("open Parquet table");
    let builder: ParquetRecordBatchReaderBuilder<fs::File> =
        ParquetRecordBatchReaderBuilder::try_new(file).expect("read Parquet metadata");
    let schema = builder.schema().clone();
    let header: Vec<String> = schema
        .fields()
        .iter()
        .map(|field| field.name().to_string())
        .collect();
    let data_types: Vec<DataType> = schema
        .fields()
        .iter()
        .map(|field| field.data_type().clone())
        .collect();
    let mut rows: Vec<Vec<String>> = Vec::new();
    let reader = builder.build().expect("build Parquet record batch reader");
    for batch_result in reader {
        let batch = batch_result.expect("read Parquet record batch");
        for row_index in 0..batch.num_rows() {
            let row: Vec<String> = batch
                .columns()
                .iter()
                .map(|column| match column.data_type() {
                    DataType::UInt64 => column
                        .as_any()
                        .downcast_ref::<UInt64Array>()
                        .expect("UInt64 Parquet column")
                        .value(row_index)
                        .to_string(),
                    DataType::Int32 => column
                        .as_any()
                        .downcast_ref::<Int32Array>()
                        .expect("Int32 Parquet column")
                        .value(row_index)
                        .to_string(),
                    DataType::Utf8 => column
                        .as_any()
                        .downcast_ref::<StringArray>()
                        .expect("Utf8 Parquet column")
                        .value(row_index)
                        .to_string(),
                    data_type => panic!("unexpected Parquet data type: {data_type:?}"),
                })
                .collect();
            rows.push(row);
        }
    }
    ParquetTable {
        header,
        data_types,
        rows,
    }
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|byte| match byte.to_ascii_uppercase() {
            b'A' => 'T',
            b'T' | b'U' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            byte => panic!("unexpected test nucleotide: {}", char::from(byte)),
        })
        .collect()
}

fn table_column_index(table: &CsvTable, name: &str) -> usize {
    table
        .header
        .iter()
        .position(|column| column == name)
        .unwrap_or_else(|| panic!("missing CSV column: {name}"))
}

fn project_linear_result_table(mut table: CsvTable, sequence_len: usize) -> CsvTable {
    let start_index: usize = table_column_index(&table, "start");
    let end_index: usize = table_column_index(&table, "end");
    for row in &mut table.rows {
        let start_rc: usize = row[start_index].parse().expect("parse RC start");
        let end_rc: usize = row[end_index].parse().expect("parse RC end");
        row[start_index] = (sequence_len - end_rc + 1).to_string();
        row[end_index] = (sequence_len - start_rc + 1).to_string();
    }
    table.rows.sort_by_key(|row| {
        (
            row[start_index].parse::<usize>().expect("parse start"),
            row[end_index].parse::<usize>().expect("parse end"),
        )
    });
    table
}

fn project_linear_family_table(mut table: CsvTable, sequence_len: usize) -> CsvTable {
    let family_index: usize = table_column_index(&table, "family_index");
    let start_index: usize = table_column_index(&table, "start");
    let end_index: usize = table_column_index(&table, "end");
    for row in &mut table.rows {
        let start_rc: usize = row[start_index].parse().expect("parse RC family start");
        let end_rc: usize = row[end_index].parse().expect("parse RC family end");
        row[start_index] = (sequence_len - end_rc + 1).to_string();
        row[end_index] = (sequence_len - start_rc + 1).to_string();
    }
    table.rows.sort_by_key(|row| {
        (
            row[start_index].parse::<usize>().expect("parse start"),
            row[end_index].parse::<usize>().expect("parse end"),
        )
    });
    for (index, row) in table.rows.iter_mut().enumerate() {
        row[family_index] = (index + 1).to_string();
    }
    table
}

fn snapshot_matching(
    snapshot: &[(String, Vec<u8>)],
    predicate: impl Fn(&str) -> bool,
) -> Vec<(String, Vec<u8>)> {
    snapshot
        .iter()
        .filter(|(name, _bytes)| predicate(name))
        .cloned()
        .collect()
}

fn fasta_scan_command(input: &Path, output_dir: &Path, mode: &str) -> Command {
    let mut command: Command = qgrs_command();
    command.args([
        "--file",
        input.to_str().expect("FASTA path must be valid UTF-8"),
        "--mode",
        mode,
        "--output-dir",
        output_dir
            .to_str()
            .expect("output directory path must be valid UTF-8"),
        "--min-tetrads",
        "4",
        "--min-score",
        "17",
        "--overlap",
    ]);
    command
}

fn run_fasta_scan(input: &Path, output_dir: &Path, mode: &str) -> Output {
    run(&mut fasta_scan_command(input, output_dir, mode))
}

fn run_fasta_revcomp_scan(input: &Path, output_dir: &Path, mode: &str) -> Output {
    let mut command: Command = fasta_scan_command(input, output_dir, mode);
    command.arg("--revcomp");
    run(&mut command)
}

fn run_circular_fasta_revcomp_scan(input: &Path, output_dir: &Path, mode: &str) -> Output {
    let mut command: Command = fasta_scan_command(input, output_dir, mode);
    command.args(["--circular", "--revcomp"]);
    run(&mut command)
}

fn assert_linear_revcomp_matches_explicit_scan(
    workspace: &TestWorkspace,
    source: &str,
    label: &str,
) -> CsvTable {
    let forward_path: PathBuf = workspace.path(&format!("{label}.csv"));
    let actual_output: Output = run(qgrs_command().args([
        "--sequence",
        source,
        "--min-tetrads",
        "4",
        "--output",
        forward_path
            .to_str()
            .expect("output path must be valid UTF-8"),
        "--overlap",
        "--revcomp",
    ]));
    assert_success(&actual_output);

    let explicit_rc_path: PathBuf = workspace.path(&format!("{label}-explicit-rc.csv"));
    let explicit_rc: String = reverse_complement(source);
    let reference_output: Output = run(qgrs_command().args([
        "--sequence",
        &explicit_rc,
        "--min-tetrads",
        "4",
        "--output",
        explicit_rc_path
            .to_str()
            .expect("explicit RC output path must be valid UTF-8"),
        "--overlap",
    ]));
    assert_success(&reference_output);

    let expected_primary: CsvTable =
        project_linear_result_table(read_csv_table(&explicit_rc_path), source.len());
    let expected_overlap: CsvTable = project_linear_result_table(
        read_csv_table(&workspace.path(&format!("{label}-explicit-rc.overlap.csv"))),
        source.len(),
    );
    let expected_family: CsvTable = project_linear_family_table(
        read_csv_table(&workspace.path(&format!("{label}-explicit-rc.family.csv"))),
        source.len(),
    );
    let actual_primary: CsvTable = read_csv_table(&workspace.path(&format!("{label}.revcomp.csv")));
    assert_eq!(actual_primary, expected_primary);
    assert_eq!(
        read_csv_table(&workspace.path(&format!("{label}.revcomp.overlap.csv"))),
        expected_overlap
    );
    assert_eq!(
        read_csv_table(&workspace.path(&format!("{label}.revcomp.family.csv"))),
        expected_family
    );
    actual_primary
}

#[test]
fn default_scan_does_not_create_revcomp_output() {
    let workspace: TestWorkspace = TestWorkspace::new("default");
    let forward_path: PathBuf = workspace.path("result.csv");
    let reverse_path: PathBuf = workspace.path("result.revcomp.csv");

    let output: Output = run(qgrs_command().args([
        "--sequence",
        ASYMMETRIC_G4_SOURCE,
        "--min-tetrads",
        "4",
        "--output",
        forward_path
            .to_str()
            .expect("output path must be valid UTF-8"),
    ]));

    assert_success(&output);
    assert!(forward_path.is_file());
    assert!(!reverse_path.exists());
}

#[test]
fn inline_revcomp_requires_output_and_maps_asymmetric_hit() {
    let workspace: TestWorkspace = TestWorkspace::new("inline");

    let missing_output: Output = run(qgrs_command().args([
        "--sequence",
        ASYMMETRIC_G4_SOURCE,
        "--min-tetrads",
        "4",
        "--revcomp",
    ]));
    assert_failure(&missing_output);
    assert!(
        String::from_utf8_lossy(&missing_output.stderr)
            .contains("--revcomp requires a file path via --output")
    );

    let stdout_sentinel: Output = run(qgrs_command().current_dir(&workspace.root).args([
        "--sequence",
        ASYMMETRIC_G4_SOURCE,
        "--min-tetrads",
        "4",
        "--output",
        "-",
        "--revcomp",
    ]));
    assert_failure(&stdout_sentinel);
    assert!(
        String::from_utf8_lossy(&stdout_sentinel.stderr)
            .contains("--revcomp requires a file path via --output")
    );
    assert!(!workspace.path("-").exists());

    let asymmetric: CsvTable =
        assert_linear_revcomp_matches_explicit_scan(&workspace, ASYMMETRIC_G4_SOURCE, "result");
    assert_eq!(asymmetric.rows.len(), 1);
    let columns: Vec<&str> = asymmetric.rows[0].iter().map(String::as_str).collect();
    assert_eq!(columns.len(), 9);
    assert_eq!(&columns[0..4], ["3", "24", "22", "4"]);
    assert_eq!(&columns[4..7], ["1", "2", "3"]);
    assert_eq!(columns[7], "82");
    assert_eq!(columns[8], ASYMMETRIC_G4);

    let overlapping: CsvTable =
        assert_linear_revcomp_matches_explicit_scan(&workspace, OVERLAPPING_G4_SOURCE, "overlap");
    let overlapping_raw: CsvTable = read_csv_table(&workspace.path("overlap.revcomp.overlap.csv"));
    assert!(
        overlapping_raw.rows.len() > overlapping.rows.len(),
        "overlapping fixture must exercise raw-hit family consolidation"
    );
}

#[test]
fn base_c_revcomp_emits_negative_strand_i_motif() {
    let workspace: TestWorkspace = TestWorkspace::new("base-c");
    let forward_path: PathBuf = workspace.path("result.csv");
    let reverse_path: PathBuf = workspace.path("result.revcomp.csv");

    let output: Output = run(qgrs_command().args([
        "--sequence",
        ASYMMETRIC_I_MOTIF_SOURCE,
        "--base",
        "c",
        "--min-tetrads",
        "4",
        "--output",
        forward_path
            .to_str()
            .expect("output path must be valid UTF-8"),
        "--revcomp",
    ]));
    assert_success(&output);

    let reverse_csv: String =
        fs::read_to_string(&reverse_path).expect("read i-motif revcomp CSV output");
    let columns: Vec<&str> = csv_columns(&reverse_csv);
    assert_eq!(&columns[0..4], ["3", "24", "22", "4"]);
    assert_eq!(&columns[4..7], ["1", "2", "3"]);
    assert_eq!(columns[8], ASYMMETRIC_I_MOTIF);
}

#[test]
fn parquet_revcomp_outputs_match_csv_schema_and_rows() {
    let workspace: TestWorkspace = TestWorkspace::new("parquet");
    let forward_path: PathBuf = workspace.path("result.parquet");
    let csv_path: PathBuf = workspace.path("reference.csv");

    let csv_output: Output = run(qgrs_command().args([
        "--sequence",
        ASYMMETRIC_G4_SOURCE,
        "--min-tetrads",
        "4",
        "--output",
        csv_path
            .to_str()
            .expect("CSV output path must be valid UTF-8"),
        "--overlap",
        "--revcomp",
    ]));
    assert_success(&csv_output);

    let output: Output = run(qgrs_command().args([
        "--sequence",
        ASYMMETRIC_G4_SOURCE,
        "--min-tetrads",
        "4",
        "--format",
        "parquet",
        "--output",
        forward_path
            .to_str()
            .expect("output path must be valid UTF-8"),
        "--overlap",
        "--revcomp",
    ]));
    assert_success(&output);

    let result_types: Vec<DataType> = vec![
        DataType::UInt64,
        DataType::UInt64,
        DataType::UInt64,
        DataType::UInt64,
        DataType::Int32,
        DataType::Int32,
        DataType::Int32,
        DataType::Int32,
        DataType::Utf8,
    ];
    let family_types: Vec<DataType> = vec![DataType::UInt64, DataType::UInt64, DataType::UInt64];
    for (parquet_name, csv_name, expected_types) in [
        (
            "result.revcomp.parquet",
            "reference.revcomp.csv",
            &result_types,
        ),
        (
            "result.revcomp.overlap.parquet",
            "reference.revcomp.overlap.csv",
            &result_types,
        ),
        (
            "result.revcomp.family.parquet",
            "reference.revcomp.family.csv",
            &family_types,
        ),
    ] {
        let parquet: ParquetTable = read_parquet_table(&workspace.path(parquet_name));
        let csv: CsvTable = read_csv_table(&workspace.path(csv_name));
        assert_eq!(
            parquet.header, csv.header,
            "schema names for {parquet_name}"
        );
        assert_eq!(
            &parquet.data_types, expected_types,
            "schema types for {parquet_name}"
        );
        assert_eq!(parquet.rows, csv.rows, "row values for {parquet_name}");
    }
}

#[test]
fn fasta_revcomp_outputs_match_across_modes_and_preserve_forward_files() {
    let workspace: TestWorkspace = TestWorkspace::new("fasta-modes");
    let plain_path: PathBuf = workspace.path("input.fa");
    let gzip_path: PathBuf = workspace.path("input.fa.gz");
    let block_sequence_len: usize = 70_000;
    let reverse_block_start: usize = block_sequence_len - REVERSE_COMPLEMENT_BLOCK_SIZE;
    let long_run_len: usize = 96;
    let long_run_start: usize = reverse_block_start - 512;
    let motif_start: usize = reverse_block_start - 13;
    let between_len: usize = motif_start - long_run_start - long_run_len;
    let block_suffix_len: usize = block_sequence_len - motif_start - OVERLAPPING_G4_SOURCE.len();
    assert!(long_run_len > 10, "fixture must exceed the default max run");
    assert!(
        long_run_start + long_run_len < reverse_block_start,
        "long target-base run must be fully scanned before the block boundary"
    );
    assert!(motif_start < reverse_block_start);
    assert!(
        motif_start + OVERLAPPING_G4_SOURCE.len() > reverse_block_start,
        "asymmetric motif must cross the reverse-read block boundary"
    );
    let block_boundary_sequence: String = format!(
        "{}{}{}{OVERLAPPING_G4_SOURCE}{}",
        "A".repeat(long_run_start),
        "C".repeat(long_run_len),
        "A".repeat(between_len),
        "A".repeat(block_suffix_len),
    );
    assert_eq!(block_boundary_sequence.len(), block_sequence_len);
    let fasta: String = format!(
        ">negative_only\n{ASYMMETRIC_G4_SOURCE}\n\
         >both_strands\n{ASYMMETRIC_G4}{}{ASYMMETRIC_G4_SOURCE}\n\
         >block_boundary\n{block_boundary_sequence}\n",
        "A".repeat(50),
    );
    fs::write(&plain_path, fasta.as_bytes()).expect("write plain FASTA");
    write_gzip(&gzip_path, fasta.as_bytes());

    let plain_mmap_dir: PathBuf = workspace.path("plain-mmap");
    let plain_stream_dir: PathBuf = workspace.path("plain-stream");
    let gzip_mmap_dir: PathBuf = workspace.path("gzip-mmap");
    let gzip_stream_dir: PathBuf = workspace.path("gzip-stream");

    for (input, output_dir, mode) in [
        (&plain_path, &plain_mmap_dir, "mmap"),
        (&plain_path, &plain_stream_dir, "stream"),
        (&gzip_path, &gzip_mmap_dir, "mmap"),
        (&gzip_path, &gzip_stream_dir, "stream"),
    ] {
        let output: Output = run_fasta_revcomp_scan(input, output_dir, mode);
        assert_success(&output);
    }

    let plain_mmap_without_revcomp: PathBuf = workspace.path("plain-mmap-forward-only");
    let plain_stream_without_revcomp: PathBuf = workspace.path("plain-stream-forward-only");
    assert_success(&run_fasta_scan(
        &plain_path,
        &plain_mmap_without_revcomp,
        "mmap",
    ));
    assert_success(&run_fasta_scan(
        &plain_path,
        &plain_stream_without_revcomp,
        "stream",
    ));

    let plain_mmap_snapshot: Vec<(String, Vec<u8>)> = directory_snapshot(&plain_mmap_dir);
    let plain_stream_snapshot: Vec<(String, Vec<u8>)> = directory_snapshot(&plain_stream_dir);
    let gzip_mmap_snapshot: Vec<(String, Vec<u8>)> = directory_snapshot(&gzip_mmap_dir);
    let gzip_stream_snapshot: Vec<(String, Vec<u8>)> = directory_snapshot(&gzip_stream_dir);
    let expected_revcomp: Vec<(String, Vec<u8>)> =
        snapshot_matching(&plain_mmap_snapshot, |name| name.contains(".revcomp."));
    assert_eq!(
        expected_revcomp,
        snapshot_matching(&plain_stream_snapshot, |name| name.contains(".revcomp."))
    );
    assert_eq!(
        expected_revcomp,
        snapshot_matching(&gzip_mmap_snapshot, |name| name.contains(".revcomp."))
    );
    assert_eq!(
        expected_revcomp,
        snapshot_matching(&gzip_stream_snapshot, |name| name.contains(".revcomp."))
    );

    let mmap_forward_with_revcomp: Vec<(String, Vec<u8>)> =
        snapshot_matching(&plain_mmap_snapshot, |name| !name.contains(".revcomp."));
    let stream_forward_with_revcomp: Vec<(String, Vec<u8>)> =
        snapshot_matching(&plain_stream_snapshot, |name| !name.contains(".revcomp."));
    assert_eq!(
        mmap_forward_with_revcomp,
        directory_snapshot(&plain_mmap_without_revcomp),
        "enabling revcomp changed mmap forward outputs"
    );
    assert_eq!(
        stream_forward_with_revcomp,
        directory_snapshot(&plain_stream_without_revcomp),
        "enabling revcomp changed stream forward outputs"
    );

    let names: Vec<&str> = plain_mmap_snapshot
        .iter()
        .map(|(name, _bytes)| name.as_str())
        .collect();
    for chromosome in ["negative_only", "both_strands", "block_boundary"] {
        for suffix in [
            "g4.csv",
            "g4.overlap.csv",
            "g4.family.csv",
            "g4.revcomp.csv",
            "g4.revcomp.overlap.csv",
            "g4.revcomp.family.csv",
        ] {
            let expected_name: String = format!("{chromosome}.{suffix}");
            assert!(
                names.contains(&expected_name.as_str()),
                "missing output {expected_name}; observed {names:?}"
            );
        }
    }

    let block_primary: &Vec<u8> = &plain_mmap_snapshot
        .iter()
        .find(|(name, _bytes)| name == "block_boundary.g4.revcomp.csv")
        .expect("find block-boundary revcomp output")
        .1;
    let block_overlap: &Vec<u8> = &plain_mmap_snapshot
        .iter()
        .find(|(name, _bytes)| name == "block_boundary.g4.revcomp.overlap.csv")
        .expect("find block-boundary revcomp overlap output")
        .1;
    assert!(csv_row_count(block_primary) > 0);
    assert!(
        csv_row_count(block_overlap) > csv_row_count(block_primary),
        "long C-rich family must exercise overlapping reverse-complement hits"
    );
}

#[test]
fn circular_revcomp_keeps_expanded_reference_coordinates() {
    let workspace: TestWorkspace = TestWorkspace::new("circular");
    let forward_path: PathBuf = workspace.path("result.csv");
    let reverse_path: PathBuf = workspace.path("result.revcomp.csv");
    let sequence: &str = "CCCCCCCTCCCCTCCCCTC";

    let output: Output = run(qgrs_command().args([
        "--sequence",
        sequence,
        "--min-tetrads",
        "4",
        "--circular",
        "--output",
        forward_path
            .to_str()
            .expect("output path must be valid UTF-8"),
        "--revcomp",
    ]));
    assert_success(&output);

    let reverse_csv: String =
        fs::read_to_string(&reverse_path).expect("read circular revcomp CSV output");
    let row: Vec<&str> = reverse_csv
        .lines()
        .skip(1)
        .map(|line| line.split(',').collect::<Vec<&str>>())
        .find(|columns| columns.get(8) == Some(&"GGGGAGGGGAGGGGAGGGG"))
        .expect("find expected cross-origin revcomp result");
    assert_eq!(row[0], "4");
    assert_eq!(row[1], "22");
    assert_eq!(row[2], "19");
    assert!(
        row[1].parse::<usize>().expect("parse expanded end") > sequence.len(),
        "circular revcomp result must retain an expanded end coordinate"
    );

    let fasta_path: PathBuf = workspace.path("circular.fa");
    fs::write(&fasta_path, format!(">circular\n{sequence}\n")).expect("write circular FASTA");
    let mmap_dir: PathBuf = workspace.path("circular-mmap");
    let stream_dir: PathBuf = workspace.path("circular-stream");
    assert_success(&run_circular_fasta_revcomp_scan(
        &fasta_path,
        &mmap_dir,
        "mmap",
    ));
    assert_success(&run_circular_fasta_revcomp_scan(
        &fasta_path,
        &stream_dir,
        "stream",
    ));
    let mmap_revcomp: Vec<(String, Vec<u8>)> =
        snapshot_matching(&directory_snapshot(&mmap_dir), |name| {
            name.contains(".revcomp.")
        });
    let stream_revcomp: Vec<(String, Vec<u8>)> =
        snapshot_matching(&directory_snapshot(&stream_dir), |name| {
            name.contains(".revcomp.")
        });
    assert_eq!(mmap_revcomp, stream_revcomp);
}

#[test]
fn stream_revcomp_cleans_tmpdir_on_success_and_invalid_iupac_error() {
    let workspace: TestWorkspace = TestWorkspace::new("tmp-cleanup");
    let spool_dir: PathBuf = workspace.path("tmp");
    fs::create_dir(&spool_dir).expect("create stream temporary directory");

    let valid_fasta: PathBuf = workspace.path("valid.fa");
    fs::write(
        &valid_fasta,
        format!(">valid_chr\nACGTRYSWKMBDHVNU{ASYMMETRIC_G4_SOURCE}\n"),
    )
    .expect("write valid FASTA");
    let valid_output_dir: PathBuf = workspace.path("valid-output");
    let valid_output: Output = run(qgrs_command().env("TMPDIR", &spool_dir).args([
        "--file",
        valid_fasta
            .to_str()
            .expect("valid FASTA path must be valid UTF-8"),
        "--mode",
        "stream",
        "--output-dir",
        valid_output_dir
            .to_str()
            .expect("valid output path must be valid UTF-8"),
        "--min-tetrads",
        "4",
        "--revcomp",
    ]));
    assert_success(&valid_output);
    assert_directory_empty(&spool_dir);

    let invalid_fasta: PathBuf = workspace.path("invalid.fa");
    fs::write(
        &invalid_fasta,
        b">bad_chr\nAZACGTRYSWKMBDHVNUGGGGAGGGGAGGGG\n",
    )
    .expect("write invalid FASTA");
    let invalid_output_dir: PathBuf = workspace.path("invalid-output");
    let invalid_output: Output = run(qgrs_command().env("TMPDIR", &spool_dir).args([
        "--file",
        invalid_fasta
            .to_str()
            .expect("invalid FASTA path must be valid UTF-8"),
        "--mode",
        "stream",
        "--output-dir",
        invalid_output_dir
            .to_str()
            .expect("invalid output path must be valid UTF-8"),
        "--min-tetrads",
        "4",
        "--revcomp",
    ]));
    assert_failure(&invalid_output);

    let stderr: String = String::from_utf8_lossy(&invalid_output.stderr).into_owned();
    assert!(stderr.contains("bad_chr"), "missing chromosome: {stderr}");
    assert!(
        stderr.contains("position 2"),
        "missing 1-based position: {stderr}"
    );
    assert!(
        stderr.to_ascii_lowercase().contains('z'),
        "missing invalid character: {stderr}"
    );
    assert_directory_empty(&spool_dir);
}
