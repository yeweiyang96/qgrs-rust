use std::collections::HashMap;
use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use qgrs_rust::qgrs::{
    self, DEFAULT_MAX_G4_LENGTH, DEFAULT_MAX_RUN, G4, InputMode, QuartetBase, ScanLimits,
    SequenceTopology,
};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

fn main() {
    // Initialize Rayon global thread pool to match machine CPU count.
    // This makes parallelism deterministic across runs and avoids relying on
    // the environment variable `RAYON_NUM_THREADS`.
    let threads = num_cpus::get();
    let _ = ThreadPoolBuilder::new().num_threads(threads).build_global();

    if let Err(err) = run_env(env::args().skip(1)) {
        eprintln!("Error: {err}");
        std::process::exit(1);
    }
}

fn run_env<I>(mut args: I) -> Result<(), String>
where
    I: Iterator<Item = String>,
{
    let mut sequence_arg: Option<String> = None;
    let mut file_arg: Option<PathBuf> = None;
    let mut min_tetrads: usize = 2;
    let mut min_score: i32 = 17;
    let mut max_run: usize = DEFAULT_MAX_RUN;
    let mut max_g4_length: usize = DEFAULT_MAX_G4_LENGTH;
    let mut format = OutputFormat::Csv;
    let mut output_path: Option<PathBuf> = None;
    let mut output_dir: Option<PathBuf> = None;
    let mut mode = InputMode::Mmap;
    let mut include_overlap = false;
    let mut include_revcomp = false;
    let mut circular = false;
    let mut target_base = QuartetBase::G;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--sequence" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --sequence"))?;
                sequence_arg = Some(value);
            }
            "--file" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --file"))?;
                file_arg = Some(PathBuf::from(value));
            }
            "--min-tetrads" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --min-tetrads"))?
                    .parse::<usize>()
                    .map_err(|_| usage("--min-tetrads must be a positive integer"))?;
                if value == 0 {
                    return Err(usage("--min-tetrads must be > 0"));
                }
                min_tetrads = value;
            }
            "--min-score" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --min-score"))?
                    .parse::<i32>()
                    .map_err(|_| usage("--min-score must be an integer"))?;
                min_score = value;
            }
            "--format" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --format"))?;
                format = value.try_into()?;
            }
            "--mode" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --mode"))?;
                mode = parse_mode(&value)?;
            }
            "--base" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --base"))?;
                target_base = parse_base(&value)?;
            }
            "--max-run" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --max-run"))?
                    .parse::<usize>()
                    .map_err(|_| usage("--max-run must be a positive integer"))?;
                if value == 0 {
                    return Err(usage("--max-run must be > 0"));
                }
                max_run = value;
            }
            "--max-g-run" => {
                return Err(usage("--max-g-run was replaced by --max-run"));
            }
            "--max-g4-length" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --max-g4-length"))?
                    .parse::<usize>()
                    .map_err(|_| usage("--max-g4-length must be a positive integer"))?;
                if value == 0 {
                    return Err(usage("--max-g4-length must be > 0"));
                }
                max_g4_length = value;
            }
            "--output" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --output"))?;
                output_path = Some(PathBuf::from(value));
            }
            "--output-dir" => {
                let value = args
                    .next()
                    .ok_or_else(|| usage("missing value for --output-dir"))?;
                output_dir = Some(PathBuf::from(value));
            }
            "--overlap" => {
                include_overlap = true;
            }
            "--revcomp" => {
                include_revcomp = true;
            }
            "--circular" => {
                circular = true;
            }
            "--help" | "-h" => return Err(usage("")),
            other => {
                return Err(usage(&format!("unknown argument '{other}'")));
            }
        }
    }

    let input = match (sequence_arg, file_arg) {
        (Some(_), Some(_)) => {
            return Err(usage("cannot provide both --sequence and --file"));
        }
        (Some(seq), None) => InputSpec::Inline(seq),
        (None, Some(path)) => InputSpec::File(path),
        (None, None) => return Err(usage("must provide --sequence or --file")),
    };

    let min_required_length = min_tetrads
        .checked_mul(4)
        .ok_or_else(|| usage("--min-tetrads is too large"))?;
    if max_run < min_tetrads {
        return Err(usage("--max-run must be ≥ --min-tetrads"));
    }
    if max_g4_length < min_required_length {
        return Err(usage("--max-g4-length must be ≥ 4 * --min-tetrads"));
    }

    let limits = ScanLimits::new(max_g4_length, max_run);
    let topology = if circular {
        SequenceTopology::Circular
    } else {
        SequenceTopology::Linear
    };
    let scan = ScanConfig::new(min_tetrads, min_score, limits, topology, target_base);

    match input {
        InputSpec::Inline(seq) => {
            if output_dir.is_some() {
                return Err(usage("--output-dir can only be used with --file"));
            }
            process_inline_sequence(
                seq,
                format,
                output_path,
                scan,
                include_overlap,
                include_revcomp,
            )?;
        }
        InputSpec::File(path) => {
            if output_path.is_some() {
                return Err(usage(
                    "--output is only valid with --sequence; use --output-dir for --file",
                ));
            }
            process_fasta_file(
                path,
                mode,
                format,
                scan,
                output_dir,
                include_overlap,
                include_revcomp,
            )?;
        }
    }
    Ok(())
}

fn usage(reason: &str) -> String {
    let mut msg = String::new();
    if !reason.is_empty() {
        msg.push_str(reason);
        msg.push('\n');
    }
    msg.push_str("Usage: cargo run --bin qgrs -- [--sequence <SEQ> | --file <PATH>] [options]\n");
    msg.push_str("Options:\n");
    msg.push_str("  --sequence <SEQ>     Inline DNA/RNA sequence to scan\n");
    msg.push_str(
        "  --file <PATH>        Read sequences from FASTA/FASTA.gz (chromosomes split independently)\n",
    );
    msg.push_str("  --min-tetrads <N>    Minimum tetrads to seed (default 2)\n");
    msg.push_str("  --min-score <S>      Minimum score (default 17)\n");
    msg.push_str(
        "  --base <g|c>         Tetrad base to scan: g for G4, c for i-motif (default g)\n",
    );
    msg.push_str("  --max-run <N>        Maximum allowed target-base run length (default 10)\n");
    msg.push_str("  --max-g4-length <N>  Maximum allowed G4 length in bp (default 45)\n");
    msg.push_str("  --format <csv|parquet>  Output format (default csv)\n");
    msg.push_str(
        "  --output <PATH>     Destination file when using --sequence (required for parquet)\n",
    );
    msg.push_str("  --output-dir <DIR>  Directory for per-chromosome exports when using --file\n");
    msg.push_str("  --mode <mmap|stream> Input mode when using --file (default mmap)\n");
    msg.push_str(
        "  --overlap            Emit raw hits (.overlap.<format>) and family ranges (.family.<format>)\n",
    );
    msg.push_str(
        "  --revcomp            Also scan the reverse-complement strand into .revcomp.<format>\n",
    );
    msg.push_str("  --circular           Treat each sequence/chromosome as circular\n");
    msg.push_str("  --help               Show this message\n");
    msg
}

fn parse_mode(value: &str) -> Result<InputMode, String> {
    match value {
        "mmap" => Ok(InputMode::Mmap),
        "stream" => Ok(InputMode::Stream),
        _ => Err(usage("--mode must be either 'mmap' or 'stream'")),
    }
}

fn parse_base(value: &str) -> Result<QuartetBase, String> {
    if value.len() != 1 {
        return Err(usage("--base must be exactly one character: g or c"));
    }
    match value.as_bytes()[0].to_ascii_lowercase() {
        b'g' => Ok(QuartetBase::G),
        b'c' => Ok(QuartetBase::C),
        _ => Err(usage("--base must be either 'g' for G4 or 'c' for i-motif")),
    }
}

enum InputSpec {
    Inline(String),
    File(PathBuf),
}

#[derive(Clone, Copy)]
struct ScanConfig {
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
}

impl ScanConfig {
    fn new(
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        target_base: QuartetBase,
    ) -> Self {
        Self {
            min_tetrads,
            min_score,
            limits,
            topology,
            target_base,
        }
    }

    fn min_tetrads(self) -> usize {
        self.min_tetrads
    }

    fn min_score(self) -> i32 {
        self.min_score
    }

    fn limits(self) -> ScanLimits {
        self.limits
    }

    fn topology(self) -> SequenceTopology {
        self.topology
    }

    fn target_base(self) -> QuartetBase {
        self.target_base
    }
}

fn process_inline_sequence(
    sequence: String,
    format: OutputFormat,
    output_path: Option<PathBuf>,
    scan: ScanConfig,
    include_overlap: bool,
    include_revcomp: bool,
) -> Result<(), String> {
    let mut normalized = sequence.into_bytes();
    normalized.make_ascii_lowercase();
    let sequence_len = normalized.len();
    if include_overlap && output_path.is_none() {
        return Err(usage("--overlap requires --output when using --sequence"));
    }
    let revcomp_output_is_missing =
        output_path.as_deref().is_none() || output_path.as_deref() == Some(Path::new("-"));
    if include_revcomp && revcomp_output_is_missing {
        return Err(usage(
            "--revcomp requires a file path via --output when using --sequence",
        ));
    }

    let (results, family_ranges, raw_hits) = run_scan_for_export(
        Arc::new(normalized.clone()),
        scan,
        include_overlap,
        sequence_len,
    );
    let revcomp_results = if include_revcomp {
        Some(run_revcomp_scan_for_export(
            &normalized,
            "inline sequence",
            scan,
            include_overlap,
        )?)
    } else {
        None
    };
    write_primary_output(
        output_path.as_deref(),
        format,
        &results,
        scan.topology(),
        sequence_len,
    )?;

    if include_overlap {
        let base = output_path
            .as_ref()
            .expect("overlap outputs require an explicit --output path");
        write_overlap_exports(
            base,
            format,
            raw_hits.as_ref().unwrap(),
            &family_ranges,
            scan.topology(),
            sequence_len,
        )?;
    }
    if let Some((hits, ranges, raw_hits)) = revcomp_results {
        let base = output_path
            .as_ref()
            .expect("revcomp output requires an explicit --output path");
        write_revcomp_exports(
            base,
            format,
            &hits,
            &ranges,
            raw_hits.as_deref(),
            scan.topology(),
            sequence_len,
        )?;
    }

    Ok(())
}

fn process_fasta_file(
    path: PathBuf,
    mode: InputMode,
    format: OutputFormat,
    scan: ScanConfig,
    output_dir: Option<PathBuf>,
    include_overlap: bool,
    include_revcomp: bool,
) -> Result<(), String> {
    let dir = output_dir.ok_or_else(|| usage("--output-dir is required when --file is used"))?;
    fs::create_dir_all(&dir).map_err(|err| format!("failed to create {dir:?}: {err}"))?;
    let mut name_counts: HashMap<String, usize> = HashMap::new();
    match mode {
        InputMode::Mmap => {
            let sequences = qgrs::load_sequences_from_path(&path, InputMode::Mmap)
                .map_err(|err| format!("failed to read {path:?}: {err}"))?;
            if sequences.is_empty() {
                return Err(format!("no sequences found in {path:?}"));
            }
            let mut chrom_outputs = Vec::with_capacity(sequences.len());
            for chrom in sequences {
                let filename = next_output_filename(
                    chrom.name(),
                    format,
                    scan.target_base(),
                    &mut name_counts,
                );
                chrom_outputs.push((chrom, dir.join(filename)));
            }
            chrom_outputs.into_par_iter().try_for_each(
                |(chrom, filepath)| -> Result<(), String> {
                    let (name, sequence) = chrom.into_parts();
                    let sequence_len = sequence.len();
                    let (results, family_ranges, raw_hits) =
                        run_scan_for_export(sequence.clone(), scan, include_overlap, sequence_len);
                    let revcomp_results = if include_revcomp {
                        Some(run_revcomp_scan_for_export(
                            sequence.as_slice(),
                            &name,
                            scan,
                            include_overlap,
                        )?)
                    } else {
                        None
                    };
                    write_results_to_path(
                        &filepath,
                        format,
                        &results,
                        scan.topology(),
                        sequence_len,
                    )?;
                    if include_overlap {
                        let raw_hits = raw_hits
                            .as_ref()
                            .expect("raw hits must be captured when overlap is requested");
                        write_overlap_exports(
                            &filepath,
                            format,
                            raw_hits,
                            &family_ranges,
                            scan.topology(),
                            sequence_len,
                        )?;
                    }
                    if let Some((hits, ranges, raw_hits)) = revcomp_results {
                        write_revcomp_exports(
                            &filepath,
                            format,
                            &hits,
                            &ranges,
                            raw_hits.as_deref(),
                            scan.topology(),
                            sequence_len,
                        )?;
                    }
                    Ok(())
                },
            )?;
        }
        InputMode::Stream => {
            let mut processed = 0usize;
            if include_revcomp && include_overlap {
                qgrs::stream::process_fasta_stream_bidirectional_with_limits_overlap_topology_and_len_with_base(
                    &path,
                    scan.min_tetrads(),
                    scan.min_score(),
                    scan.limits(),
                    scan.topology(),
                    scan.target_base(),
                    |name, mut forward, mut reverse_rc, sequence_len| {
                        processed += 1;
                        let filename = next_output_filename(
                            &name,
                            format,
                            scan.target_base(),
                            &mut name_counts,
                        );
                        let filepath = dir.join(&filename);
                        write_results_to_path(
                            &filepath,
                            format,
                            &forward.hits,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        let forward_raw_hits = forward
                            .raw_hits
                            .take()
                            .expect("raw hits missing from forward bidirectional stream results");
                        write_overlap_exports(
                            &filepath,
                            format,
                            &forward_raw_hits,
                            &forward.family_ranges,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;

                        let reverse_raw_hits = reverse_rc
                            .raw_hits
                            .take()
                            .expect("raw hits missing from reverse bidirectional stream results");
                        write_revcomp_exports(
                            &filepath,
                            format,
                            &reverse_rc.hits,
                            &reverse_rc.family_ranges,
                            Some(&reverse_raw_hits),
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        Ok(())
                    },
                )
                .map_err(|err| format!("failed to process {path:?}: {err}"))?;
            } else if include_revcomp {
                qgrs::stream::process_fasta_stream_bidirectional_with_limits_topology_and_len_with_base(
                    &path,
                    scan.min_tetrads(),
                    scan.min_score(),
                    scan.limits(),
                    scan.topology(),
                    scan.target_base(),
                    |name, forward_hits, reverse_hits_rc, sequence_len| {
                        processed += 1;
                        let filename = next_output_filename(
                            &name,
                            format,
                            scan.target_base(),
                            &mut name_counts,
                        );
                        let filepath = dir.join(&filename);
                        write_results_to_path(
                            &filepath,
                            format,
                            &forward_hits,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        write_revcomp_exports(
                            &filepath,
                            format,
                            &reverse_hits_rc,
                            &[],
                            None,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        Ok(())
                    },
                )
                .map_err(|err| format!("failed to process {path:?}: {err}"))?;
            } else if include_overlap {
                qgrs::stream::process_fasta_stream_with_limits_overlap_topology_and_len_with_base(
                    &path,
                    scan.min_tetrads(),
                    scan.min_score(),
                    scan.limits(),
                    scan.topology(),
                    scan.target_base(),
                    |name, mut stream_results, sequence_len| {
                        processed += 1;
                        let filename = next_output_filename(
                            &name,
                            format,
                            scan.target_base(),
                            &mut name_counts,
                        );
                        let filepath = dir.join(&filename);
                        write_results_to_path(
                            &filepath,
                            format,
                            &stream_results.hits,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        let raw_hits = stream_results
                            .raw_hits
                            .take()
                            .expect("raw hits missing from overlap stream results");

                        write_overlap_exports(
                            &filepath,
                            format,
                            &raw_hits,
                            &stream_results.family_ranges,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        Ok(())
                    },
                )
                .map_err(|err| format!("failed to process {path:?}: {err}"))?;
            } else {
                qgrs::stream::process_fasta_stream_with_limits_topology_and_len_with_base(
                    &path,
                    scan.min_tetrads(),
                    scan.min_score(),
                    scan.limits(),
                    scan.topology(),
                    scan.target_base(),
                    |name, results, sequence_len| {
                        processed += 1;
                        let filename = next_output_filename(
                            &name,
                            format,
                            scan.target_base(),
                            &mut name_counts,
                        );
                        let filepath = dir.join(&filename);
                        write_results_to_path(
                            &filepath,
                            format,
                            &results,
                            scan.topology(),
                            sequence_len,
                        )
                        .map_err(io::Error::other)?;
                        Ok(())
                    },
                )
                .map_err(|err| format!("failed to process {path:?}: {err}"))?;
            }
            if processed == 0 {
                return Err(format!("no sequences found in {path:?}"));
            }
        }
    }
    Ok(())
}

fn next_output_filename(
    name: &str,
    format: OutputFormat,
    target_base: QuartetBase,
    counts: &mut HashMap<String, usize>,
) -> String {
    let sanitized = sanitize_name(name);
    // 处理同名染色体的重复输出(万一)
    let entry = counts.entry(sanitized.clone()).or_insert(0);
    let suffix = if *entry == 0 {
        String::new()
    } else {
        format!("_{}", entry)
    };
    *entry += 1;
    format!(
        "{}{suffix}.{}.{}",
        sanitized,
        output_motif_label(target_base),
        format.extension()
    )
}

fn output_motif_label(target_base: QuartetBase) -> &'static str {
    match target_base {
        QuartetBase::G => "g4",
        QuartetBase::C => "i-motif",
    }
}

fn sanitize_name(raw: &str) -> String {
    // let mut sanitized = String::new();
    // for ch in raw.chars() {
    //     if ch.is_ascii_alphanumeric() || matches!(ch, '-' | '_') {
    //         sanitized.push(ch);
    //     } else {
    //         sanitized.push('_');
    //     }
    // }
    // if sanitized.is_empty() {
    //     "chromosome".to_string()
    // } else {
    //     sanitized
    // }
    if raw.is_empty() {
        "chromosome".to_string()
    } else {
        raw.to_string()
    }
}

type ConsolidatedResults = (Vec<G4>, Vec<(usize, usize)>, Option<Vec<G4>>);

fn consolidate_for_export(
    raw: Vec<G4>,
    capture_raw: bool,
    topology: SequenceTopology,
    sequence_len: usize,
) -> ConsolidatedResults {
    if capture_raw {
        let raw_copy = raw.clone();
        let (hits, ranges) = qgrs::consolidate_g4s_with_topology(raw, topology, sequence_len);
        (hits, ranges, Some(raw_copy))
    } else {
        let (hits, ranges) = qgrs::consolidate_g4s_with_topology(raw, topology, sequence_len);
        (hits, ranges, None)
    }
}

fn run_scan_for_export(
    sequence: Arc<Vec<u8>>,
    scan: ScanConfig,
    capture_raw: bool,
    sequence_len: usize,
) -> ConsolidatedResults {
    let raw = qgrs::find_owned_bytes_with_topology_and_base(
        sequence,
        scan.min_tetrads(),
        scan.min_score(),
        scan.limits(),
        scan.topology(),
        scan.target_base(),
    );
    consolidate_for_export(raw, capture_raw, scan.topology(), sequence_len)
}

fn run_revcomp_scan_for_export(
    forward_sequence: &[u8],
    sequence_name: &str,
    scan: ScanConfig,
    capture_raw: bool,
) -> Result<ConsolidatedResults, String> {
    let sequence_len = forward_sequence.len();
    let revcomp = reverse_complement_lowercase(forward_sequence, sequence_name)?;
    Ok(run_scan_for_export(
        Arc::new(revcomp),
        scan,
        capture_raw,
        sequence_len,
    ))
}

fn reverse_complement_lowercase(sequence: &[u8], sequence_name: &str) -> Result<Vec<u8>, String> {
    let mut revcomp = Vec::with_capacity(sequence.len());
    for (reverse_index, byte) in sequence.iter().rev().copied().enumerate() {
        let original_position = sequence.len() - reverse_index;
        revcomp.push(complement_iupac_lowercase(
            byte,
            sequence_name,
            original_position,
        )?);
    }
    Ok(revcomp)
}

fn complement_iupac_lowercase(
    byte: u8,
    sequence_name: &str,
    original_position: usize,
) -> Result<u8, String> {
    let complement = match byte.to_ascii_lowercase() {
        b'a' => b't',
        b't' | b'u' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'r' => b'y',
        b'y' => b'r',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        b's' => b's',
        b'w' => b'w',
        b'n' => b'n',
        _ => {
            return Err(format!(
                "invalid nucleotide byte 0x{byte:02X} ('{}') in {sequence_name} at 1-based position {original_position}",
                char::from(byte).escape_default()
            ));
        }
    };
    Ok(complement)
}

fn project_revcomp_hits_to_forward(
    hits: &[G4],
    topology: SequenceTopology,
    sequence_len: usize,
) -> Vec<G4> {
    let mut projected = hits.to_vec();
    for hit in &mut projected {
        let (start, end) = map_revcomp_interval(hit.start, hit.end, topology, sequence_len);
        hit.start = start;
        hit.end = end;
    }
    projected.sort_by(|left, right| (left.start, left.end).cmp(&(right.start, right.end)));
    projected
}

fn project_revcomp_ranges_to_forward(
    ranges: &[(usize, usize)],
    topology: SequenceTopology,
    sequence_len: usize,
) -> Vec<(usize, usize)> {
    let mut projected = ranges
        .iter()
        .map(|(start, end)| map_revcomp_interval(*start, *end, topology, sequence_len))
        .collect::<Vec<_>>();
    projected.sort_unstable();
    projected
}

fn map_revcomp_interval(
    start_rc: usize,
    end_rc: usize,
    topology: SequenceTopology,
    sequence_len: usize,
) -> (usize, usize) {
    if sequence_len == 0 {
        return (start_rc, end_rc);
    }
    let end_anchor = if topology.is_circular() {
        ((end_rc - 1) % sequence_len) + 1
    } else {
        end_rc
    };
    let start = sequence_len - end_anchor + 1;
    let end = start + (end_rc - start_rc);
    (start, end)
}

fn write_primary_output(
    output_path: Option<&Path>,
    format: OutputFormat,
    results: &[G4],
    _topology: SequenceTopology,
    _sequence_len: usize,
) -> Result<(), String> {
    match format {
        OutputFormat::Csv => {
            let csv = qgrs::render_csv_results(results);
            if let Some(path) = output_path {
                fs::write(path, csv).map_err(|err| format!("failed to write {path:?}: {err}"))?;
            } else {
                print!("{csv}");
            }
            Ok(())
        }
        OutputFormat::Parquet => {
            let path =
                output_path.ok_or_else(|| usage("--output is required when --format parquet"))?;
            write_results_to_path(path, format, results, _topology, _sequence_len)
        }
    }
}

fn write_results_to_path(
    path: &Path,
    format: OutputFormat,
    results: &[G4],
    _topology: SequenceTopology,
    _sequence_len: usize,
) -> Result<(), String> {
    match format {
        OutputFormat::Csv => {
            let csv = qgrs::render_csv_results(results);
            fs::write(path, csv).map_err(|err| format!("failed to write {path:?}: {err}"))?;
        }
        OutputFormat::Parquet => {
            let file = fs::File::create(path)
                .map_err(|err| format!("failed to create {path:?}: {err}"))?;
            qgrs::write_parquet_results(results, file)
                .map_err(|err| format!("failed to write parquet {path:?}: {err}"))?;
        }
    }
    Ok(())
}

fn write_overlap_exports(
    base: &Path,
    format: OutputFormat,
    raw_hits: &[G4],
    family_ranges: &[(usize, usize)],
    _topology: SequenceTopology,
    _sequence_len: usize,
) -> Result<(), String> {
    let overlap_path = overlap_path(base, format);
    let family_path = family_path(base, format);
    match format {
        OutputFormat::Csv => {
            let overlap_csv = qgrs::render_csv_results(raw_hits);
            fs::write(&overlap_path, overlap_csv)
                .map_err(|err| format!("failed to write {overlap_path:?}: {err}"))?;

            let family_csv = qgrs::render_family_ranges_csv(family_ranges);
            fs::write(&family_path, family_csv)
                .map_err(|err| format!("failed to write {family_path:?}: {err}"))?;
        }
        OutputFormat::Parquet => {
            let overlap_file = fs::File::create(&overlap_path)
                .map_err(|err| format!("failed to create {overlap_path:?}: {err}"))?;
            qgrs::write_parquet_results(raw_hits, overlap_file)
                .map_err(|err| format!("failed to write parquet {overlap_path:?}: {err}"))?;

            let family_file = fs::File::create(&family_path)
                .map_err(|err| format!("failed to create {family_path:?}: {err}"))?;
            qgrs::write_parquet_family_ranges(family_ranges, family_file)
                .map_err(|err| format!("failed to write parquet {family_path:?}: {err}"))?;
        }
    }
    Ok(())
}

fn write_revcomp_exports(
    primary_path: &Path,
    format: OutputFormat,
    hits: &[G4],
    family_ranges: &[(usize, usize)],
    raw_hits: Option<&[G4]>,
    topology: SequenceTopology,
    sequence_len: usize,
) -> Result<(), String> {
    let revcomp_path = revcomp_output_path(primary_path, format);
    let projected_hits = project_revcomp_hits_to_forward(hits, topology, sequence_len);
    write_results_to_path(
        &revcomp_path,
        format,
        &projected_hits,
        topology,
        sequence_len,
    )?;
    if let Some(raw_hits) = raw_hits {
        let projected_raw_hits = project_revcomp_hits_to_forward(raw_hits, topology, sequence_len);
        let projected_ranges =
            project_revcomp_ranges_to_forward(family_ranges, topology, sequence_len);
        write_overlap_exports(
            &revcomp_path,
            format,
            &projected_raw_hits,
            &projected_ranges,
            topology,
            sequence_len,
        )?;
    }
    Ok(())
}

fn revcomp_output_path(base: &Path, format: OutputFormat) -> PathBuf {
    append_output_suffix(base, ".revcomp", format)
}

fn overlap_path(base: &Path, format: OutputFormat) -> PathBuf {
    append_output_suffix(base, ".overlap", format)
}

fn family_path(base: &Path, format: OutputFormat) -> PathBuf {
    append_output_suffix(base, ".family", format)
}

fn append_output_suffix(path: &Path, suffix: &str, format: OutputFormat) -> PathBuf {
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let stem = path
        .file_stem()
        .and_then(|s| s.to_str())
        .filter(|s| !s.is_empty())
        .unwrap_or("chromosome");
    let name = format!("{stem}{suffix}.{}", format.extension());
    parent.join(name)
}

#[derive(Clone, Copy)]
enum OutputFormat {
    Csv,
    Parquet,
}

impl TryFrom<String> for OutputFormat {
    type Error = String;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        match value.as_str() {
            "csv" => Ok(OutputFormat::Csv),
            "parquet" => Ok(OutputFormat::Parquet),
            _ => Err(usage("--format must be either 'csv' or 'parquet'")),
        }
    }
}

impl OutputFormat {
    fn extension(&self) -> &'static str {
        match self {
            OutputFormat::Csv => "csv",
            OutputFormat::Parquet => "parquet",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

    use flate2::Compression;
    use flate2::write::GzEncoder;

    #[test]
    fn default_limits_are_valid() {
        let limits = ScanLimits::default();
        assert_eq!(limits.max_g4_length, DEFAULT_MAX_G4_LENGTH);
        assert_eq!(limits.max_run, DEFAULT_MAX_RUN);
        assert!(limits.max_run >= 2);
        assert!(limits.max_g4_length >= 8);
    }

    #[test]
    fn usage_fails_on_invalid_limits() {
        let err = run_with_args([
            "--sequence",
            "GGGG",
            "--min-tetrads",
            "4",
            "--max-g4-length",
            "12",
            "--max-run",
            "4",
        ]);
        assert!(err.is_err());
        let msg = err.unwrap_err().to_string();
        assert!(msg.contains("max-g4-length"));
    }

    #[test]
    fn overlap_requires_output_for_inline() {
        let err = run_with_args(["--sequence", "GGGG", "--overlap"]);
        assert!(err.is_err());
        let msg = err.unwrap_err();
        assert!(msg.contains("--overlap requires --output"));
    }

    #[test]
    fn old_max_g_run_is_rejected_with_migration_guidance() {
        let err = run_with_args(["--sequence", "GGGG", "--max-g-run", "4"]);
        assert!(err.is_err());
        let msg = err.unwrap_err();
        assert!(msg.contains("--max-g-run was replaced by --max-run"));
    }

    #[test]
    fn invalid_base_values_are_rejected() {
        for value in ["a", "t", "gc", ""] {
            let err = run_with_owned_args(vec![
                "--sequence".to_string(),
                "GGGG".to_string(),
                "--base".to_string(),
                value.to_string(),
            ]);
            assert!(err.is_err(), "base value {value:?} should fail");
        }
    }

    #[test]
    fn output_filename_includes_motif_label() {
        let mut counts = HashMap::new();
        assert_eq!(
            next_output_filename("chr1", OutputFormat::Parquet, QuartetBase::G, &mut counts),
            "chr1.g4.parquet"
        );
        assert_eq!(
            next_output_filename("chr1", OutputFormat::Parquet, QuartetBase::G, &mut counts),
            "chr1_1.g4.parquet"
        );

        let mut counts = HashMap::new();
        assert_eq!(
            next_output_filename("chr2", OutputFormat::Csv, QuartetBase::C, &mut counts),
            "chr2.i-motif.csv"
        );
    }

    #[test]
    fn circular_flag_is_supported_for_inline_scan() {
        let result = run_with_args([
            "--sequence",
            "GAGGGGAGGGGAGGGGGGG",
            "--min-tetrads",
            "4",
            "--min-score",
            "17",
            "--circular",
        ]);
        assert!(result.is_ok());
    }

    #[test]
    fn circular_cli_outputs_keep_expanded_coordinates() {
        let base = unique_test_path("qgrs_circular_cli");
        let output = base.with_extension("csv");
        let output_str = output.to_string_lossy().into_owned();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            "GAGGGGAGGGGAGGGGGGG".to_string(),
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--circular".to_string(),
            "--overlap".to_string(),
            "--output".to_string(),
            output_str.clone(),
        ]);
        assert!(result.is_ok());

        let csv = fs::read_to_string(&output).expect("main output");
        assert!(csv.contains("\n17,35,19,4,1,1,1,84,GGGGAGGGGAGGGGAGGGG\n"));

        let overlap =
            fs::read_to_string(overlap_path(&output, OutputFormat::Csv)).expect("overlap output");
        for line in overlap.lines().skip(1) {
            let mut cols = line.split(',');
            let start = cols.next().unwrap().parse::<usize>().unwrap();
            let end = cols.next().unwrap().parse::<usize>().unwrap();
            assert!(start <= 19);
            assert!(end >= start);
        }
        assert!(overlap.lines().skip(1).any(|line| {
            let mut cols = line.split(',');
            let _start = cols.next().unwrap().parse::<usize>().unwrap();
            let end = cols.next().unwrap().parse::<usize>().unwrap();
            end > 19
        }));

        let family =
            fs::read_to_string(family_path(&output, OutputFormat::Csv)).expect("family output");
        let family_line = family.lines().nth(1).expect("family row");
        let mut cols = family_line.split(',');
        assert_eq!(cols.next(), Some("1"));
        let start = cols.next().unwrap().parse::<usize>().unwrap();
        let end = cols.next().unwrap().parse::<usize>().unwrap();
        assert!(start <= 19);
        assert!(end > 19);
        assert!(end >= start);

        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(overlap_path(&output, OutputFormat::Csv));
        let _ = fs::remove_file(family_path(&output, OutputFormat::Csv));
    }

    #[test]
    fn parquet_overlap_and_family_follow_format() {
        let base = unique_test_path("qgrs_parquet_sidecars");
        let output = base.with_extension("parquet");
        let output_str = output.to_string_lossy().into_owned();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            "GGGGAGGGGAGGGGAGGGG".to_string(),
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--format".to_string(),
            "parquet".to_string(),
            "--overlap".to_string(),
            "--output".to_string(),
            output_str,
        ]);
        assert!(result.is_ok());

        let overlap = overlap_path(&output, OutputFormat::Parquet);
        let family = family_path(&output, OutputFormat::Parquet);
        let overlap_meta = fs::metadata(&overlap).expect("overlap parquet output");
        let family_meta = fs::metadata(&family).expect("family parquet output");
        assert!(overlap_meta.len() > 0);
        assert!(family_meta.len() > 0);

        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(overlap);
        let _ = fs::remove_file(family);
    }

    #[test]
    fn base_c_inline_outputs_i_motif_hits_on_original_sequence() {
        let base = unique_test_path("qgrs_base_c_inline");
        let output = base.with_extension("csv");
        let output_str = output.to_string_lossy().into_owned();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            "AAAAAAACCCCTCCCCTCCCCTCCCCTT".to_string(),
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--base".to_string(),
            "c".to_string(),
            "--output".to_string(),
            output_str,
        ]);
        assert!(result.is_ok());

        let csv = fs::read_to_string(&output).expect("base c output");
        assert!(csv.starts_with("start,end,length,tetrads,y1,y2,y3,score,sequence\n"));
        assert!(csv.contains("\n8,26,19,4,1,1,1,84,CCCCTCCCCTCCCCTCCCC\n"));
        assert!(!csv.contains("GGGGAGGGGAGGGGAGGGG"));

        let _ = fs::remove_file(&output);
    }

    #[test]
    fn base_c_circular_inline_outputs_expanded_coordinates() {
        let base = unique_test_path("qgrs_base_c_circular");
        let output = base.with_extension("csv");
        let output_str = output.to_string_lossy().into_owned();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            "CACCCCACCCCACCCCCCC".to_string(),
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--base".to_string(),
            "c".to_string(),
            "--circular".to_string(),
            "--output".to_string(),
            output_str,
        ]);
        assert!(result.is_ok());

        let csv = fs::read_to_string(&output).expect("base c circular output");
        assert!(csv.contains("\n17,35,19,4,1,1,1,84,CCCCACCCCACCCCACCCC\n"));

        let _ = fs::remove_file(&output);
    }

    #[test]
    fn base_c_with_overlap_outputs_primary_sidecars() {
        let base = unique_test_path("qgrs_base_c_overlap");
        let output = base.with_extension("csv");
        let output_str = output.to_string_lossy().into_owned();
        let result = run_with_owned_args(vec![
            "--sequence".to_string(),
            "AAAAAAACCCCTCCCCTCCCCTCCCCTT".to_string(),
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--base".to_string(),
            "C".to_string(),
            "--overlap".to_string(),
            "--output".to_string(),
            output_str,
        ]);
        assert!(result.is_ok());

        let overlap_path = overlap_path(&output, OutputFormat::Csv);
        let family_path = family_path(&output, OutputFormat::Csv);
        assert!(fs::metadata(&overlap_path).is_ok());
        assert!(fs::metadata(&family_path).is_ok());
        let overlap = fs::read_to_string(&overlap_path).expect("overlap output");
        assert!(overlap.contains("CCCCTCCCCTCCCCTCCCC"));
        assert!(!overlap.contains("GGGGAGGGGAGGGGAGGGG"));

        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(overlap_path);
        let _ = fs::remove_file(family_path);
    }

    #[test]
    fn circular_file_outputs_match_between_mmap_and_stream() {
        let fasta = unique_test_path("qgrs_circular_modes").with_extension("fa");
        fs::write(
            &fasta,
            b">chr1\nGAGGGGAGGGGAGGGGGGG\n>chr2\nGGGCGGGGAGGGGAGGGGAG\n",
        )
        .unwrap();
        let mmap_dir = unique_test_path("qgrs_mmap_out");
        let stream_dir = unique_test_path("qgrs_stream_out");
        fs::create_dir_all(&mmap_dir).unwrap();
        fs::create_dir_all(&stream_dir).unwrap();

        let fasta_str = fasta.to_string_lossy().into_owned();
        let mmap_dir_str = mmap_dir.to_string_lossy().into_owned();
        let stream_dir_str = stream_dir.to_string_lossy().into_owned();

        let mmap_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_str.clone(),
            "--mode".to_string(),
            "mmap".to_string(),
            "--output-dir".to_string(),
            mmap_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--circular".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(mmap_result.is_ok());

        let stream_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_str,
            "--mode".to_string(),
            "stream".to_string(),
            "--output-dir".to_string(),
            stream_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--circular".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(stream_result.is_ok());

        for filename in [
            "chr1.g4.csv",
            "chr1.g4.overlap.csv",
            "chr1.g4.family.csv",
            "chr2.g4.csv",
            "chr2.g4.overlap.csv",
            "chr2.g4.family.csv",
        ] {
            let mmap_contents = fs::read_to_string(mmap_dir.join(filename)).unwrap();
            let stream_contents = fs::read_to_string(stream_dir.join(filename)).unwrap();
            assert_eq!(mmap_contents, stream_contents, "mismatch for {filename}");
        }

        let _ = fs::remove_file(&fasta);
        let _ = fs::remove_dir_all(&mmap_dir);
        let _ = fs::remove_dir_all(&stream_dir);
    }

    #[test]
    fn base_c_file_outputs_match_between_mmap_and_stream() {
        let fasta = unique_test_path("qgrs_base_c_modes").with_extension("fa");
        fs::write(
            &fasta,
            b">chr1\nAAACCCTCCCCTCCCCTCCCCAAA\n>chr2\nTTTTCCCCCTCCCCTCCCCTCCCCGG\n",
        )
        .unwrap();
        let mmap_dir = unique_test_path("qgrs_base_c_mmap_out");
        let stream_dir = unique_test_path("qgrs_base_c_stream_out");
        fs::create_dir_all(&mmap_dir).unwrap();
        fs::create_dir_all(&stream_dir).unwrap();

        let fasta_str = fasta.to_string_lossy().into_owned();
        let mmap_dir_str = mmap_dir.to_string_lossy().into_owned();
        let stream_dir_str = stream_dir.to_string_lossy().into_owned();

        let mmap_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_str.clone(),
            "--mode".to_string(),
            "mmap".to_string(),
            "--output-dir".to_string(),
            mmap_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--base".to_string(),
            "c".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(mmap_result.is_ok());

        let stream_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_str,
            "--mode".to_string(),
            "stream".to_string(),
            "--output-dir".to_string(),
            stream_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--base".to_string(),
            "c".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(stream_result.is_ok());

        for filename in [
            "chr1.i-motif.csv",
            "chr1.i-motif.overlap.csv",
            "chr1.i-motif.family.csv",
            "chr2.i-motif.csv",
            "chr2.i-motif.overlap.csv",
            "chr2.i-motif.family.csv",
        ] {
            let mmap_contents = fs::read_to_string(mmap_dir.join(filename)).unwrap();
            let stream_contents = fs::read_to_string(stream_dir.join(filename)).unwrap();
            assert_eq!(mmap_contents, stream_contents, "mismatch for {filename}");
        }
        let chr2 = fs::read_to_string(mmap_dir.join("chr2.i-motif.csv")).unwrap();
        assert!(chr2.contains("CCCCTCCCCTCCCCTCCCC"));

        let _ = fs::remove_file(&fasta);
        let _ = fs::remove_dir_all(&mmap_dir);
        let _ = fs::remove_dir_all(&stream_dir);
    }

    #[test]
    fn gzip_file_outputs_match_between_default_mmap_and_stream() {
        let base = unique_test_path("qgrs_gzip_modes");
        let fasta = base.with_extension("fa");
        let fasta_gz = base.with_extension("fna.data");
        let fasta_bytes = b">chr1\nGAGGGGAGGGGAGGGGGGG\n>chr2\nGGGCGGGGAGGGGAGGGGAG\n";
        fs::write(&fasta, fasta_bytes).unwrap();
        write_gzip(&fasta_gz, fasta_bytes);

        let mmap_dir = unique_test_path("qgrs_gzip_mmap_out");
        let stream_dir = unique_test_path("qgrs_gzip_stream_out");
        fs::create_dir_all(&mmap_dir).unwrap();
        fs::create_dir_all(&stream_dir).unwrap();

        let fasta_gz_str = fasta_gz.to_string_lossy().into_owned();
        let mmap_dir_str = mmap_dir.to_string_lossy().into_owned();
        let stream_dir_str = stream_dir.to_string_lossy().into_owned();

        let mmap_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_gz_str.clone(),
            "--output-dir".to_string(),
            mmap_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(mmap_result.is_ok());

        let stream_result = run_with_owned_args(vec![
            "--file".to_string(),
            fasta_gz_str,
            "--mode".to_string(),
            "stream".to_string(),
            "--output-dir".to_string(),
            stream_dir_str,
            "--min-tetrads".to_string(),
            "4".to_string(),
            "--min-score".to_string(),
            "17".to_string(),
            "--overlap".to_string(),
        ]);
        assert!(stream_result.is_ok());

        for filename in [
            "chr1.g4.csv",
            "chr1.g4.overlap.csv",
            "chr1.g4.family.csv",
            "chr2.g4.csv",
            "chr2.g4.overlap.csv",
            "chr2.g4.family.csv",
        ] {
            let mmap_contents = fs::read_to_string(mmap_dir.join(filename)).unwrap();
            let stream_contents = fs::read_to_string(stream_dir.join(filename)).unwrap();
            assert_eq!(mmap_contents, stream_contents, "mismatch for {filename}");
        }

        let _ = fs::remove_file(&fasta);
        let _ = fs::remove_file(&fasta_gz);
        let _ = fs::remove_dir_all(&mmap_dir);
        let _ = fs::remove_dir_all(&stream_dir);
    }

    fn run_with_args<const N: usize>(args: [&'static str; N]) -> Result<(), String> {
        let args = args.iter().map(|arg| arg.to_string()).collect::<Vec<_>>();
        run_with_owned_args(args)
    }

    fn run_with_owned_args(args: Vec<String>) -> Result<(), String> {
        let mut argv = vec![String::from("qgrs")];
        argv.extend(args);
        let original = env::args_os().collect::<Vec<_>>();
        let _ = original;
        run_env(argv.into_iter().skip(1))
    }

    fn unique_test_path(prefix: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time before unix epoch")
            .as_nanos();
        env::temp_dir().join(format!("{prefix}_{}_{}", std::process::id(), nonce))
    }

    fn write_gzip(path: &Path, bytes: &[u8]) {
        let file = fs::File::create(path).expect("create gzip file");
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(bytes).expect("write gzip data");
        encoder.finish().expect("finish gzip");
    }
}
