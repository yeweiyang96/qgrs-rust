# QGRS-Rust <img alt="DNA icon" align="right" width="90" src="https://img.shields.io/badge/G4-scan-blueviolet?logo=databricks&logoColor=white" />

<div align="center">
 
![Rust 1.75+](https://img.shields.io/badge/Rust-1.75%2B-orange?logo=rust&logoColor=white)

</div>

QGRS-Rust is a ground-up Rust rewrite of the [freezer333/qgrs-cpp](https://github.com/freezer333/qgrs-cpp) project. The refactor preserves the legacy scoring math while unlocking modern performance tricks: zero-copy SIMD-friendly scanners, Rayon-powered chromosome fan-out, and dual ingestion paths (`mmap` or streaming) that keep throughput high and memory usage predictable. New functionality includes Parquet export, strict CLI validation, and a true streaming FASTA parser capable of processing multi-hundred-gigabyte `.fa` archives without ever materializing an entire chromosome in RAM. The single `qgrs` CLI now covers inline sequences, mmap batches, and the ultra-large streaming mode in one coherent interface.

## Table of contents

- [✨ Features](#-features)
- [🧰 Requirements](#-requirements)
- [⚙️ Build](#️-build)
- [🧪 Usage](#-usage)
  - [Quick recipes](#quick-recipes)
  - [CLI reference](#cli-reference)
  - [How `--max-g4-length` works](#how---max-g4-length-works)
  - [Output schema](#output-schema)
  - [Reverse-complement exports](#reverse-complement-exports---revcomp)
- [🚢 Release notes](#-release-notes)
- [✅ Testing & QA](#-testing--qa)
- [📊 Benchmarking tips](#-benchmarking-tips)

## ✨ Features

- Rust scanner mirrors the legacy scoring heuristics while benefiting from zero-copy iterators and Rayon parallelism.
- Memory-mapped (`mmap`) and streaming (`stream`) readers let you pick the best strategy per dataset.
- Optional `--base g|c` selects G4 (`g`) or i-motif (`c`) tetrad runs on the original input sequence.
- Optional `--revcomp` scans the selected motif class on the reverse-complement strand and writes separate `.revcomp.<format>` outputs.
- Optional `--circular` topology support treats each sequence/chromosome as a ring for wrap-around motif detection.
- CSV/Parquet exporters always report 1-based, inclusive coordinates for genome-browser compatibility.
- CLI validation enforces sane tetrad, loop, and window settings to avoid silent misconfiguration.
- Optional `--overlap` flag writes both the sorted raw hits and post-consolidation family ranges alongside your primary export, following `--format` (`.csv` or `.parquet`).
- FASTA outputs include the motif class in the filename: `{seqid}.g4.<format>` for `--base g` and `{seqid}.i-motif.<format>` for `--base c`.

## 🚢 Release notes

The current CLI supports direct base-selectable scans and optional reverse-complement output. See [`RELEASE.md`](RELEASE.md) for the version history.

Highlights:

- `--base g|c` selects the target tetrad base. The default remains `--base g` for G4 scans.
- `--revcomp` adds a true reverse-complement scan without changing the existing forward output or schema.
- `--max-g-run` has been replaced by `--max-run`.
- CSV and Parquet output use `score` instead of `gscore`.
- FASTA output files are motif-labeled: `{seqid}.g4.*` or `{seqid}.i-motif.*`.

## 🧰 Requirements

- Rust toolchain 1.75 or newer (install via [`rustup`](https://rustup.rs)).
- `cargo` is required for building, running, and testing.
- macOS and Linux are tested; Windows should work via WSL2.

## 🧱 Architecture

The core of QGRS-Rust lives in `src/qgrs/`, where each module maps to a distinct stage of the search pipeline:

- `data.rs`: Defines zero-copy data containers such as `ChromSequence`, `SequenceData`, and `ScanLimits`.
- `search.rs`: Implements target-base run scanning, BFS candidate expansion, scoring, and raw `G4` construction.
- `chunks.rs`: Computes windows and overlaps from `ScanLimits`, dispatches `find_raw_*`, and merges Rayon results.
- `consolidation.rs`: Deduplicates and clusters raw hits, keeping the highest `score` in each overlap family.
- `stream.rs`: Implements `StreamChromosome`/`StreamChunkScheduler` for incremental parsing of huge FASTA files.
- `loaders.rs`: Wraps mmap and regular file loaders for CLI reuse in batch mode.
- `export.rs`: Provides CSV/Parquet renderers and error types with consistent 1-based coordinate output.
- `tests/`: Centralizes unit and integration tests to ensure chunk/stream mode consistency.

`src/lib.rs` only re-exports the public API, while `src/bin/qgrs.rs` maps CLI options to the modules above to keep the entrypoint clean.

## ⚙️ Build

```bash
# clone and enter the workspace
git clone https://github.com/<your-org>/QGRS-Rust.git
cd QGRS-Rust

# dev build (fast iteration)
cargo build --bin qgrs

# optimized binary for large genomes
cargo build --release --bin qgrs

# the optimized binary lives here after a release build
target/release/qgrs --help
```

> Tip: use `cargo install --path . --bin qgrs` if you want the binary on your `~/.cargo/bin` for reuse across projects.

## 🧪 Usage

`qgrs` accepts either an inline sequence (`--sequence`) or an input file (`--file`). FASTA inputs (plain text or gzip-compressed `.gz`) are split per chromosome header, and each slice is processed independently. If you provide a file, choose either the memory-mapped (`mmap`) or buffered streaming (`stream`) pipeline with `--mode`. Pass `--base c` to scan i-motif C tetrads instead of the default G4 G tetrads. Pass `--revcomp` to retain the normal forward output and add a separate scan of the selected motif class on the reverse-complement strand. Pass `--circular` when the sequence/chromosome should be scanned as a circular molecule (wrap-around hits allowed). All examples below assume you already built the release binary (`target/release/qgrs`) or installed it as `qgrs`; use `cargo run --release --bin qgrs -- …` only when iterating locally. The banner below comes straight from `src/bin/qgrs.rs` so it always matches the binary.

```
Usage: qgrs -- [--sequence <SEQ> | --file <PATH>] [options]
Options:
   --sequence <SEQ>       Inline DNA/RNA sequence to scan
   --file <PATH>          Read sequences from FASTA/FASTA.gz (chromosomes split independently)
   --min-tetrads <N>      Minimum tetrads to seed (default 2)
   --min-score <S>        Minimum score (default 17)
   --base <g|c>           Tetrad base to scan: g for G4, c for i-motif (default g)
   --max-run <N>          Maximum allowed target-base run length (default 10)
   --max-g4-length <N>    Maximum allowed G4 length in bp (default 45)
   --format <csv|parquet> Output format (default csv)
   --output <PATH>        Destination file when using --sequence (required for parquet)
   --output-dir <DIR>     Directory for per-chromosome exports when using --file
   --mode <mmap|stream>   Input mode when using --file (default mmap)
   --overlap              Also emit raw hits and family ranges beside each primary output
   --revcomp              Also scan reverse-complement and emit .revcomp.<format>
   --circular             Treat each sequence/chromosome as circular
   --help                 Show this message
```

### Quick recipes

```bash
# 1. Quick sanity check against a short inline sequence
target/release/qgrs \
   --sequence GGGGAGGGGAGGGGAGGGG \
   --min-tetrads 4 \
   --min-score 17 \
   --format csv

# 2. Process a FASTA file using mmap and emit CSV files per chromosome
target/release/qgrs \
   --file data/genome.fa \
   --mode mmap \
   --min-tetrads 4 \
   --min-score 20 \
   --output-dir ./qgrs_csv

# 3. Stream extremely large FASTA (plain or .gz) with Parquet output
target/release/qgrs \
   --file hg38.fa.gz \
   --mode stream \
   --min-tetrads 3 \
   --max-run 12 \
   --max-g4-length 60 \
   --format parquet \
   --output-dir ./qgrs_parquet

# 4. Clamp loop and run lengths for custom heuristics while writing to stdout
target/release/qgrs \
   --sequence GGGGTTTTGGGGTTTTGGGGTTTTGGGG \
   --max-run 8 \
   --max-g4-length 45 \
   --output -

# 5. Keep raw hits and family ranges for post-processing
target/release/qgrs \
   --file data/genome.fa \
   --mode stream \
   --output-dir ./qgrs_debug \
   --overlap

# 6. Enable circular topology (wrap-around hits)
target/release/qgrs \
   --sequence GAGGGGAGGGGAGGGGGGG \
   --min-tetrads 4 \
   --circular

# 7. Scan i-motif C tetrads on the original input sequence
target/release/qgrs \
   --file data/genome.fa \
   --mode stream \
   --output-dir ./qgrs_i_motif \
   --format parquet \
   --base c \
   --overlap

# 8. Emit forward and reverse-strand G4 results
target/release/qgrs \
   --file data/genome.fa \
   --mode stream \
   --output-dir ./qgrs_both_strands \
   --format parquet \
   --revcomp \
   --overlap
```

### CLI reference

| Flag                      | Description                                                                                | Default                  |
| ------------------------- | ------------------------------------------------------------------------------------------ | ------------------------ |
| `--sequence <SEQ>`        | Inline DNA sequence to scan (mutually exclusive with `--file`).                            | _none_                   |
| `--file <PATH>`           | FASTA input path (plain text or gzip-compressed `.gz`) containing one or more sequences.   | _none_                   |
| `--mode <mmap\|stream>`   | File ingestion strategy; `mmap` favors fast disks, `stream` lowers RAM.                    | `mmap`                   |
| `--min-tetrads <INT>`     | Minimum number of stacked tetrads required for a hit.                                      | `2`                      |
| `--min-score <INT>`       | Minimum score threshold.                                                                   | `17`                     |
| `--base <g\|c>`           | Tetrad base to scan: `g` for G4 or `c` for i-motif.                                        | `g`                      |
| `--max-run <INT>`         | Upper bound for contiguous target-base run length (must be ≥ `min-tetrads`).               | `10`                     |
| `--max-g4-length <INT>`   | Upper bound for the full quadruplex length (must be ≥ `4 * min_tetrads`).                  | `45`                     |
| `--format <csv\|parquet>` | Output encoding. CSV defaults to stdout for inline sequences; Parquet requires a file/dir. | `csv`                    |
| `--output <FILE\|- >`     | Single output file (or `-` for stdout) when scanning inline sequences.                     | stdout for CSV           |
| `--output-dir <DIR>`      | Directory for per-chromosome files when reading FASTA/plain inputs. File names are `{seqid}.g4.<format>` or `{seqid}.i-motif.<format>`. | _required with `--file`_ |
| `--overlap`               | Emit `{seqid}.{motif}.overlap.<format>` (raw hits) and `{seqid}.{motif}.family.<format>` (family ranges) per FASTA output file. | off                      |
| `--revcomp`               | Also scan the reverse-complement strand and emit `{seqid}.{motif}.revcomp.<format>`. | off                      |
| `--circular`              | Treat each sequence/chromosome as circular; wrap-around hits keep expanded coordinates in output, so `end` may exceed chromosome length `N`. | off                      |

The CLI aborts with a descriptive error if incompatible parameters are provided (e.g., `--mode stream` without `--file`, `--base a`, or `--max-run < min-tetrads`). When scanning files you must pass `--output-dir`; when `--overlap` is enabled for inline scans, `--output` is required so sidecar files can be named deterministically. Inline `--revcomp` also requires a real output path and cannot write its two result streams to stdout.

### How `--max-g4-length` works

`--max-g4-length` affects more than final hit filtering. It participates in candidate seeding, loop expansion, viability checks, score calculation, chunk overlap, and circular wrap-around buffering.

- Candidate seeding limits tetrads to `min(max_run, floor(max_g4_length / 4))`, so smaller values can eliminate high-tetrad candidates before BFS expansion starts.
- Each candidate does not use `max_g4_length` directly. Instead, it uses `min(legacy_cap, max_g4_length)`, where `legacy_cap = 30` for `tetrads < 3` and `legacy_cap = 45` for `tetrads >= 3`.
- As a result, increasing `--max-g4-length` above `30` does not further relax 2-tetrad scoring/length checks, and increasing it above `45` does not further relax 3+-tetrad scoring/length checks.
- Decreasing `--max-g4-length` below those legacy caps reduces the allowed total motif length, narrows the loop search space, lowers the score ceiling, and can remove candidates entirely.

With the default `--max-g4-length 45`, the effective per-candidate limits are:

| tetrads | effective candidate `max_length` | practical effect |
| --- | --- | --- |
| `< 3` | `30` | 2-tetrad motifs are still scored and filtered against a 30 bp legacy cap |
| `>= 3` | `45` | 3+-tetrad motifs use the 45 bp legacy cap |

If you lower the flag, replace those values with `min(legacy_cap, --max-g4-length)`. For example, with `--max-g4-length 32`, 2-tetrad candidates still use `30`, while 3+-tetrad candidates drop from `45` to `32`.

### Output schema

Both exporters emit the same fields (see `render_csv` and `write_parquet_from_results` in `src/qgrs/mod.rs`):

| Column           | Meaning                                                                                 |
| ---------------- | --------------------------------------------------------------------------------------- |
| `start`          | 1-based inclusive start coordinate of the hit within the processed sequence/chromosome. |
| `end`            | 1-based inclusive end coordinate. In circular mode, wrap-around hits keep expanded coordinates, so `end` may be larger than the chromosome length `N`. |
| `length`         | Total number of bases spanned by the quadruplex.                                        |
| `tetrads`        | Count of stacked tetrads contributing to the hit.                                       |
| `y1`, `y2`, `y3` | Loop lengths between successive target-base runs (0 means no spacer).                  |
| `score`          | Score used for filtering and ranking candidates.                                       |
| `sequence`       | Exact motif sequence extracted from the input.                                         |

CSV output always includes the header `start,end,length,tetrads,y1,y2,y3,score,sequence`. When scanning FASTA inputs, each chromosome is written to its own motif-labeled file such as `chr1.g4.csv` or `chr1.i-motif.csv` (so the filename, not a column, captures the chromosome name and motif class). Parquet exports contain the same columns using Arrow types (`UInt64` for coordinates/lengths, `Int32` for loop lengths and score, and UTF-8 for sequences). In circular mode, CLI exports keep the same expanded-coordinate representation used internally, so wrap-around motifs can appear with `end > N`.

### Reverse-complement exports (`--revcomp`)

`--revcomp` leaves every forward output unchanged and adds `{base}.revcomp.<format>`. The reverse output uses the same schema, thresholds, selected `--base`, topology, and consolidation logic as the forward scan.

- `start` and `end` are mapped back to 1-based coordinates on the original input. Linear hits use `start = N - end_rc + 1` and `end = N - start_rc + 1`.
- Circular hits keep expanded coordinates after mapping, so a cross-origin reverse hit can also report `end > N`.
- `sequence` remains in reverse-complement scan direction (negative-strand 5′→3′).
- With `--overlap`, the reverse scan also writes `{base}.revcomp.overlap.<format>` and `{base}.revcomp.family.<format>`.
- `--base g --revcomp` reports negative-strand G4 candidates; `--base c --revcomp` reports negative-strand i-motif candidates.

In stream mode, each normalized chromosome is spooled to an automatically deleted temporary file and read backward in bounded blocks for the reverse scan. RAM remains bounded, while temporary disk usage peaks at approximately the size of the largest chromosome. Set the standard `TMPDIR` environment variable when those files must reside on a specific volume.

### Overlap exports (`--overlap`)

Pass `--overlap` to retain additional debugging artifacts for every output file:

- **Raw hits**: `{seqid}.{motif}.overlap.<format>` mirrors the primary result schema but contains the full pre-consolidation hit list. This lets you diff against other implementations or inspect families before winners are picked.
- **Family ranges**: `{seqid}.{motif}.family.<format>` lists `family_index,start,end` for each consolidated family, using the same 1-based inclusive coordinates. The index column reflects the order in which families were discovered.

For inline scans you must also supply `--output`, because the overlap files reuse that explicit base path. When scanning FASTA files, each chromosome inherits the motif-labeled filename that would have been written normally (for example, `chr2.i-motif.parquet` also writes `chr2.i-motif.overlap.parquet` and `chr2.i-motif.family.parquet`). In streaming mode the extra files are flushed as soon as each chromosome finishes, so the memory footprint stays bounded even for gigantic inputs.

## Testing & QA

```bash
# unit tests
cargo test

# lint + formatting (optional but recommended before sending patches)
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
```

## Benchmarking tips

```bash
time target/release/qgrs --file aaa.fa --mode mmap   --max-g4-length 32 --max-run 8  --output-dir out-mmap
time target/release/qgrs --file aaa.fa --mode stream --max-g4-length 45 --max-run 10 --output-dir out-stream
```

Track `real` time, CPU%, and RSS with your preferred profiler to decide whether `mmap` or `stream` is better for your environment. Always benchmark with `--release` builds to enable full optimizations.

### compare_modes consistency tester

`target/release/compare_modes` (defined in `src/bin/compare_modes.rs`) benchmarks and cross-checks the two ingestion pipelines against the same FASTA input. It scans every chromosome once with the mmap batch loader and once with the streaming reader, reports per-mode timings and hit counts, then diff-checks every field (`start`, `end`, `length`, loops, tetrads, score, sequence) to ensure both paths stay bit-for-bit aligned. The process exits with code `0` on success and `1` with detailed mismatch logs when discrepancies are detected.

**Before you run it**

- Build the binary at least once so that `target/release/compare_modes` exists (typically from a prior release build).
- Provide a readable FASTA file containing one or more chromosomes; the tool operates read-only and prints summaries to stdout.

**CLI arguments**

| Positional argument | Description                         | Default |
| ------------------- | ----------------------------------- | ------- |
| `<fasta>`           | Path to the FASTA file to compare.  | _none_  |
| `[min_tetrads]`     | Minimum tetrads threshold per scan. | `2`     |
| `[min_score]`       | Minimum score threshold per scan.   | `17`    |

**Example commands** (using the compiled release binary; adjust the path if you install it elsewhere):

```bash
# Compare dme.fa with default thresholds
target/release/compare_modes dme.fa

# Tighten heuristics to 3 tetrads / score 20
target/release/compare_modes dme.fa 3 20

# Point at a chromosome subset file
target/release/compare_modes output/chromosome-2L.fa 4 30
```

During a run you will see individual sections for the mmap phase, the stream phase, a speed comparison, and the final consistency verdict. An error summary (up to 10 detailed mismatches) is printed before the program returns a non-zero exit status, which makes the tool suitable for automated regression checks.
