# Release Notes

## v1.8.0 - Reverse-complement strand exports

- Restored optional `--revcomp` scans for both `--base g` and `--base c` without changing the forward CSV/Parquet schema.
- Reverse hits, raw overlap hits, and family ranges are written to separate `.revcomp.*` files with coordinates mapped back to the original input and motif sequence retained in negative-strand 5′→3′ direction.
- Streaming scans spool only the current chromosome to an automatically removed temporary file, then read it backward in bounded blocks so reverse-complement output does not require chromosome-sized RAM.
- Reverse-complement validation supports ASCII IUPAC DNA symbols plus `U` and reports unsupported bytes with chromosome and 1-based input position.

## v1.7.0 - Bounded target-base run expansion

This release fixes a memory blow-up in target-base scans on pathological contigs that contain very long uninterrupted G or C runs.

### Fixed issue

- Fixed unbounded candidate growth when `--base c` scans i-motif candidates across a long C-rich contig, such as a 140.5 kb sequence with a 34,681 bp continuous C run.
- The same failure mode could affect `--base g` on G-rich inputs with very long uninterrupted G runs.
- The bug came from applying `--max-run` only to the tetrad seed size. The search still enumerated every start offset inside a long target-base run and allowed loop expansion to keep candidates whose full motif sequence contained more than `--max-run` consecutive target bases.

### Behavior change

- `--max-run` now also limits the longest continuous selected target-base run inside each partial and complete candidate motif.
- Candidates whose covered motif interval already contains more than `--max-run` consecutive target bases are pruned during expansion.
- Complete candidates are checked again before export, so reported G4/i-motif hits now match the documented target-base run limit.
- Existing results can shrink when previous outputs contained motifs with over-limit continuous G or C runs. Motifs whose G/C runs are within the configured limit are unchanged.

### Validation covered

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
git diff --check
```

## v1.6.4 - CI warning-gate compatibility

This patch release keeps the `v1.6.3` base-selectable G4/i-motif behavior and updates the codebase for the current strict Clippy warning gate used by CI.

### Fixes

- Replaced tuple comparison `sort_by` calls with `sort_by_key` in search, chunk, and stream result ordering paths.
- Derived `Default` for `QuartetBase` and marked `G` as the default variant.
- Grouped chunk window bounds into `RawSearchWindow` so the raw window scanner stays under the Clippy argument-count limit without relaxing warnings.

### Validation covered

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
git diff --check
```

## v1.6.3 - Base-selectable G4/i-motif scanning

This release replaces the old reverse-complement workflow with direct target-base scanning on the original input sequence. It is a breaking CLI and output-schema update.

### CLI changes

- Added `--base <g|c>`.
  - `--base g` scans G4 motifs by searching G tetrad runs.
  - `--base c` scans i-motif-like motifs by searching C tetrad runs directly on the original input sequence.
  - The default is `--base g`, so existing G4 scans keep the same search behavior unless they depended on renamed fields or filenames.
- Removed `--revcomp`.
  - Reverse-complement sequence generation, reverse-complement coordinate remapping, and `.revcomp.*` sidecar output are no longer part of the CLI.
  - Use `--base c` when the intended workflow is i-motif discovery.
- Replaced `--max-g-run` with `--max-run`.
  - The limit now describes the selected target base, not only G runs.
  - The old `--max-g-run` flag is rejected with migration guidance.

### Search behavior

- The existing tetrad, loop, scoring, overlap, consolidation, chunking, streaming, and circular-coordinate logic is reused for both bases.
- `--base c` does not compute a reverse complement. It scans C-rich motifs on the input strand and reports sequence text from that original input.
- Circular mode still preserves expanded coordinates, so wrap-around motifs can report `end > sequence_length`.
- mmap and stream paths use the same target-base configuration, including `--overlap` outputs.

### Output schema changes

- CSV header changed from:

```text
start,end,length,tetrads,y1,y2,y3,gscore,sequence
```

to:

```text
start,end,length,tetrads,y1,y2,y3,score,sequence
```

- Parquet field `gscore` is now `score`.
- Internal and public result field `G4.gscore` is now `G4.score`.
- The result type name remains `G4` for now.

### Output filename changes

FASTA file outputs now include the motif class in the filename:

- `--base g`: `{seqid}.g4.csv` or `{seqid}.g4.parquet`
- `--base c`: `{seqid}.i-motif.csv` or `{seqid}.i-motif.parquet`

Overlap and family sidecars follow the same motif-labeled base path:

- `{seqid}.g4.overlap.csv`
- `{seqid}.g4.family.csv`
- `{seqid}.i-motif.overlap.parquet`
- `{seqid}.i-motif.family.parquet`

Inline scans still use the explicit `--output` path provided by the caller; `--overlap` sidecars are derived from that path.

### Migration examples

G4 scan with the new run-limit name:

```bash
target/release/qgrs \
  --file genome.fa \
  --base g \
  --max-run 10 \
  --output-dir out_g4
```

i-motif-like scan on the original strand:

```bash
target/release/qgrs \
  --file genome.fa \
  --base c \
  --min-tetrads 4 \
  --min-score 17 \
  --output-dir out_i_motif \
  --overlap
```

### Validation covered

- Default scans and explicit `--base g` keep G4 search behavior while using the new `score` schema.
- `--base c` finds C tetrad motifs such as `CCCCTCCCCTCCCCTCCCC` on the original input sequence.
- Invalid `--base` values are rejected.
- `--max-g-run` is rejected with migration guidance.
- mmap and stream outputs match for C-base scans, including overlap/family outputs.
- Circular C-base scans preserve expanded coordinates.
- Verification commands used:

```bash
cargo fmt --check
cargo test
git diff --check
```
