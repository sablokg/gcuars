# GCUA-rs

A Rust port of the core analysis engine from **GCUA (General Codon Usage
Analysis)**

```
General Codon Usage Analysis (Rust port)

Usage: gcua <COMMAND>

Commands:
  list-codes  List all supported NCBI genetic codes
  analyze     Run a full codon-usage analysis on a FASTA file and write TSV reports
  optimize    Rewrite every sequence in a FASTA file using optimal codons
  help        Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

```

## What's included

- All 25 NCBI genetic code / translation tables (Standard, vertebrate/yeast/
  invertebrate/etc. mitochondrial codes, ciliate/euplotid/etc. nuclear
  codes), matching the Python `GENETIC_CODE_MAP`.
- A small, dependency-free FASTA parser (duplicate IDs are auto-suffixed,
  same as the original).
- Per-gene codon counts, amino acid usage, and base composition
  (GC, GC1, GC2, GC3, GC3s).
- **RSCU** (Relative Synonymous Codon Usage).
- **ENC** (Effective Number of Codons, Wright's method).
- **SCUO** (Synonymous Codon Usage Order, information-theoretic).
- **Optimal codon** identification (frequency-based, optionally against a
  reference gene subset).
- **CAI** (Codon Adaptation Index) and **Fop** (Frequency of Optimal
  Codons).
- Sequence optimization: rewrite a FASTA file so every codon is replaced
  by the optimal synonymous codon for its amino acid.
- TSV/plain-text export of all of the above.
- SVG plots via the [`kuva`](https://crates.io/crates/kuva) scientific
  plotting crate, written to `<out-dir>/plots/` by `analyze` (skip with
  `--no-plots`):
  - `gc_vs_gc3.svg` — GC vs GC3 scatter with an OLS regression line
  - `enc_wright_plot.svg` — Wright's plot (observed ENC vs GC3s, with the
    theoretical no-selection curve overlaid)
  - `cai_distribution.svg` — histogram of CAI values
  - `rscu_heatmap.svg` — genes × codons RSCU heatmap (≤60 genes), or a
    2-row mean/std-dev summary heatmap for larger datasets

## Build

**Toolchain requirement:** `kuva` needs **rustc ≥ 1.87** (it uses
`u32::is_multiple_of`, stabilized in 1.87). Check with `rustc --version`;
if you're on an older toolchain, update via `rustup update stable` (or
`rustup install stable` if you don't have rustup — see
https://rustup.rs).

```bash
cargo build --release
```

The binary is at `target/release/gcua`. This has been built and run
end-to-end in a clean environment (rustc/cargo 1.93, installed via
`apt`), including generating and visually inspecting all four plot types.

## Usage

```bash
# List supported genetic codes
gcua list-codes

# Run a full analysis
gcua analyze --fasta genes.fasta --genetic-code 11 --out-dir results

# ...optionally restrict CAI/Fop/optimal-codon calculation to a reference
# gene set (one gene name per line)
gcua analyze --fasta genes.fasta --genetic-code 11 \
    --out-dir results --reference highly_expressed.txt

# Rewrite every sequence in a FASTA file using optimal codons
gcua optimize --fasta genes.fasta --genetic-code 11 --out-dir results
```

`analyze` writes, into `--out-dir`:

- `comprehensive_metrics.tsv` — per-gene length, GC stats, ENC, CAI, Fop, SCUO
- `codon_usage.tsv` — raw codon counts, genes × 64 codons
- `rscu_values.tsv` — RSCU values, genes × 61 non-stop codons
- `optimal_codons.tsv` — one optimal codon per amino acid
- `plots/*.svg` — the four plots listed above (unless `--no-plots`)

`optimize` writes `optimized.fasta` plus the `optimal_codons.tsv` it used.

## Tested

Built and run end-to-end (`cargo build --release`, rustc/cargo 1.93)
against synthetic FASTA files:

- A 5-gene set with genetic code 11 — sane ENC (~20–22), CAI (~0.89–0.98),
  Fop (~0.88–0.95); a codon-optimized FASTA that still translates to the
  same protein sequences; all four SVG plots rendered to PNG (via
  `rsvg-convert`) and visually verified — correct regression line, correct
  Wright's-plot theoretical curve shape (peak ENC=61 at GC3s=50%), correct
  histogram, correct per-gene RSCU heatmap.
- An 80-gene random-sequence set — exercises the summary-heatmap code path
  (mean/std-dev rows instead of one row per gene); mean RSCU correctly
  centers near 1.0 as expected for unbiased random sequences.

Rendered example output from the 5-gene run is in `examples/` (both `.svg`
and rasterized `.png` for quick viewing without an SVG-capable viewer).

Gaurav Sablok \
gsablok@proton.me
