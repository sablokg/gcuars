mod data;
mod fasta;
mod genetic_code;
mod plots;
use anyhow::Result;
use clap::{Parser, Subcommand};
use data::GcuaData;
use genetic_code::{GENETIC_CODES, GENETIC_CODE_NAMES};
use std::fs;
use std::path::PathBuf;

/*
Gaurav Sablok
gsablok@proton.me
 */

const VERSION: &str = "0.1.0";
const BANNER: &str = r#"
*******************************************************
 GCUA-rs: General Codon Usage Analysis (Rust port)
 Gaurav Sablok gsablok@proton.me
 Original: McInerney JO. Bioinformatics. 1998;14(4):372-3.
*******************************************************
"#;

#[derive(Parser)]
#[command(name = "gcua", version = VERSION, about = "General Codon Usage Analysis (Rust port)")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// List all supported NCBI genetic codes
    ListCodes,

    /// Run a full codon-usage analysis on a FASTA file and write TSV reports
    Analyze {
        /// Path to a FASTA file of protein-coding DNA sequences
        #[arg(short, long)]
        fasta: PathBuf,
        /// NCBI genetic code / translation table id (see `list-codes`)
        #[arg(short, long, default_value_t = 1)]
        genetic_code: u8,
        /// Output directory for TSV reports
        #[arg(short, long, default_value = "gcua_outputs")]
        out_dir: PathBuf,
        /// Optional file with one reference gene name per line (used for
        /// CAI / Fop / optimal-codon calculation). Defaults to all genes.
        #[arg(short, long)]
        reference: Option<PathBuf>,
        /// Skip SVG plot generation (plots are written by default)
        #[arg(long, default_value_t = false)]
        no_plots: bool,
    },

    /// Rewrite every sequence in a FASTA file using optimal codons
    Optimize {
        #[arg(short, long)]
        fasta: PathBuf,
        #[arg(short, long, default_value_t = 1)]
        genetic_code: u8,
        #[arg(short, long, default_value = "gcua_outputs")]
        out_dir: PathBuf,
        #[arg(short, long)]
        reference: Option<PathBuf>,
    },
}

fn load_reference(path: &Option<PathBuf>) -> Result<Option<Vec<String>>> {
    match path {
        None => Ok(None),
        Some(p) => {
            let content = fs::read_to_string(p)?;
            Ok(Some(
                content
                    .lines()
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect(),
            ))
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::ListCodes => {
            println!("{BANNER}");
            let mut ids: Vec<&u8> = GENETIC_CODES.keys().collect();
            ids.sort();
            for id in ids {
                println!(
                    "  [{:>2}] {}",
                    id,
                    GENETIC_CODE_NAMES.get(id).copied().unwrap_or("Unknown")
                );
            }
        }

        Commands::Analyze {
            fasta,
            genetic_code,
            out_dir,
            reference,
            no_plots,
        } => {
            println!("{BANNER}");
            fs::create_dir_all(&out_dir)?;

            let mut data = GcuaData::new(genetic_code)?;
            println!(
                "Genetic code: [{}] {}",
                data.genetic_code,
                data.genetic_code_name()
            );
            println!("Loading {}...", fasta.display());
            data.load_fasta(fasta.to_str().unwrap())?;
            println!("Loaded {} sequences.", data.gene_names.len());

            let reference_genes = load_reference(&reference)?;
            let reference_slice = reference_genes.as_deref();
            if let Some(r) = reference_slice {
                println!(
                    "Using {} reference genes for CAI/Fop/optimal codons.",
                    r.len()
                );
            }

            println!("Calculating RSCU...");
            data.calc_rscu();
            println!("Calculating amino acid usage...");
            data.calc_aa_usage();
            println!("Calculating ENC...");
            data.calc_enc();
            println!("Calculating SCUO...");
            data.calc_scuo();
            println!("Identifying optimal codons...");
            data.calc_optimal_codons(reference_slice);
            println!("Calculating CAI...");
            data.calc_cai(reference_slice);
            println!("Calculating Fop...");
            data.calc_fop(reference_slice);

            data.export_comprehensive(out_dir.join("comprehensive_metrics.tsv").to_str().unwrap())?;
            data.export_rscu(out_dir.join("rscu_values.tsv").to_str().unwrap())?;
            data.export_codon_usage(out_dir.join("codon_usage.tsv").to_str().unwrap())?;
            data.export_optimal_codons(out_dir.join("optimal_codons.tsv").to_str().unwrap())?;

            if !no_plots {
                println!("Generating plots...");
                let plots_dir = out_dir.join("plots");
                let written = plots::write_all_plots(&data, &plots_dir)?;
                for path in &written {
                    println!("  wrote {}", path.display());
                }
            }

            // Quick summary to stdout
            if let Some(enc) = &data.enc {
                let vals: Vec<f64> = enc.values().copied().collect();
                print_summary("ENC", &vals);
            }
            if let Some(cai) = &data.cai {
                let vals: Vec<f64> = cai.values().copied().collect();
                print_summary("CAI", &vals);
            }
            if let Some(fop) = &data.fop {
                let vals: Vec<f64> = fop.values().copied().collect();
                print_summary("Fop", &vals);
            }

            println!("\nDone. Results written to {}", out_dir.display());
        }

        Commands::Optimize {
            fasta,
            genetic_code,
            out_dir,
            reference,
        } => {
            println!("{BANNER}");
            fs::create_dir_all(&out_dir)?;

            let mut data = GcuaData::new(genetic_code)?;
            data.load_fasta(fasta.to_str().unwrap())?;
            println!("Loaded {} sequences.", data.gene_names.len());

            let reference_genes = load_reference(&reference)?;
            data.calc_optimal_codons(reference_genes.as_deref());

            let out_path = out_dir.join("optimized.fasta");
            let mut out = String::new();
            let mut n_optimized = 0usize;
            for gene in &data.gene_names {
                if let Some(seq) = data.optimize_gene(gene) {
                    out.push_str(&format!(">{gene}_optimized\n{seq}\n"));
                    n_optimized += 1;
                }
            }
            fs::write(&out_path, out)?;
            data.export_optimal_codons(out_dir.join("optimal_codons.tsv").to_str().unwrap())?;

            println!(
                "Optimized {n_optimized} sequences. Written to {}",
                out_path.display()
            );
        }
    }

    Ok(())
}

fn print_summary(name: &str, values: &[f64]) {
    if values.is_empty() {
        return;
    }
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("  {name}: mean={mean:.4}  min={min:.4}  max={max:.4}");
}
