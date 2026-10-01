//! Core codon-usage analysis engine. This is the Rust equivalent of the
//! `GCUAData` class in the Python original: it loads sequences, counts
//! codons, and derives all of the standard codon-usage-bias metrics.

/*
Gaurav Sablok
gsablok@proton.me
 */

use crate::fasta::read_fasta;
use crate::genetic_code::{aa_to_codons, dna_codon, CODONS, GENETIC_CODES, GENETIC_CODE_NAMES};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::io::Write;

#[derive(Default, Clone, Debug)]
pub struct BaseComposition {
    pub length: usize,
    pub aa_count: usize,
    pub a: u32,
    pub t: u32,
    pub g: u32,
    pub c: u32,
    pub gc: f64,
    pub gc1: f64,
    pub gc2: f64,
    pub gc3: f64,
    pub gc3s: f64,
}

/// Dinucleotide families whose third position is 4-fold degenerate
/// (used to compute GC3s, mirroring the Python `fourfold_codons` set).
const FOURFOLD_FAMILIES: [&str; 7] = ["CU", "CC", "CG", "AC", "GU", "GC", "GG"];

pub struct GcuaData {
    pub genetic_code: u8,
    pub code: HashMap<&'static str, &'static str>,
    pub aa_map: HashMap<String, Vec<String>>, // amino acid -> synonymous codons
    pub gene_names: Vec<String>,
    pub sequences: HashMap<String, String>, // raw (uppercase) DNA, kept for optimization
    pub codon_counts: HashMap<String, HashMap<String, u32>>,
    pub base_composition: HashMap<String, BaseComposition>,
    pub rscu: HashMap<String, HashMap<String, f64>>,
    pub aa_usage: HashMap<String, HashMap<String, u32>>,
    pub enc: Option<HashMap<String, f64>>,
    pub cai: Option<HashMap<String, f64>>,
    pub fop: Option<HashMap<String, f64>>,
    pub scuo: Option<HashMap<String, f64>>,
    pub optimal_codons: Option<HashMap<String, String>>,
}

impl GcuaData {
    pub fn new(genetic_code: u8) -> Result<Self> {
        let code = GENETIC_CODES
            .get(&genetic_code)
            .cloned()
            .with_context(|| format!("Unknown genetic code id: {genetic_code}"))?;
        let aa_map = aa_to_codons(&code);
        Ok(Self {
            genetic_code,
            code,
            aa_map,
            gene_names: Vec::new(),
            sequences: HashMap::new(),
            codon_counts: HashMap::new(),
            base_composition: HashMap::new(),
            rscu: HashMap::new(),
            aa_usage: HashMap::new(),
            enc: None,
            cai: None,
            fop: None,
            scuo: None,
            optimal_codons: None,
        })
    }

    pub fn genetic_code_name(&self) -> &'static str {
        GENETIC_CODE_NAMES
            .get(&self.genetic_code)
            .copied()
            .unwrap_or("Unknown")
    }

    // -------------------------------------------------------------------
    // Loading & basic per-sequence stats
    // -------------------------------------------------------------------

    pub fn load_fasta(&mut self, path: &str) -> Result<()> {
        let records = read_fasta(path)?;
        if records.is_empty() {
            anyhow::bail!("No sequences found in {path}");
        }
        for rec in records {
            self.gene_names.push(rec.id.clone());
            let upper = rec.sequence.to_uppercase();
            self.process_sequence(&rec.id, &upper);
            self.sequences.insert(rec.id, upper);
        }
        Ok(())
    }

    fn process_sequence(&mut self, gene_id: &str, seq: &str) {
        let dna = seq.as_bytes();
        let usable_len = dna.len() - (dna.len() % 3);

        let mut counts: HashMap<String, u32> = HashMap::new();
        let mut base_counts = [0u32; 4]; // A T G C
        let mut pos1 = [0u32; 4];
        let mut pos2 = [0u32; 4];
        let mut pos3 = [0u32; 4];
        let mut pos3_4fold = [0u32; 4];
        let mut aa_count = 0usize;

        let idx_of = |b: u8| -> Option<usize> {
            match b {
                b'A' => Some(0),
                b'T' => Some(1),
                b'G' => Some(2),
                b'C' => Some(3),
                _ => None,
            }
        };

        let mut i = 0;
        while i + 3 <= usable_len {
            let codon_dna = &dna[i..i + 3];
            i += 3;
            if codon_dna.contains(&b'N') {
                continue;
            }
            for (pos, &b) in codon_dna.iter().enumerate() {
                if let Some(bi) = idx_of(b) {
                    base_counts[bi] += 1;
                    match pos {
                        0 => pos1[bi] += 1,
                        1 => pos2[bi] += 1,
                        2 => {
                            pos3[bi] += 1;
                            let family = [codon_dna[0] as char, codon_dna[1] as char]
                                .iter()
                                .collect::<String>();
                            if FOURFOLD_FAMILIES.contains(&family.as_str()) {
                                pos3_4fold[bi] += 1;
                            }
                        }
                        _ => {}
                    }
                }
            }
            let codon_str: String = codon_dna.iter().map(|&b| b as char).collect();
            let codon_rna = codon_str.replace('T', "U");
            if self.code.contains_key(codon_rna.as_str()) {
                *counts.entry(codon_rna).or_insert(0) += 1;
                aa_count += 1;
            }
        }

        let total = base_counts.iter().sum::<u32>() as f64;
        let gc = pct(base_counts[2] + base_counts[3], total as u32);
        let gc1 = pct(pos1[2] + pos1[3], pos1.iter().sum());
        let gc2 = pct(pos2[2] + pos2[3], pos2.iter().sum());
        let gc3 = pct(pos3[2] + pos3[3], pos3.iter().sum());
        let gc3s = pct(pos3_4fold[2] + pos3_4fold[3], pos3_4fold.iter().sum());

        let bc = BaseComposition {
            length: usable_len,
            aa_count,
            a: base_counts[0],
            t: base_counts[1],
            g: base_counts[2],
            c: base_counts[3],
            gc,
            gc1,
            gc2,
            gc3,
            gc3s,
        };

        self.codon_counts.insert(gene_id.to_string(), counts);
        self.base_composition.insert(gene_id.to_string(), bc);
    }

    fn gene_codon_count(&self, gene: &str, codon: &str) -> u32 {
        self.codon_counts
            .get(gene)
            .and_then(|m| m.get(codon))
            .copied()
            .unwrap_or(0)
    }

    // -------------------------------------------------------------------
    // Metrics
    // -------------------------------------------------------------------

    /// Relative Synonymous Codon Usage.
    pub fn calc_rscu(&mut self) {
        let genes = self.gene_names.clone();
        for gene in genes {
            let mut gene_rscu: HashMap<String, f64> = HashMap::new();
            for (aa, codons) in self.aa_map.iter() {
                if aa == "STOP" {
                    continue;
                }
                let n = codons.len() as f64;
                let total: u32 = codons.iter().map(|c| self.gene_codon_count(&gene, c)).sum();
                if total == 0 {
                    continue;
                }
                for c in codons {
                    let count = self.gene_codon_count(&gene, c);
                    gene_rscu.insert(c.clone(), count as f64 * n / total as f64);
                }
            }
            self.rscu.insert(gene, gene_rscu);
        }
    }

    /// Amino acid usage (sum of synonymous codon counts per amino acid).
    pub fn calc_aa_usage(&mut self) {
        let genes = self.gene_names.clone();
        for gene in genes {
            let mut usage: HashMap<String, u32> = HashMap::new();
            for (aa, codons) in self.aa_map.iter() {
                let total: u32 = codons.iter().map(|c| self.gene_codon_count(&gene, c)).sum();
                usage.insert(aa.clone(), total);
            }
            self.aa_usage.insert(gene, usage);
        }
    }

    /// Effective Number of Codons (Wright 1990).
    pub fn calc_enc(&mut self) {
        let mut result = HashMap::new();
        for gene in &self.gene_names {
            let mut by_degeneracy: HashMap<usize, Vec<f64>> = HashMap::new();
            for (aa, codons) in self.aa_map.iter() {
                if aa == "STOP" {
                    continue;
                }
                let degeneracy = codons.len();
                if degeneracy <= 1 {
                    continue;
                }
                let total: u32 = codons.iter().map(|c| self.gene_codon_count(gene, c)).sum();
                if total == 0 {
                    continue;
                }
                let homozygosity: f64 = codons
                    .iter()
                    .map(|c| {
                        let f = self.gene_codon_count(gene, c) as f64 / total as f64;
                        f * f
                    })
                    .sum();
                by_degeneracy
                    .entry(degeneracy)
                    .or_default()
                    .push(homozygosity);
            }
            let f_for = |d: usize| -> f64 {
                by_degeneracy
                    .get(&d)
                    .map(|v| v.iter().sum::<f64>() / v.len() as f64)
                    .unwrap_or(0.0)
            };
            let (f2, f3, f4, f6) = (f_for(2), f_for(3), f_for(4), f_for(6));
            let mut enc = 2.0;
            if f2 > 0.0 {
                enc += 9.0 / f2;
            }
            if f3 > 0.0 {
                enc += 1.0 / f3;
            }
            if f4 > 0.0 {
                enc += 5.0 / f4;
            }
            if f6 > 0.0 {
                enc += 3.0 / f6;
            }
            result.insert(gene.clone(), enc.min(61.0));
        }
        self.enc = Some(result);
    }

    /// Synonymous Codon Usage Order (information-theoretic bias measure).
    pub fn calc_scuo(&mut self) {
        let mut result = HashMap::new();
        for gene in &self.gene_names {
            let mut aa_scuo: HashMap<String, f64> = HashMap::new();
            let mut aa_weight: HashMap<String, u32> = HashMap::new();
            for (aa, codons) in self.aa_map.iter() {
                if aa == "STOP" || codons.len() <= 1 {
                    continue;
                }
                let total: u32 = codons.iter().map(|c| self.gene_codon_count(gene, c)).sum();
                if total == 0 {
                    continue;
                }
                let entropy: f64 = codons
                    .iter()
                    .map(|c| {
                        let f = self.gene_codon_count(gene, c) as f64 / total as f64;
                        if f > 0.0 {
                            -f * f.ln()
                        } else {
                            0.0
                        }
                    })
                    .sum();
                let max_entropy = (codons.len() as f64).ln();
                let scuo = if max_entropy > 0.0 {
                    1.0 - entropy / max_entropy
                } else {
                    0.0
                };
                aa_scuo.insert(aa.clone(), scuo);
                aa_weight.insert(aa.clone(), total);
            }
            let total_weight: u32 = aa_weight.values().sum();
            let scuo_val = if total_weight > 0 {
                aa_scuo
                    .iter()
                    .map(|(aa, v)| v * (*aa_weight.get(aa).unwrap() as f64))
                    .sum::<f64>()
                    / total_weight as f64
            } else {
                0.0
            };
            result.insert(gene.clone(), scuo_val);
        }
        self.scuo = Some(result);
    }

    /// Identify the most frequently used codon for each amino acid,
    /// optionally restricted to a reference gene set.
    pub fn calc_optimal_codons(&mut self, reference: Option<&[String]>) {
        let genes: Vec<String> = reference
            .map(|r| r.to_vec())
            .unwrap_or_else(|| self.gene_names.clone());

        let mut ref_counts: HashMap<&str, u32> = HashMap::new();
        for codon in CODONS.iter() {
            let total: u32 = genes.iter().map(|g| self.gene_codon_count(g, codon)).sum();
            ref_counts.insert(codon, total);
        }

        let mut optimal = HashMap::new();
        for (aa, codons) in self.aa_map.iter() {
            if aa == "STOP" {
                continue;
            }
            if let Some(best) = codons
                .iter()
                .max_by_key(|c| ref_counts.get(c.as_str()).copied().unwrap_or(0))
            {
                optimal.insert(aa.clone(), best.clone());
            }
        }
        self.optimal_codons = Some(optimal);
    }

    /// Codon Adaptation Index (Sharp & Li 1987).
    pub fn calc_cai(&mut self, reference: Option<&[String]>) {
        let genes: Vec<String> = reference
            .map(|r| r.to_vec())
            .unwrap_or_else(|| self.gene_names.clone());

        let mut ref_counts: HashMap<&str, u32> = HashMap::new();
        for codon in CODONS.iter() {
            let total: u32 = genes.iter().map(|g| self.gene_codon_count(g, codon)).sum();
            ref_counts.insert(codon, total);
        }

        let mut weights: HashMap<&str, f64> = HashMap::new();
        for (aa, codons) in self.aa_map.iter() {
            if aa == "STOP" {
                continue;
            }
            let max_count = codons
                .iter()
                .map(|c| ref_counts.get(c.as_str()).copied().unwrap_or(0))
                .max()
                .unwrap_or(0);
            for c in codons {
                let count = ref_counts.get(c.as_str()).copied().unwrap_or(0);
                let w = if count > 0 && max_count > 0 {
                    count as f64 / max_count as f64
                } else {
                    0.001
                };
                weights.insert(CODONS.iter().find(|x| **x == c.as_str()).unwrap(), w);
            }
        }

        let mut result = HashMap::new();
        for gene in &self.gene_names {
            let mut weighted_sum = 0.0;
            let mut codon_count = 0u32;
            for codon in CODONS.iter() {
                let aa = self.code.get(codon).copied().unwrap_or("");
                if aa == "STOP" {
                    continue;
                }
                let count = self.gene_codon_count(gene, codon);
                if count > 0 {
                    if let Some(&w) = weights.get(codon) {
                        if w > 0.0 {
                            weighted_sum += count as f64 * w.ln();
                            codon_count += count;
                        }
                    }
                }
            }
            let cai_val = if codon_count > 0 {
                (weighted_sum / codon_count as f64).exp()
            } else {
                0.0
            };
            result.insert(gene.clone(), cai_val);
        }
        self.cai = Some(result);
    }

    /// Frequency of Optimal Codons.
    pub fn calc_fop(&mut self, reference: Option<&[String]>) {
        if self.optimal_codons.is_none() {
            self.calc_optimal_codons(reference);
        }
        let optimal = self.optimal_codons.clone().unwrap();

        let mut result = HashMap::new();
        for gene in &self.gene_names {
            let mut opt_count = 0u32;
            let mut syn_count = 0u32;
            for (aa, codons) in self.aa_map.iter() {
                if aa == "STOP" || codons.len() <= 1 {
                    continue;
                }
                for c in codons {
                    let count = self.gene_codon_count(gene, c);
                    syn_count += count;
                    if optimal.get(aa) == Some(c) {
                        opt_count += count;
                    }
                }
            }
            let fop_val = if syn_count > 0 {
                opt_count as f64 / syn_count as f64
            } else {
                0.0
            };
            result.insert(gene.clone(), fop_val);
        }
        self.fop = Some(result);
    }

    // -------------------------------------------------------------------
    // Sequence optimization
    // -------------------------------------------------------------------

    /// Rewrite one gene's sequence, replacing each codon with the optimal
    /// synonymous codon for its amino acid (stop codons are left as-is).
    pub fn optimize_gene(&self, gene: &str) -> Option<String> {
        let optimal = self.optimal_codons.as_ref()?;
        let seq = self.sequences.get(gene)?;
        let bytes = seq.as_bytes();
        let mut out = String::with_capacity(bytes.len());
        let mut i = 0;
        while i + 3 <= bytes.len() {
            let codon_dna: String = bytes[i..i + 3].iter().map(|&b| b as char).collect();
            if codon_dna.contains('N') {
                out.push_str(&codon_dna);
                i += 3;
                continue;
            }
            let codon_rna = codon_dna.replace('T', "U");
            match self.code.get(codon_rna.as_str()).copied() {
                Some("STOP") => out.push_str(&codon_dna),
                Some(aa) => {
                    if let Some(opt_rna) = optimal.get(aa) {
                        out.push_str(&dna_codon(opt_rna));
                    } else {
                        out.push_str(&codon_dna);
                    }
                }
                None => out.push_str(&codon_dna),
            }
            i += 3;
        }
        if i < bytes.len() {
            out.push_str(&String::from_utf8_lossy(&bytes[i..]));
        }
        Some(out)
    }

    // -------------------------------------------------------------------
    // Export
    // -------------------------------------------------------------------

    pub fn export_comprehensive(&self, path: &str) -> Result<()> {
        let mut f = fs::File::create(path)?;
        writeln!(f, "# GCUA-rs comprehensive metrics")?;
        writeln!(
            f,
            "# Genetic code: [{}] {}",
            self.genetic_code,
            self.genetic_code_name()
        )?;
        writeln!(f, "# Total genes: {}", self.gene_names.len())?;
        write!(f, "Gene\tLength\tGC\tGC1\tGC2\tGC3\tGC3s")?;
        if self.enc.is_some() {
            write!(f, "\tENC")?;
        }
        if self.cai.is_some() {
            write!(f, "\tCAI")?;
        }
        if self.fop.is_some() {
            write!(f, "\tFop")?;
        }
        if self.scuo.is_some() {
            write!(f, "\tSCUO")?;
        }
        writeln!(f)?;

        for gene in &self.gene_names {
            let bc = &self.base_composition[gene];
            write!(
                f,
                "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                gene, bc.length, bc.gc, bc.gc1, bc.gc2, bc.gc3, bc.gc3s
            )?;
            if let Some(enc) = &self.enc {
                write!(f, "\t{:.3}", enc[gene])?;
            }
            if let Some(cai) = &self.cai {
                write!(f, "\t{:.4}", cai[gene])?;
            }
            if let Some(fop) = &self.fop {
                write!(f, "\t{:.4}", fop[gene])?;
            }
            if let Some(scuo) = &self.scuo {
                write!(f, "\t{:.4}", scuo[gene])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }

    pub fn export_codon_usage(&self, path: &str) -> Result<()> {
        let mut f = fs::File::create(path)?;
        write!(f, "Gene")?;
        for codon in CODONS.iter() {
            write!(f, "\t{codon}")?;
        }
        writeln!(f)?;
        for gene in &self.gene_names {
            write!(f, "{gene}")?;
            for codon in CODONS.iter() {
                write!(f, "\t{}", self.gene_codon_count(gene, codon))?;
            }
            writeln!(f)?;
        }
        Ok(())
    }

    pub fn export_rscu(&self, path: &str) -> Result<()> {
        let mut f = fs::File::create(path)?;
        let non_stop: Vec<&str> = CODONS
            .iter()
            .filter(|c| self.code.get(*c).copied() != Some("STOP"))
            .copied()
            .collect();
        write!(f, "Gene")?;
        for codon in &non_stop {
            write!(f, "\t{codon}")?;
        }
        writeln!(f)?;
        for gene in &self.gene_names {
            write!(f, "{gene}")?;
            let gene_rscu = self.rscu.get(gene);
            for codon in &non_stop {
                let val = gene_rscu
                    .and_then(|m| m.get(*codon))
                    .copied()
                    .unwrap_or(0.0);
                write!(f, "\t{val:.4}")?;
            }
            writeln!(f)?;
        }
        Ok(())
    }

    pub fn export_optimal_codons(&self, path: &str) -> Result<()> {
        let Some(optimal) = &self.optimal_codons else {
            anyhow::bail!("Optimal codons have not been calculated yet");
        };
        let mut f = fs::File::create(path)?;
        writeln!(f, "# GCUA-rs optimal codons")?;
        writeln!(
            f,
            "# Genetic code: [{}] {}",
            self.genetic_code,
            self.genetic_code_name()
        )?;
        writeln!(f, "AA\tRNA_Codon\tDNA_Codon")?;
        let mut aas: Vec<&String> = optimal.keys().collect();
        aas.sort();
        for aa in aas {
            let rna = &optimal[aa];
            writeln!(f, "{aa}\t{rna}\t{}", dna_codon(rna))?;
        }
        Ok(())
    }
}

fn pct(part: u32, total: u32) -> f64 {
    if total > 0 {
        part as f64 / total as f64 * 100.0
    } else {
        0.0
    }
}
