//! Minimal FASTA reader with duplicate-ID handling, mirroring the behaviour
//! of the Python GCUA loader (no external bioinformatics crate required).

/*
Gaurav Sablok
gsablok@proton.me
 */

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;

pub struct FastaRecord {
    pub id: String,
    pub sequence: String,
}

pub fn read_fasta(path: &str) -> Result<Vec<FastaRecord>> {
    let content =
        fs::read_to_string(path).with_context(|| format!("Failed to read FASTA file: {path}"))?;

    let mut records = Vec::new();
    let mut id_counts: HashMap<String, usize> = HashMap::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for raw_line in content.lines() {
        let line = raw_line.trim_end();
        if line.starts_with('>') {
            if let Some(id) = current_id.take() {
                records.push(FastaRecord {
                    id,
                    sequence: std::mem::take(&mut current_seq),
                });
            }
            let header = line[1..].trim();
            let base_id = header
                .split_whitespace()
                .next()
                .unwrap_or("unnamed")
                .to_string();
            let count = id_counts.entry(base_id.clone()).or_insert(0);
            let id = if *count > 0 {
                format!("{base_id}_{count}")
            } else {
                base_id.clone()
            };
            *count += 1;
            current_id = Some(id);
        } else if current_id.is_some() {
            current_seq.push_str(line.trim());
        }
    }

    if let Some(id) = current_id {
        records.push(FastaRecord {
            id,
            sequence: current_seq,
        });
    }

    Ok(records)
}
