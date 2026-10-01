//! NCBI genetic code (translation table) definitions, mirroring the tables
//! used by the original GCUA Python program.

/*
Gaurav Sablok
gsablok@proton.me
 */

use once_cell::sync::Lazy;
use std::collections::HashMap;

/// The 64 RNA codons in a fixed, stable order (used for matrix-style output).
pub const CODONS: [&str; 64] = [
    "UUU", "UUC", "UUA", "UUG", "UCU", "UCC", "UCA", "UCG", "UAU", "UAC", "UAA", "UAG", "UGU",
    "UGC", "UGA", "UGG", "CUU", "CUC", "CUA", "CUG", "CCU", "CCC", "CCA", "CCG", "CAU", "CAC",
    "CAA", "CAG", "CGU", "CGC", "CGA", "CGG", "AUU", "AUC", "AUA", "AUG", "ACU", "ACC", "ACA",
    "ACG", "AAU", "AAC", "AAA", "AAG", "AGU", "AGC", "AGA", "AGG", "GUU", "GUC", "GUA", "GUG",
    "GCU", "GCC", "GCA", "GCG", "GAU", "GAC", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG",
];

/// Convert an RNA codon (with U) into its DNA form (with T).
pub fn dna_codon(rna: &str) -> String {
    rna.replace('U', "T")
}

/// The Standard / Universal genetic code (NCBI translation table 1).
pub static STANDARD_CODE: Lazy<HashMap<&'static str, &'static str>> = Lazy::new(|| {
    [
        ("UUU", "Phe"),
        ("UUC", "Phe"),
        ("UUA", "Leu"),
        ("UUG", "Leu"),
        ("UCU", "Ser"),
        ("UCC", "Ser"),
        ("UCA", "Ser"),
        ("UCG", "Ser"),
        ("UAU", "Tyr"),
        ("UAC", "Tyr"),
        ("UAA", "STOP"),
        ("UAG", "STOP"),
        ("UGU", "Cys"),
        ("UGC", "Cys"),
        ("UGA", "STOP"),
        ("UGG", "Trp"),
        ("CUU", "Leu"),
        ("CUC", "Leu"),
        ("CUA", "Leu"),
        ("CUG", "Leu"),
        ("CCU", "Pro"),
        ("CCC", "Pro"),
        ("CCA", "Pro"),
        ("CCG", "Pro"),
        ("CAU", "His"),
        ("CAC", "His"),
        ("CAA", "Gln"),
        ("CAG", "Gln"),
        ("CGU", "Arg"),
        ("CGC", "Arg"),
        ("CGA", "Arg"),
        ("CGG", "Arg"),
        ("AUU", "Ile"),
        ("AUC", "Ile"),
        ("AUA", "Ile"),
        ("AUG", "Met"),
        ("ACU", "Thr"),
        ("ACC", "Thr"),
        ("ACA", "Thr"),
        ("ACG", "Thr"),
        ("AAU", "Asn"),
        ("AAC", "Asn"),
        ("AAA", "Lys"),
        ("AAG", "Lys"),
        ("AGU", "Ser"),
        ("AGC", "Ser"),
        ("AGA", "Arg"),
        ("AGG", "Arg"),
        ("GUU", "Val"),
        ("GUC", "Val"),
        ("GUA", "Val"),
        ("GUG", "Val"),
        ("GCU", "Ala"),
        ("GCC", "Ala"),
        ("GCA", "Ala"),
        ("GCG", "Ala"),
        ("GAU", "Asp"),
        ("GAC", "Asp"),
        ("GAA", "Glu"),
        ("GAG", "Glu"),
        ("GGU", "Gly"),
        ("GGC", "Gly"),
        ("GGA", "Gly"),
        ("GGG", "Gly"),
    ]
    .into_iter()
    .collect()
});

fn with_overrides(
    base: &HashMap<&'static str, &'static str>,
    overrides: &[(&'static str, &'static str)],
) -> HashMap<&'static str, &'static str> {
    let mut m = base.clone();
    for (k, v) in overrides {
        m.insert(k, v);
    }
    m
}

/// All supported NCBI genetic codes, keyed by their translation table id.
pub static GENETIC_CODES: Lazy<HashMap<u8, HashMap<&'static str, &'static str>>> =
    Lazy::new(|| {
        let base = &*STANDARD_CODE;
        let mut codes: HashMap<u8, HashMap<&'static str, &'static str>> = HashMap::new();

        codes.insert(1, base.clone()); // Standard
        codes.insert(
            2, // Vertebrate Mitochondrial
            with_overrides(
                base,
                &[
                    ("AUA", "Met"),
                    ("UGA", "Trp"),
                    ("AGA", "STOP"),
                    ("AGG", "STOP"),
                ],
            ),
        );
        codes.insert(
            3, // Yeast Mitochondrial
            with_overrides(
                base,
                &[
                    ("AUA", "Met"),
                    ("UGA", "Trp"),
                    ("CUU", "Thr"),
                    ("CUC", "Thr"),
                    ("CUA", "Thr"),
                    ("CUG", "Thr"),
                ],
            ),
        );
        codes.insert(4, with_overrides(base, &[("UGA", "Trp")])); // Mold/Protozoan/Coelenterate Mito & Mycoplasma
        codes.insert(
            5, // Invertebrate Mitochondrial
            with_overrides(
                base,
                &[
                    ("AUA", "Met"),
                    ("UGA", "Trp"),
                    ("AGA", "Ser"),
                    ("AGG", "Ser"),
                ],
            ),
        );
        codes.insert(6, with_overrides(base, &[("UAA", "Gln"), ("UAG", "Gln")])); // Ciliate/Dasycladacean/Hexamita Nuclear
        codes.insert(
            9, // Echinoderm and Flatworm Mitochondrial
            with_overrides(
                base,
                &[
                    ("AAA", "Asn"),
                    ("AGA", "Ser"),
                    ("AGG", "Ser"),
                    ("UGA", "Trp"),
                    ("AUA", "Ile"),
                ],
            ),
        );
        codes.insert(10, with_overrides(base, &[("UGA", "Cys")])); // Euplotid Nuclear
        codes.insert(11, base.clone()); // Bacterial, Archaeal, Plant Plastid
        codes.insert(12, with_overrides(base, &[("CUG", "Ser")])); // Alternative Yeast Nuclear
        codes.insert(
            13, // Ascidian Mitochondrial
            with_overrides(
                base,
                &[
                    ("AUA", "Met"),
                    ("UGA", "Trp"),
                    ("AGA", "Gly"),
                    ("AGG", "Gly"),
                ],
            ),
        );
        codes.insert(
            14, // Alternative Flatworm Mitochondrial
            with_overrides(
                base,
                &[
                    ("AAA", "Asn"),
                    ("AGA", "Ser"),
                    ("AGG", "Ser"),
                    ("UAA", "Tyr"),
                    ("UGA", "Trp"),
                ],
            ),
        );
        codes.insert(16, with_overrides(base, &[("UAG", "Leu")])); // Chlorophycean Mitochondrial
        codes.insert(
            21, // Trematode Mitochondrial
            with_overrides(
                base,
                &[
                    ("UGA", "Trp"),
                    ("AUA", "Met"),
                    ("AAA", "Asn"),
                    ("AGA", "Ser"),
                    ("AGG", "Ser"),
                ],
            ),
        );
        codes.insert(22, with_overrides(base, &[("UCA", "STOP"), ("UAG", "Leu")])); // Scenedesmus obliquus Mito
        codes.insert(23, with_overrides(base, &[("UUA", "STOP")])); // Thraustochytrium Mitochondrial
        codes.insert(
            24, // Pterobranchia Mitochondrial
            with_overrides(base, &[("AGA", "Ser"), ("AGG", "Lys"), ("UGA", "Trp")]),
        );
        codes.insert(25, with_overrides(base, &[("UGA", "Gly")])); // Candidate Division SR1 and Gracilibacteria
        codes.insert(26, with_overrides(base, &[("CUG", "Ala")])); // Pachysolen tannophilus Nuclear
        codes.insert(
            27,
            with_overrides(base, &[("UAA", "Gln"), ("UAG", "Gln"), ("UGA", "Trp")]),
        ); // Karyorelict Nuclear
        codes.insert(
            28,
            with_overrides(base, &[("UAA", "Gln"), ("UAG", "Gln"), ("UGA", "Trp")]),
        ); // Condylostoma Nuclear
        codes.insert(29, with_overrides(base, &[("UAA", "Tyr"), ("UAG", "Tyr")])); // Mesodinium Nuclear
        codes.insert(30, with_overrides(base, &[("UAA", "Glu"), ("UAG", "Glu")])); // Peritrich Nuclear
        codes.insert(
            31,
            with_overrides(base, &[("UGA", "Trp"), ("UAG", "Glu"), ("UAA", "Glu")]),
        ); // Blastocrithidia Nuclear
        codes.insert(
            33, // Cephalodiscidae Mitochondrial UAA-Tyr
            with_overrides(
                base,
                &[
                    ("UAA", "Tyr"),
                    ("UGA", "Trp"),
                    ("AGA", "Ser"),
                    ("AGG", "Lys"),
                ],
            ),
        );

        codes
    });

/// Human-readable names for each genetic code id.
pub static GENETIC_CODE_NAMES: Lazy<HashMap<u8, &'static str>> = Lazy::new(|| {
    HashMap::from([
        (1, "Standard (Universal)"),
        (2, "Vertebrate Mitochondrial"),
        (3, "Yeast Mitochondrial"),
        (
            4,
            "Mold, Protozoan, Coelenterate Mitochondrial & Mycoplasma/Spiroplasma",
        ),
        (5, "Invertebrate Mitochondrial"),
        (6, "Ciliate, Dasycladacean and Hexamita Nuclear"),
        (9, "Echinoderm and Flatworm Mitochondrial"),
        (10, "Euplotid Nuclear"),
        (11, "Bacterial, Archaeal, and Plant Plastid"),
        (12, "Alternative Yeast Nuclear"),
        (13, "Ascidian Mitochondrial"),
        (14, "Alternative Flatworm Mitochondrial"),
        (16, "Chlorophycean Mitochondrial"),
        (21, "Trematode Mitochondrial"),
        (22, "Scenedesmus obliquus Mitochondrial"),
        (23, "Thraustochytrium Mitochondrial"),
        (24, "Pterobranchia Mitochondrial"),
        (25, "Candidate Division SR1 and Gracilibacteria"),
        (26, "Pachysolen tannophilus Nuclear"),
        (27, "Karyorelict Nuclear"),
        (28, "Condylostoma Nuclear"),
        (29, "Mesodinium Nuclear"),
        (30, "Peritrich Nuclear"),
        (31, "Blastocrithidia Nuclear"),
        (33, "Cephalodiscidae Mitochondrial UAA-Tyr"),
    ])
});

/// Build an amino-acid -> synonymous-codons map for a given genetic code table.
pub fn aa_to_codons(code: &HashMap<&'static str, &'static str>) -> HashMap<String, Vec<String>> {
    let mut m: HashMap<String, Vec<String>> = HashMap::new();
    for codon in CODONS.iter() {
        let aa = code.get(codon).copied().unwrap_or("Xaa");
        m.entry(aa.to_string())
            .or_default()
            .push((*codon).to_string());
    }
    m
}
