//! SVG plot generation for GCUA-rs, built on the `kuva` scientific
//! plotting crate. Mirrors (a subset of) the Python version's Plotly
//! visualizations: GC vs GC3, Wright's ENC-vs-GC3s plot, a CAI distribution
//! histogram, and an RSCU heatmap.

/*
Gaurav Sablok
gsablok@proton.me
 */

use crate::data::GcuaData;
use crate::genetic_code::CODONS;
use anyhow::Result;
use kuva::backend::svg::SvgBackend;
use kuva::plot::{Heatmap, Histogram, LinePlot, LineStyle, ScatterPlot};
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use kuva::render::render::render_multiple;
use std::fs;

fn svg(plots: Vec<Plot>, layout: Layout) -> String {
    SvgBackend.render_scene(&render_multiple(plots, layout))
}

/// Ordinary least-squares fit of y = slope * x + intercept.
fn linear_regression(points: &[(f64, f64)]) -> Option<(f64, f64)> {
    let n = points.len() as f64;
    if n < 2.0 {
        return None;
    }
    let mean_x = points.iter().map(|p| p.0).sum::<f64>() / n;
    let mean_y = points.iter().map(|p| p.1).sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for &(x, y) in points {
        num += (x - mean_x) * (y - mean_y);
        den += (x - mean_x) * (x - mean_x);
    }
    if den == 0.0 {
        return None;
    }
    let slope = num / den;
    let intercept = mean_y - slope * mean_x;
    Some((slope, intercept))
}

/// GC vs GC3 scatter plot with an overlaid OLS regression line
/// (mirrors the Python `GCContentVisualization`).
pub fn gc_vs_gc3(data: &GcuaData) -> Option<String> {
    let points: Vec<(f64, f64)> = data
        .gene_names
        .iter()
        .filter_map(|g| data.base_composition.get(g))
        .map(|bc| (bc.gc, bc.gc3))
        .collect();
    if points.len() < 2 {
        return None;
    }

    let scatter = ScatterPlot::new()
        .with_data(points.clone())
        .with_color("steelblue")
        .with_legend("genes");

    let mut plots = vec![Plot::Scatter(scatter)];

    if let Some((slope, intercept)) = linear_regression(&points) {
        let min_x = points.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
        let max_x = points.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
        let line = LinePlot::new()
            .with_data(vec![
                (min_x, slope * min_x + intercept),
                (max_x, slope * max_x + intercept),
            ])
            .with_color("crimson")
            .with_stroke_width(2.0)
            .with_line_style(LineStyle::Dashed)
            .with_legend(format!("y = {slope:.3}x + {intercept:.3}"));
        plots.push(Plot::Line(line));
    }

    let layout = Layout::auto_from_plots(&plots)
        .with_title("GC3 vs GC Content")
        .with_x_label("GC content (%)")
        .with_y_label("GC3 content (%)");

    Some(svg(plots, layout))
}

/// Wright's plot: observed ENC vs GC3s, with the theoretical no-selection
/// curve ENC = 2 + s + 29 / (s^2 + (1-s)^2) overlaid.
pub fn enc_wright_plot(data: &GcuaData) -> Option<String> {
    let enc = data.enc.as_ref()?;

    let points: Vec<(f64, f64)> = data
        .gene_names
        .iter()
        .filter_map(|g| {
            let bc = data.base_composition.get(g)?;
            let e = enc.get(g)?;
            Some((bc.gc3s, *e))
        })
        .collect();
    if points.is_empty() {
        return None;
    }

    let curve: Vec<(f64, f64)> = (0..=100)
        .map(|i| {
            let gc3s = i as f64;
            let s = gc3s / 100.0;
            let denom = s * s + (1.0 - s) * (1.0 - s);
            let expected = if denom > 0.0 {
                (2.0 + s + 29.0 / denom).min(61.0)
            } else {
                61.0
            };
            (gc3s, expected)
        })
        .collect();

    let scatter = ScatterPlot::new()
        .with_data(points)
        .with_color("steelblue")
        .with_legend("Observed ENC");
    let line = LinePlot::new()
        .with_data(curve)
        .with_color("crimson")
        .with_stroke_width(2.0)
        .with_legend("Expected (no selection)");

    let plots = vec![Plot::Scatter(scatter), Plot::Line(line)];
    let layout = Layout::auto_from_plots(&plots)
        .with_title("ENC vs GC3s (Wright's Plot)")
        .with_x_label("GC3s (%)")
        .with_y_label("ENC");

    Some(svg(plots, layout))
}

/// Histogram of CAI values across all genes.
pub fn cai_distribution(data: &GcuaData) -> Option<String> {
    let cai = data.cai.as_ref()?;
    let values: Vec<f64> = data
        .gene_names
        .iter()
        .filter_map(|g| cai.get(g).copied())
        .collect();
    if values.is_empty() {
        return None;
    }

    let hist = Histogram::new()
        .with_data(values)
        .with_bins(20)
        .with_color("steelblue");

    let plots = vec![Plot::Histogram(hist)];
    let layout = Layout::auto_from_plots(&plots)
        .with_title("Distribution of CAI Values")
        .with_x_label("CAI")
        .with_y_label("Frequency");

    Some(svg(plots, layout))
}

/// RSCU heatmap. For small gene sets, renders genes x codons directly;
/// for large ones, renders a 2-row (mean / std-dev) summary instead
/// (mirrors the Python version's large-dataset fallback).
pub fn rscu_heatmap(data: &GcuaData) -> Option<String> {
    if data.rscu.is_empty() {
        return None;
    }

    let non_stop: Vec<&str> = CODONS
        .iter()
        .filter(|c| data.code.get(*c).copied() != Some("STOP"))
        .copied()
        .collect();

    const MAX_GENES_FOR_FULL_HEATMAP: usize = 60;

    let (matrix, row_labels): (Vec<Vec<f64>>, Vec<String>) =
        if data.gene_names.len() <= MAX_GENES_FOR_FULL_HEATMAP {
            let matrix: Vec<Vec<f64>> = data
                .gene_names
                .iter()
                .map(|g| {
                    let gene_rscu = data.rscu.get(g);
                    non_stop
                        .iter()
                        .map(|c| gene_rscu.and_then(|m| m.get(*c)).copied().unwrap_or(0.0))
                        .collect()
                })
                .collect();
            (matrix, data.gene_names.clone())
        } else {
            let n = data.gene_names.len() as f64;
            let mut means = vec![0.0f64; non_stop.len()];
            let mut sq = vec![0.0f64; non_stop.len()];
            for g in &data.gene_names {
                if let Some(gene_rscu) = data.rscu.get(g) {
                    for (i, c) in non_stop.iter().enumerate() {
                        let v = gene_rscu.get(*c).copied().unwrap_or(0.0);
                        means[i] += v;
                        sq[i] += v * v;
                    }
                }
            }
            for m in means.iter_mut() {
                *m /= n;
            }
            let stds: Vec<f64> = sq
                .iter()
                .zip(means.iter())
                .map(|(&s, &m)| ((s / n) - m * m).max(0.0).sqrt())
                .collect();
            (
                vec![means, stds],
                vec!["Mean RSCU".to_string(), "Std Dev".to_string()],
            )
        };

    let col_labels: Vec<String> = non_stop.iter().map(|s| s.to_string()).collect();

    let heatmap = Heatmap::new()
        .with_data(matrix)
        .with_labels(row_labels.clone(), col_labels.clone());

    let plots = vec![Plot::Heatmap(heatmap)];
    let title = if row_labels.len() == 2 {
        format!("RSCU Summary ({} genes)", data.gene_names.len())
    } else {
        "RSCU Heatmap".to_string()
    };
    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_x_categories(col_labels)
        .with_y_categories(row_labels);

    Some(svg(plots, layout))
}

/// Generate all available plots and write them as SVG files into `dir`.
/// Returns the list of files actually written (some plots are skipped if
/// their prerequisite metric hasn't been computed or there isn't enough
/// data).
pub fn write_all_plots(data: &GcuaData, dir: &std::path::Path) -> Result<Vec<std::path::PathBuf>> {
    fs::create_dir_all(dir)?;
    let mut written = Vec::new();

    let jobs: Vec<(&str, Option<String>)> = vec![
        ("gc_vs_gc3.svg", gc_vs_gc3(data)),
        ("enc_wright_plot.svg", enc_wright_plot(data)),
        ("cai_distribution.svg", cai_distribution(data)),
        ("rscu_heatmap.svg", rscu_heatmap(data)),
    ];

    for (filename, maybe_svg) in jobs {
        if let Some(svg_content) = maybe_svg {
            let path = dir.join(filename);
            fs::write(&path, svg_content)?;
            written.push(path);
        }
    }

    Ok(written)
}
