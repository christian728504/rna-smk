use clap::{value_parser, Arg, Command};
use std::{
    error::Error, path::PathBuf
};
use polars::prelude::*;
use plotters::prelude::*;
use std::collections::HashMap;
use ndarray::{Array1, Axis, stack};
use ndarray_stats::CorrelationExt;
use std::fs;

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    values[values.len() / 2]
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut qc_metrics = HashMap::new();

    let matches = Command::new("madqc")
        .about("MAD QC CLI")
        .about("Calculates the MAD of log ratios given two RSEM gene quantification tsv files")
            .arg(
                Arg::new("tsv1")
                    .long("tsv1")
                    .help("Input gene quantification tsv from rsem")
                    .value_parser(value_parser!(PathBuf))
                    .required(true)
                    .index(1)
            )
            .arg(
                Arg::new("tsv2")
                    .long("tsv2")
                    .help("Other input gene quantification tsv from rsem")
                    .value_parser(value_parser!(PathBuf))
                    .required(true)
                    .index(2)
            ).get_matches();
    
    let tsv1 = matches.get_one::<PathBuf>("tsv1").ok_or("missing tsv1")?.to_owned();
    let tsv2 = matches.get_one::<PathBuf>("tsv2").ok_or("missing tsv2")?.to_owned();

    let df1 = CsvReadOptions::default()
        .with_has_header(true)
        .with_infer_schema_length(None)
        .with_parse_options(CsvParseOptions::default().with_separator(b'\t'))
        .try_into_reader_with_file_path(Some(tsv1))?
        .finish()?;

     let df2 = CsvReadOptions::default()
        .with_has_header(true)
        .with_infer_schema_length(None)
        .with_parse_options(CsvParseOptions::default().with_separator(b'\t'))
        .try_into_reader_with_file_path(Some(tsv2))?
        .finish()?;

    let df = df1.inner_join(&df2, ["gene_id"], ["gene_id"])?.select(["gene_id", "TPM", "FPKM", "TPM_right", "FPKM_right"])?;
    let mask = col("FPKM").eq(0.0).or(col("FPKM_right").eq(0.0)).not();
    let df = df.lazy().filter(mask).collect()?;
    
    let expr1 = df.column("FPKM")?.f64()?;
    let expr1_vec: Vec<f64> = expr1.into_no_null_iter().collect();
    let logvec1 = Array1::from(expr1_vec).mapv(|x| x.log2());

    let expr2 = df.column("FPKM_right")?.f64()?;
    let expr2_vec: Vec<f64> = expr2.into_no_null_iter().collect();
    let logvec2 = Array1::from(expr2_vec).mapv(|x| x.log2());

    let a= (&logvec1 + &logvec2) / 2.0;
    let m = &logvec1 - &logvec2;
    let a_mask = a.mapv(|x| x > 0.0);
    let m_abs = m.mapv(|x| x.abs());
    let mut m_filtered: Vec<f64> = m_abs.iter()
        .zip(&a_mask)
        .filter_map(|(&val, &keep)| if keep { Some(val) } else { None })
        .collect();
    let mad = median(&mut m_filtered) * 1.4826;

    let m_filtered: Vec<f64> = m.iter()
        .zip(&a_mask)
        .filter_map(|(&val, &keep)| if keep { Some(val) } else { None })
        .collect();
    let m_filtered_borrow = m_filtered.as_slice();
    let m_power_2: Vec<f64> = m_filtered_borrow.iter().map(|x| x.powf(2.0)).collect();
    let m_mean = m_power_2.iter().sum::<f64>() / m_power_2.len() as f64;
    let sd_of_log_ratios = m_mean.sqrt();

    let logvec1_filtered: Vec<f64> = logvec1.iter()
        .zip(&a_mask)
        .filter_map(|(&val, &keep)| if keep { Some(val) } else { None })
        .collect();
    let logvec2_filtered: Vec<f64> = logvec2.iter()
        .zip(&a_mask)
        .filter_map(|(&val, &keep)| if keep { Some(val) } else { None })
        .collect();
    
    let logarr1 = Array1::from(logvec1_filtered);
    let logarr2 = Array1::from(logvec2_filtered);
    let stacked = stack![Axis(0), logarr1, logarr2];
    let corr = stacked.pearson_correlation()?[[0, 1]];

    let root = BitMapBackend::new("mad_qc_plot.png", (1400, 1400)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("MAD QC Plot", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-10f64..20f64, -10f64..10f64)?;

    chart
        .configure_mesh()
        .x_desc("A")
        .y_desc("M")
        .draw()?;

    chart.draw_series(
        a.iter().zip(m.iter())
            .map(|(&x, &y)| Circle::new((x, y), 3, BLACK.filled()))
    )?;
    root.present()?;

    println!("MAD of log ratios: {mad}");
    println!("Standard deviation of log ratios: {sd_of_log_ratios}");
    println!("Pearson correlation: {corr}");

    qc_metrics.insert("MAD".to_string(), mad);
    qc_metrics.insert("Standard_deviation".to_string(), sd_of_log_ratios);
    qc_metrics.insert("Pearson_correlation".to_string(), corr);

    let qc_metrics_json = serde_json::to_string_pretty(&qc_metrics)?;
    fs::write("mad_qc.json", qc_metrics_json)?;

    Ok(())
}
