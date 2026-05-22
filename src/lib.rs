use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use noodles::sam::alignment::record::cigar::op::Kind;
use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::{Interval, IntervalIndex, IntervalSet};

pub fn peak_counts(
    bam_path: &Path,
    bed_path: &Path,
    output: &mut dyn Write,
    min_mapq: u8,
    workers: NonZero<usize>,
) -> Result<u64> {
    let peaks = load_bed(bed_path)?;
    let set: IntervalSet = peaks.iter().cloned().collect();
    let index = IntervalIndex::build(&set);
    let mut counts: HashMap<(String, u64, u64), u64> = peaks
        .iter()
        .map(|p| ((p.chrom.clone(), p.start, p.end), 0))
        .collect();

    let mut reader = rsomics_bamio::open_with_workers(bam_path, workers)?;
    let header = reader
        .read_header()
        .map_err(|e| RsomicsError::InvalidInput(format!("reading header: {e}")))?;
    let refs: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();

    for result in reader.records() {
        let record =
            result.map_err(|e| RsomicsError::InvalidInput(format!("reading record: {e}")))?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        if record.mapping_quality().is_some_and(|q| q.get() < min_mapq) {
            continue;
        }
        let Some(Ok(tid)) = record.reference_sequence_id() else {
            continue;
        };
        let chrom = &refs[tid];

        // 0-based reference span [start, end) from the alignment start + CIGAR,
        // so a read counts toward a peak it overlaps by >=1bp (matches
        // bedtools multicov), not just one whose start falls inside.
        let Some(Ok(pos)) = record.alignment_start() else {
            continue;
        };
        let start = usize::from(pos) as u64 - 1;
        let ref_len: u64 = record
            .cigar()
            .iter()
            .filter_map(std::result::Result::ok)
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::Match
                        | Kind::Deletion
                        | Kind::Skip
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                )
            })
            .map(|op| op.len() as u64)
            .sum();
        let end = start + ref_len.max(1);

        index.for_each_overlap(chrom, start, end, |peak| {
            *counts
                .entry((chrom.clone(), peak.start, peak.end))
                .or_insert(0) += 1;
        });
    }

    let mut out = BufWriter::new(output);
    let mut total = 0u64;
    for p in &peaks {
        let c = counts
            .get(&(p.chrom.clone(), p.start, p.end))
            .copied()
            .unwrap_or(0);
        writeln!(out, "{}\t{}\t{}\t{c}", p.chrom, p.start, p.end).map_err(RsomicsError::Io)?;
        total += c;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(total)
}

fn load_bed(path: &Path) -> Result<Vec<Interval>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    rsomics_intervals::bed::read(BufReader::new(file))
        .map_err(|e| RsomicsError::InvalidInput(format!("reading BED: {e}")))
}
