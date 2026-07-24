use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::mpsc::{self, Receiver, Sender, TryRecvError};

use rayon::spawn;
use tempfile::tempfile;

use super::{
    G4, QuartetBase, ScanLimits, SequenceTopology, chunk_size_for_limits, compute_chunk_overlap,
    consolidate_g4s_with_topology, find_raw_bytes_no_chunking, input::open_input_reader,
    parse_chrom_name, retain_circular_raw_hits, shift_g4,
};

const REVERSE_COMPLEMENT_BLOCK_SIZE: usize = 64 * 1024;

pub struct StreamChromosomeResults {
    pub hits: Vec<G4>,
    pub family_ranges: Vec<(usize, usize)>,
    pub raw_hits: Option<Vec<G4>>,
}

pub fn process_fasta_stream<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_fasta_stream_with_limits_topology(
        path,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_overlap<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    process_fasta_stream_with_limits_overlap_topology(
        path,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_fasta_stream_with_limits_topology(
        path,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_topology<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_topology(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_topology_and_len<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>, usize) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_topology_and_len(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_topology_and_len_with_base<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>, usize) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_topology_and_len_with_base(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        target_base,
        &mut on_chromosome,
    )
}

/// Streams forward and reverse-complement scans for each chromosome.
///
/// The second hit vector passed to `on_chromosome` uses 1-based coordinates in
/// the reverse-complement sequence. Callers that export reference coordinates
/// must project those intervals using the reported chromosome length. The
/// current chromosome is spooled to an automatically removed temporary file.
pub fn process_fasta_stream_bidirectional_with_limits_topology_and_len_with_base<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>, Vec<G4>, usize) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_bidirectional_with_limits_topology_and_len_with_base(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        target_base,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_topology_and_sequence<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, Vec<G4>, Vec<u8>) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_topology_and_sequence(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    process_fasta_stream_with_limits_overlap_topology(
        path,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap_topology<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_overlap_topology(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap_topology_and_len<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults, usize) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_overlap_topology_and_len(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap_topology_and_len_with_base<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults, usize) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_overlap_topology_and_len_with_base(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        target_base,
        &mut on_chromosome,
    )
}

/// Streams forward and reverse-complement scans with raw hits and family ranges.
///
/// The reverse `StreamChromosomeResults` passed to `on_chromosome` remains in
/// reverse-complement coordinate space. The current chromosome is spooled to an
/// automatically removed temporary file before the reverse scan is dispatched.
pub fn process_fasta_stream_bidirectional_with_limits_overlap_topology_and_len_with_base<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults, StreamChromosomeResults, usize) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_bidirectional_with_limits_overlap_topology_and_len_with_base(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        target_base,
        &mut on_chromosome,
    )
}

pub fn process_fasta_stream_with_limits_overlap_topology_and_sequence<F>(
    path: &Path,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    mut on_chromosome: F,
) -> io::Result<usize>
where
    F: FnMut(String, StreamChromosomeResults, Vec<u8>) -> io::Result<()>,
{
    let reader = open_input_reader(path)?;
    process_reader_with_limits_overlap_topology_and_sequence(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        &mut on_chromosome,
    )
}

pub fn process_reader<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_reader_with_limits_topology(
        reader,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_reader_with_overlap<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    process_reader_with_limits_overlap_topology(
        reader,
        min_tetrads,
        min_score,
        ScanLimits::default(),
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_reader_with_limits<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    process_reader_with_limits_topology(
        reader,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_reader_with_limits_topology<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>) -> io::Result<()>,
{
    let mut line = String::new();
    let mut chrom_index = 0usize;
    let mut current: Option<StreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            if let Some(chrom) = current.take() {
                let (name, results) = chrom.finish();
                on_chromosome(name, results)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
            ));
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte.to_ascii_lowercase());
            }
        }
    }

    if let Some(chrom) = current {
        let (name, results) = chrom.finish();
        on_chromosome(name, results)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_with_limits_topology_and_len<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>, usize) -> io::Result<()>,
{
    process_reader_with_limits_topology_and_len_with_base(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        QuartetBase::G,
        on_chromosome,
    )
}

fn process_reader_with_limits_topology_and_len_with_base<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>, usize) -> io::Result<()>,
{
    let mut line = String::new();
    let mut chrom_index = 0usize;
    let mut current: Option<StreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            if let Some(chrom) = current.take() {
                let (name, results, sequence_len) = chrom.finish_with_sequence_len();
                on_chromosome(name, results, sequence_len)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new_with_base(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new_with_base(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            ));
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte.to_ascii_lowercase());
            }
        }
    }

    if let Some(chrom) = current {
        let (name, results, sequence_len) = chrom.finish_with_sequence_len();
        on_chromosome(name, results, sequence_len)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_with_limits_topology_and_sequence<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>, Vec<u8>) -> io::Result<()>,
{
    let mut line = String::new();
    let mut chrom_index = 0usize;
    let mut current: Option<StreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            if let Some(chrom) = current.take() {
                let (name, results, sequence) = chrom.finish_with_sequence();
                on_chromosome(name, results, sequence)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new_with_sequence_capture(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
                true,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new_with_sequence_capture(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
                true,
            ));
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte.to_ascii_lowercase());
            }
        }
    }

    if let Some(chrom) = current {
        let (name, results, sequence) = chrom.finish_with_sequence();
        on_chromosome(name, results, sequence)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

pub fn process_reader_with_limits_overlap<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    process_reader_with_limits_overlap_topology(
        reader,
        min_tetrads,
        min_score,
        limits,
        SequenceTopology::Linear,
        on_chromosome,
    )
}

pub fn process_reader_with_limits_overlap_topology<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults) -> io::Result<()>,
{
    let mut line = String::new();
    let mut chrom_index = 0usize;
    let mut current: Option<StreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            if let Some(chrom) = current.take() {
                let (name, results) = chrom.finish_with_overlap();
                on_chromosome(name, results)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
            ));
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte.to_ascii_lowercase());
            }
        }
    }

    if let Some(chrom) = current {
        let (name, results) = chrom.finish_with_overlap();
        on_chromosome(name, results)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_with_limits_overlap_topology_and_len<R, F>(
    reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults, usize) -> io::Result<()>,
{
    process_reader_with_limits_overlap_topology_and_len_with_base(
        reader,
        min_tetrads,
        min_score,
        limits,
        topology,
        QuartetBase::G,
        on_chromosome,
    )
}

fn process_reader_with_limits_overlap_topology_and_len_with_base<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults, usize) -> io::Result<()>,
{
    let mut line = String::new();
    let mut chrom_index = 0usize;
    let mut current: Option<StreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            if let Some(chrom) = current.take() {
                let (name, results, sequence_len) = chrom.finish_with_overlap_and_sequence_len();
                on_chromosome(name, results, sequence_len)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new_with_base(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new_with_base(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            ));
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte.to_ascii_lowercase());
            }
        }
    }

    if let Some(chrom) = current {
        let (name, results, sequence_len) = chrom.finish_with_overlap_and_sequence_len();
        on_chromosome(name, results, sequence_len)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_with_limits_overlap_topology_and_sequence<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults, Vec<u8>) -> io::Result<()>,
{
    let mut line = String::new();
    let mut chrom_index = 0usize;
    let mut current: Option<StreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            if let Some(chrom) = current.take() {
                let (name, results, sequence) = chrom.finish_with_overlap_and_sequence();
                on_chromosome(name, results, sequence)?;
            }
            chrom_index += 1;
            let name = parse_chrom_name(&line, chrom_index);
            current = Some(StreamChromosome::new_with_sequence_capture(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
                true,
            ));
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(StreamChromosome::new_with_sequence_capture(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
                true,
            ));
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.bytes() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte.to_ascii_lowercase());
            }
        }
    }

    if let Some(chrom) = current {
        let (name, results, sequence) = chrom.finish_with_overlap_and_sequence();
        on_chromosome(name, results, sequence)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_bidirectional_with_limits_topology_and_len_with_base<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, Vec<G4>, Vec<G4>, usize) -> io::Result<()>,
{
    let mut line = Vec::new();
    let mut chrom_index = 0usize;
    let mut current: Option<BidirectionalStreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.starts_with(b">") {
            if let Some(chrom) = current.take() {
                let (name, forward_hits, reverse_hits_rc, sequence_len) = chrom.finish()?;
                on_chromosome(name, forward_hits, reverse_hits_rc, sequence_len)?;
            }
            chrom_index += 1;
            let name = parse_bidirectional_chrom_name(&line, chrom_index)?;
            current = Some(BidirectionalStreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            )?);
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(BidirectionalStreamChromosome::new(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            )?);
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.iter().copied() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte)?;
            }
        }
    }

    if let Some(chrom) = current {
        let (name, forward_hits, reverse_hits_rc, sequence_len) = chrom.finish()?;
        on_chromosome(name, forward_hits, reverse_hits_rc, sequence_len)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn process_reader_bidirectional_with_limits_overlap_topology_and_len_with_base<R, F>(
    mut reader: R,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    on_chromosome: &mut F,
) -> io::Result<usize>
where
    R: BufRead,
    F: FnMut(String, StreamChromosomeResults, StreamChromosomeResults, usize) -> io::Result<()>,
{
    let mut line = Vec::new();
    let mut chrom_index = 0usize;
    let mut current: Option<BidirectionalStreamChromosome> = None;

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.starts_with(b">") {
            if let Some(chrom) = current.take() {
                let (name, forward, reverse_rc, sequence_len) = chrom.finish_with_overlap()?;
                on_chromosome(name, forward, reverse_rc, sequence_len)?;
            }
            chrom_index += 1;
            let name = parse_bidirectional_chrom_name(&line, chrom_index)?;
            current = Some(BidirectionalStreamChromosome::new(
                name,
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            )?);
            continue;
        }
        if current.is_none() {
            chrom_index += 1;
            let fallback = format!("chromosome_{}", chrom_index);
            current = Some(BidirectionalStreamChromosome::new(
                fallback,
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            )?);
        }
        if let Some(chrom) = current.as_mut() {
            for byte in line.iter().copied() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                chrom.push_byte(byte)?;
            }
        }
    }

    if let Some(chrom) = current {
        let (name, forward, reverse_rc, sequence_len) = chrom.finish_with_overlap()?;
        on_chromosome(name, forward, reverse_rc, sequence_len)?;
        Ok(chrom_index.max(1))
    } else {
        Ok(0)
    }
}

fn parse_bidirectional_chrom_name(line: &[u8], chrom_index: usize) -> io::Result<String> {
    let header = std::str::from_utf8(line).map_err(|error| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid UTF-8 in FASTA header {} while preparing bidirectional stream scan: {}",
                chrom_index, error
            ),
        )
    })?;
    Ok(parse_chrom_name(header, chrom_index))
}

struct BidirectionalStreamChromosome {
    name: String,
    forward_scheduler: StreamChunkScheduler,
    spool: BufWriter<File>,
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
}

impl BidirectionalStreamChromosome {
    fn new(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        target_base: QuartetBase,
    ) -> io::Result<Self> {
        let spool_file = tempfile().map_err(|error| {
            spool_io_error(
                error,
                &name,
                "create the temporary reverse-complement spool",
            )
        })?;
        Ok(Self {
            name,
            forward_scheduler: StreamChunkScheduler::new(
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            ),
            spool: BufWriter::new(spool_file),
            min_tetrads,
            min_score,
            limits,
            topology,
            target_base,
        })
    }

    fn push_byte(&mut self, byte: u8) -> io::Result<()> {
        let original_position = self.forward_scheduler.sequence_len() + 1;
        let normalized = normalize_iupac_base(byte, &self.name, original_position)?;
        self.spool.write_all(&[normalized]).map_err(|error| {
            spool_io_error(
                error,
                &self.name,
                &format!(
                    "write original sequence position {} to the reverse-complement spool",
                    original_position
                ),
            )
        })?;
        self.forward_scheduler.push_byte(normalized);
        Ok(())
    }

    fn finish(mut self) -> io::Result<(String, Vec<G4>, Vec<G4>, usize)> {
        let sequence_len = self.forward_scheduler.sequence_len();
        let reverse_scheduler = self.build_reverse_scheduler(sequence_len)?;
        let forward_hits = self.forward_scheduler.finish();
        let reverse_hits_rc = reverse_scheduler.finish();
        Ok((self.name, forward_hits, reverse_hits_rc, sequence_len))
    }

    fn finish_with_overlap(
        mut self,
    ) -> io::Result<(
        String,
        StreamChromosomeResults,
        StreamChromosomeResults,
        usize,
    )> {
        let sequence_len = self.forward_scheduler.sequence_len();
        let reverse_scheduler = self.build_reverse_scheduler(sequence_len)?;
        let (forward_hits, forward_ranges, forward_raw_hits) =
            self.forward_scheduler.finish_with_overlap();
        let (reverse_hits, reverse_ranges, reverse_raw_hits) =
            reverse_scheduler.finish_with_overlap();
        let forward = StreamChromosomeResults {
            hits: forward_hits,
            family_ranges: forward_ranges,
            raw_hits: Some(forward_raw_hits),
        };
        let reverse_rc = StreamChromosomeResults {
            hits: reverse_hits,
            family_ranges: reverse_ranges,
            raw_hits: Some(reverse_raw_hits),
        };
        Ok((self.name, forward, reverse_rc, sequence_len))
    }

    fn build_reverse_scheduler(&mut self, sequence_len: usize) -> io::Result<StreamChunkScheduler> {
        self.spool.flush().map_err(|error| {
            spool_io_error(
                error,
                &self.name,
                "flush the temporary reverse-complement spool",
            )
        })?;
        let mut reverse_scheduler = StreamChunkScheduler::new_with_primary_window_ownership(
            self.min_tetrads,
            self.min_score,
            self.limits,
            self.topology,
            self.target_base,
        );
        feed_reverse_complement_from_spool(
            self.spool.get_mut(),
            sequence_len,
            &self.name,
            &mut reverse_scheduler,
        )?;
        Ok(reverse_scheduler)
    }
}

fn feed_reverse_complement_from_spool(
    spool: &mut File,
    sequence_len: usize,
    chromosome_name: &str,
    scheduler: &mut StreamChunkScheduler,
) -> io::Result<()> {
    let mut block = vec![0u8; REVERSE_COMPLEMENT_BLOCK_SIZE];
    let mut remaining = sequence_len;
    while remaining > 0 {
        let read_len = remaining.min(REVERSE_COMPLEMENT_BLOCK_SIZE);
        let block_start = remaining - read_len;
        let block_start_u64 = u64::try_from(block_start).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "sequence offset {} for chromosome \"{}\" exceeds the temporary spool format",
                    block_start, chromosome_name
                ),
            )
        })?;
        spool
            .seek(SeekFrom::Start(block_start_u64))
            .map_err(|error| {
                spool_io_error(
                    error,
                    chromosome_name,
                    &format!("seek to byte offset {}", block_start),
                )
            })?;
        spool.read_exact(&mut block[..read_len]).map_err(|error| {
            spool_io_error(
                error,
                chromosome_name,
                &format!(
                    "read original sequence positions {} through {}",
                    block_start + 1,
                    block_start + read_len
                ),
            )
        })?;
        for relative_index in (0..read_len).rev() {
            let original_position = block_start + relative_index + 1;
            let complement =
                complement_iupac_base(block[relative_index], chromosome_name, original_position)?;
            scheduler.push_byte(complement);
        }
        remaining = block_start;
    }
    Ok(())
}

fn normalize_iupac_base(
    byte: u8,
    chromosome_name: &str,
    original_position: usize,
) -> io::Result<u8> {
    let normalized = byte.to_ascii_lowercase();
    match normalized {
        b'a' | b't' | b'u' | b'c' | b'g' | b'r' | b'y' | b'k' | b'm' | b'b' | b'v' | b'd'
        | b'h' | b's' | b'w' | b'n' => Ok(normalized),
        _ => Err(invalid_iupac_byte_error(
            byte,
            chromosome_name,
            original_position,
        )),
    }
}

fn complement_iupac_base(
    byte: u8,
    chromosome_name: &str,
    original_position: usize,
) -> io::Result<u8> {
    match byte.to_ascii_lowercase() {
        b'a' => Ok(b't'),
        b't' | b'u' => Ok(b'a'),
        b'c' => Ok(b'g'),
        b'g' => Ok(b'c'),
        b'r' => Ok(b'y'),
        b'y' => Ok(b'r'),
        b'k' => Ok(b'm'),
        b'm' => Ok(b'k'),
        b'b' => Ok(b'v'),
        b'v' => Ok(b'b'),
        b'd' => Ok(b'h'),
        b'h' => Ok(b'd'),
        b's' => Ok(b's'),
        b'w' => Ok(b'w'),
        b'n' => Ok(b'n'),
        _ => Err(invalid_iupac_byte_error(
            byte,
            chromosome_name,
            original_position,
        )),
    }
}

fn invalid_iupac_byte_error(
    byte: u8,
    chromosome_name: &str,
    original_position: usize,
) -> io::Error {
    let escaped_byte: String = char::from(byte).escape_default().collect();
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "invalid IUPAC sequence byte for chromosome \"{}\" at original 1-based position {}: '{}' (0x{:02X}); expected A, T, U, C, G, R, Y, K, M, B, V, D, H, S, W, or N",
            chromosome_name, original_position, escaped_byte, byte
        ),
    )
}

fn spool_io_error(error: io::Error, chromosome_name: &str, action: &str) -> io::Error {
    io::Error::new(
        error.kind(),
        format!(
            "failed to {} for chromosome \"{}\": {}",
            action, chromosome_name, error
        ),
    )
}

struct StreamChromosome {
    name: String,
    scheduler: StreamChunkScheduler,
    captured_sequence: Option<Vec<u8>>,
}

impl StreamChromosome {
    fn new(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
    ) -> Self {
        Self::new_with_base(
            name,
            min_tetrads,
            min_score,
            limits,
            topology,
            QuartetBase::G,
        )
    }

    fn new_with_base(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        target_base: QuartetBase,
    ) -> Self {
        Self::new_with_sequence_capture_and_base(
            name,
            min_tetrads,
            min_score,
            limits,
            topology,
            false,
            target_base,
        )
    }

    fn new_with_sequence_capture(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        capture_sequence: bool,
    ) -> Self {
        Self::new_with_sequence_capture_and_base(
            name,
            min_tetrads,
            min_score,
            limits,
            topology,
            capture_sequence,
            QuartetBase::G,
        )
    }

    fn new_with_sequence_capture_and_base(
        name: String,
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        capture_sequence: bool,
        target_base: QuartetBase,
    ) -> Self {
        Self {
            name,
            scheduler: StreamChunkScheduler::new(
                min_tetrads,
                min_score,
                limits,
                topology,
                target_base,
            ),
            captured_sequence: capture_sequence.then(Vec::new),
        }
    }

    fn push_byte(&mut self, byte: u8) {
        if let Some(sequence) = self.captured_sequence.as_mut() {
            sequence.push(byte);
        }
        self.scheduler.push_byte(byte);
    }

    fn finish(self) -> (String, Vec<G4>) {
        let results = self.scheduler.finish();
        (self.name, results)
    }

    fn finish_with_sequence_len(self) -> (String, Vec<G4>, usize) {
        let sequence_len = self.scheduler.sequence_len();
        let results = self.scheduler.finish();
        (self.name, results, sequence_len)
    }

    fn finish_with_sequence(self) -> (String, Vec<G4>, Vec<u8>) {
        let sequence = self.captured_sequence.unwrap_or_default();
        let results = self.scheduler.finish();
        (self.name, results, sequence)
    }

    fn finish_with_overlap(self) -> (String, StreamChromosomeResults) {
        let (hits, ranges, raw_hits) = self.scheduler.finish_with_overlap();
        (
            self.name,
            StreamChromosomeResults {
                hits,
                family_ranges: ranges,
                raw_hits: Some(raw_hits),
            },
        )
    }

    fn finish_with_overlap_and_sequence_len(self) -> (String, StreamChromosomeResults, usize) {
        let sequence_len = self.scheduler.sequence_len();
        let (name, results) = self.finish_with_overlap();
        (name, results, sequence_len)
    }

    fn finish_with_overlap_and_sequence(self) -> (String, StreamChromosomeResults, Vec<u8>) {
        let sequence = self.captured_sequence.unwrap_or_default();
        let (hits, ranges, raw_hits) = self.scheduler.finish_with_overlap();
        (
            self.name,
            StreamChromosomeResults {
                hits,
                family_ranges: ranges,
                raw_hits: Some(raw_hits),
            },
            sequence,
        )
    }
}

struct StreamChunkScheduler {
    min_tetrads: usize,
    min_score: i32,
    limits: ScanLimits,
    topology: SequenceTopology,
    target_base: QuartetBase,
    chunk_size: usize,
    overlap: usize,
    buffer: VecDeque<u8>,
    offset: usize,
    sequence_len: usize,
    circular_boundary_bp: usize,
    circular_head: VecDeque<u8>,
    circular_tail: VecDeque<u8>,
    tx: Sender<Vec<G4>>,
    rx: Receiver<Vec<G4>>,
    inflight: usize,
    max_inflight: usize,
    completed_hits: Vec<G4>,
    select_chunk_hits: ChunkHitSelector,
}

type FinishParts = (Vec<G4>, Vec<(usize, usize)>, Option<Vec<G4>>);
type ChunkHitSelector = fn(Vec<G4>, usize) -> Vec<G4>;

fn keep_all_chunk_hits(hits: Vec<G4>, _primary_len: usize) -> Vec<G4> {
    hits
}

fn keep_primary_window_hits(hits: Vec<G4>, primary_len: usize) -> Vec<G4> {
    hits.into_iter()
        .filter(|g4| g4.start <= primary_len)
        .collect()
}

impl StreamChunkScheduler {
    fn new(
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        target_base: QuartetBase,
    ) -> Self {
        Self::new_with_hit_selector(
            min_tetrads,
            min_score,
            limits,
            topology,
            target_base,
            keep_all_chunk_hits,
        )
    }

    fn new_with_primary_window_ownership(
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        target_base: QuartetBase,
    ) -> Self {
        Self::new_with_hit_selector(
            min_tetrads,
            min_score,
            limits,
            topology,
            target_base,
            keep_primary_window_hits,
        )
    }

    fn new_with_hit_selector(
        min_tetrads: usize,
        min_score: i32,
        limits: ScanLimits,
        topology: SequenceTopology,
        target_base: QuartetBase,
        select_chunk_hits: ChunkHitSelector,
    ) -> Self {
        let (tx, rx) = mpsc::channel();
        let chunk_size = chunk_size_for_limits(limits);
        let overlap = compute_chunk_overlap(min_tetrads, limits);
        let capacity = chunk_size + overlap;
        let circular_boundary_bp = if topology.is_circular() {
            limits.max_g4_length.saturating_sub(1)
        } else {
            0
        };
        Self {
            min_tetrads,
            min_score,
            limits,
            topology,
            target_base,
            chunk_size,
            overlap,
            buffer: VecDeque::with_capacity(capacity),
            offset: 0,
            sequence_len: 0,
            circular_boundary_bp,
            circular_head: VecDeque::with_capacity(circular_boundary_bp),
            circular_tail: VecDeque::with_capacity(circular_boundary_bp),
            tx,
            rx,
            inflight: 0,
            max_inflight: rayon::current_num_threads().max(1),
            completed_hits: Vec::new(),
            select_chunk_hits,
        }
    }

    fn push_byte(&mut self, byte: u8) {
        self.sequence_len += 1;
        if self.circular_boundary_bp > 0 {
            if self.circular_head.len() < self.circular_boundary_bp {
                self.circular_head.push_back(byte);
            }
            self.circular_tail.push_back(byte);
            if self.circular_tail.len() > self.circular_boundary_bp {
                self.circular_tail.pop_front();
            }
        }
        self.buffer.push_back(byte);
        self.flush_ready_chunks(false);
    }

    fn flush_ready_chunks(&mut self, finishing: bool) {
        let threshold = self.chunk_size + self.overlap;
        while self.buffer.len() >= threshold {
            self.dispatch_chunk(false, threshold);
        }
        if finishing && !self.buffer.is_empty() {
            self.dispatch_chunk(true, self.buffer.len());
        }
    }

    fn dispatch_chunk(&mut self, is_last: bool, window_len: usize) {
        if self.buffer.is_empty() {
            return;
        }
        if self.inflight >= self.max_inflight {
            self.drain_one_completed_chunk();
        }
        let primary_len = if is_last {
            self.buffer.len()
        } else {
            self.chunk_size.min(self.buffer.len())
        };
        if primary_len == 0 {
            return;
        }
        let take = window_len.min(self.buffer.len());
        let mut chunk = Vec::with_capacity(take);
        let (front, back) = self.buffer.as_slices();
        if take <= front.len() {
            chunk.extend_from_slice(&front[..take]);
        } else {
            chunk.extend_from_slice(front);
            let remaining = take - front.len();
            chunk.extend_from_slice(&back[..remaining]);
        }
        // Efficiently remove the primary_len elements from the front.
        self.buffer.drain(..primary_len);
        let offset = self.offset;
        self.offset += primary_len;
        let min_tetrads = self.min_tetrads;
        let min_score = self.min_score;
        let limits = self.limits;
        let target_base = self.target_base;
        let select_chunk_hits = self.select_chunk_hits;
        let tx = self.tx.clone();
        self.inflight += 1;
        spawn(move || {
            // Use the no-chunking variant here: the scheduler already supplied
            // a window (primary + overlap) and we must not re-chunk it.
            let hits =
                find_raw_bytes_no_chunking(chunk, min_tetrads, min_score, limits, target_base);
            let mut hits = select_chunk_hits(hits, primary_len);
            for g4 in &mut hits {
                shift_g4(g4, offset);
            }
            // worker-local dedup is disabled; send raw hits to consolidator
            let _ = tx.send(hits);
        });
    }

    fn drain_one_completed_chunk(&mut self) {
        let mut hits = loop {
            match self.rx.try_recv() {
                Ok(hits) => break hits,
                Err(TryRecvError::Empty) if rayon::current_thread_index().is_some() => {
                    let _ = rayon::yield_now();
                }
                Err(TryRecvError::Empty) => {
                    break self
                        .rx
                        .recv()
                        .expect("stream chunk worker channel disconnected");
                }
                Err(TryRecvError::Disconnected) => {
                    panic!("stream chunk worker channel disconnected");
                }
            }
        };
        self.inflight = self
            .inflight
            .checked_sub(1)
            .expect("completed stream chunk requires an in-flight worker");
        self.completed_hits.append(&mut hits);
    }

    fn finish(self) -> Vec<G4> {
        let (hits, _, _) = self.finish_internal(false);
        hits
    }

    fn finish_with_overlap(self) -> (Vec<G4>, Vec<(usize, usize)>, Vec<G4>) {
        let (hits, ranges, raw) = self.finish_internal(true);
        (
            hits,
            ranges,
            raw.expect("raw hits must be captured when capture_raw is true"),
        )
    }

    fn finish_internal(mut self, capture_raw: bool) -> FinishParts {
        self.flush_ready_chunks(true);
        while self.inflight > 0 {
            self.drain_one_completed_chunk();
        }
        let mut combined = std::mem::take(&mut self.completed_hits);
        if self.topology.is_circular() {
            self.append_wraparound_hits(&mut combined);
            retain_circular_raw_hits(&mut combined, self.sequence_len);
        } else {
            combined.sort_by_key(|a| (a.start, a.end));
        }
        let raw_hits = if capture_raw {
            Some(combined.clone())
        } else {
            None
        };
        let (hits, ranges) =
            consolidate_g4s_with_topology(combined, self.topology, self.sequence_len);
        (hits, ranges, raw_hits)
    }

    fn sequence_len(&self) -> usize {
        self.sequence_len
    }

    fn append_wraparound_hits(&self, combined: &mut Vec<G4>) {
        if self.sequence_len == 0
            || self.circular_boundary_bp == 0
            || self.circular_head.is_empty()
            || self.circular_tail.is_empty()
        {
            return;
        }
        let mut boundary = Vec::with_capacity(self.circular_tail.len() + self.circular_head.len());
        boundary.extend(self.circular_tail.iter().copied());
        boundary.extend(self.circular_head.iter().copied());
        let mut hits = find_raw_bytes_no_chunking(
            boundary,
            self.min_tetrads,
            self.min_score,
            self.limits,
            self.target_base,
        );
        let offset = self.sequence_len.saturating_sub(self.circular_tail.len());
        for g4 in &mut hits {
            shift_g4(g4, offset);
        }
        hits.retain(|g4| g4.end > self.sequence_len);
        combined.extend(hits);
    }
}
