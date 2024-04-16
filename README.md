# CPBS 7712 Assignment: Contig Assembly

• [Overview](#Overview) <br>
• [Dependencies](#dependencies) <br>
• [Installation](#installation) <br>
• [Usage](#usage) <br>
• [Input Arguments](#input-arguments) <br>
• [Output Files](#output-files) <br>
• [Additional Documentation](#additional-documentation) <br>

## Overview
This is a python-based command line tool and it takes as input the set of all next-generation sequencing reads identified in a sample and an initial query sequence and returns the largest sequence contig that can be constructed from the reads that contains the initial query sequence. <br>

## Dependencies
To run this, you will need to have python 3.11 installed with the following packages

* `argparse`
* `pip`

## Installation
`$ git clone git@github.com:am794/CPBS7712.git` <br>
`$ conda env create --name contigAssembly --file ./envs/environment.yaml`


## Usage
`Python ./src/contig_assembly.py -s ./example/READS.fasta -q ./example/QUERY.fasta -k 15 -sl 7 -o "./example/"`

• Use a subset file (with first 50 reads) to check the results: <br>
`Python ./src/contig_assembly.py -s ./example/subset_READS.fasta -q ./example/QUERY.fasta -k 15 -sl 7 -o "./example/"`

## Input Arguments
This require two input files<br>
-s or --SequencerReads: Sequencer reads in FASTA format<br>
-q or --QuerySequence: Query sequence in FASTA format<br>

Additional arguments:<br>
-k or --kmerLen: K-mer length<br>
-sl or --seedLen: Seed length to create hash index for alignment <br>
-o or --outputPath: path to save output files<br>

## Output files
Two output files:<br>
a) ALLELES.fasta: fasta file of the largest constructed contig (allele) containing the initial query<br>
b) ALLELES.aln: tab-delimited text file describing alignment of sequence reads to contig(s) in ALLELES.fasta in the following format<br>

| sseqid | qseqid | sstart | send | qstart | qend |
| :----: | :----: | :----: | :--: | :----: | :--: |

`• sseqid: name of sequencing read (from READS.fastq.gz)` <br>

`• qseqid: name of contig matched (from ALLELES.fasta)` <br>

`• sstart: starting coordinate in sequencing read sseqid that matches qseq` <br>

`• qseq: ending coordinate in sequencing read sseqid that matches qseq` <br>

`• qstart: starting coordinate in contig that matches sseq` <br>

`• qend: ending coordinate in contig that matches sseq` <br>

## Additional Documentation
https://github.com/am794/CPBS7712/blob/master/Report.pdf