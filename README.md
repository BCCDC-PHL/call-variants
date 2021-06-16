# Call Variants Pipeline

This pipeline takes a directory containing paired-end illumina sequence data, plus a reference sequence,
aligns the reads against the reference and calls variants. It is based directly off of the excellent
[Snippy](https://github.com/tseemann/snippy) by Torsten Seemann ([@tseemann](https://github.com/tseemann)).

## Usage

```
nextflow run BCCDC-PHL/call-variants \
  --ref <ref.fa> \
  --gff <genes.gff> \
  --fastq_input </path/to/fastq> \
  --outdir <outdir>
```
