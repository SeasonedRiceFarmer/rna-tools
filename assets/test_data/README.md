# Test data provenance

This tiny dataset exists so the pipeline can be run end-to-end in seconds/minutes on a
laptop, instead of needing a full-size RNA-seq run and a full human genome index.

## Reads

- Source: SRA run `SRR23992075`, part of BioProject
  [PRJNA949611](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA949611)
  ("Whole blood gene expression data and matched white blood cell counts generated
  from human donors diagnosed with a variety of chronic diseases", Case Western
  Reserve University / Grant C. O'Connell lab).
- Real human **whole blood** RNA-seq, paired-end, Illumina NovaSeq 6000, public/open
  access (no dbGaP restriction).
- Subsampled to the **first 20,000 read pairs** of the run with:
  `fastq-dump -X 20000 --split-files SRR23992075.sra`
- Files: `fastq/SRR23992075_1.fastq.gz` (R1), `fastq/SRR23992075_2.fastq.gz` (R2)

## Reference

- `reference/genome_chr21.fa.gz` — human chromosome 21 only (GRCh38, Ensembl release 115)
- `reference/genes_chr21.gtf.gz` — gene annotation, filtered to chromosome 21 rows only

**Only chromosome 21 is included**, not the whole genome, to keep the STAR index small
and fast to build on a laptop. This means most reads in the test data will **not**
align — the highly-expressed whole-blood genes (e.g. hemoglobin `HBB`, `HBA1`, `HBA2`)
are on chromosome 11 and 16, not 21. A low mapping rate here is expected and is not a
bug. This test only checks that the pipeline *mechanism* works (files flow through
every step correctly) — it is not meant to produce biologically meaningful counts.
