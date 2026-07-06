## Data preparation

The statistical analyses and graphical visualizations presented in the publication **"Immune checkpoint inhibitor-associated cardiotoxicity in lung cancer patients assessed by blood transcriptome profiling"** require preprocessing of the raw RNA-seq data.

Before running the analysis scripts `new_run.R`:

1. Download the raw sequencing data from the **GEO** repository.
2. Align the reads to the reference genome using **STAR**.
3. Generate gene-level read counts with **HTSeq**.
4. Place the resulting HTSeq count files in the `htseq/` directory.

Once these steps are completed, the analysis scripts can be executed without further data preprocessing.
