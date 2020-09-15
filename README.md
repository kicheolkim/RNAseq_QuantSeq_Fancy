# QuantSeq 3'mRNAseq data analysis
Code for RNA-Seq data analysis of the paper titled "Aberrant oligodendroglial-vascular interactions disrupt the Blood Brain Barrier triggering CNS inflammation" (Jianqin Niu, Hui-Hsin Tsai, Kimberly Hoi, Nanxin Huang, Guangdan Yu, Kicheol Kim, Sergio Baranzini, Lan Xiao, Jonah R. Chan, Stephen P.J. Fancy). Nature Neuroscience (22, pages 709–718(2019)) https://www.nature.com/articles/s41593-019-0369-4

The sequencing library was prepared according to the manufacturer's instructions (Lexogen).

## Analysis workflow
#### 1. Trim Galore v0.4.4 with CutAdapt v1.14: Trimming illumina adaptor
#### 2. FastQC v0.11.6: QC of fastq file
#### 3. STAR aligner v2.5.0a: sequence alignment
- genome reference: Gencode mouse genome reference (M14), https://www.gencodegenes.org/mouse/release_M14.html
#### 4. Statistical analysis using R and Bioconductor 
- R version 3.5.1 (2018-07-02) and Bioconductor version 3.7 (BiocInstaller 1.30.0)
- Additional packages: ggplot2_3.1.0, RColorBrewer_1.1-2, pheatmap_1.0.10, DESeq2_1.20.0, readr_1.2.1                


### References
- Martin M. (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17(1):10-12.
- Andrews S. (2017). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
- Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29(1):15-21.
- Krueger F. (2017) Trim Galore! Available online at: http://www.bioinformatics.babraham.ac.uk/projects/trim_galore
- R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
- Hadley Wickham, Jim Hester and Romain Francois (2018). readr: Read Rectangular Text Data. R package version 1.2.1. https://CRAN.R-project.org/package=readr
- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15(12):550 (2014)
- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
- Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer
- Raivo Kolde (2018). pheatmap: Pretty Heatmaps. R package version 1.0.10. https://CRAN.R-project.org/package=pheatmap
- Hadley Wickham, Jim Hester and Romain Francois (2018). readr: Read Rectangular Text Data. R package version 1.2.1. https://CRAN.R-project.org/package=readr


-----
- Created by Kicheol Kim, PhD (January 2019. updated)
- Baranzini Lab. (https://github.com/baranzini-lab), Department of Neurology, UCSF
