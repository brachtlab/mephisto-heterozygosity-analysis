


---
# 1. Measuring heterozygosity as alternative allele frequency 

### 1.1 **Map with `minimap2 onto HAPLOID genome assembly of haplotype 1`**:
  - Haplotype 1 can be arbitrarily determined, and snps will represent the other haplotype.
   - Use the following `LSF` script for mapping reads with `minimap2`: (LSF stands for Load-Sharing Facility and is the job scheduler for our High-Performance Cluster, HPC). 
     ```bash
     #!/bin/bash
     #BSUB -J minimap2
     #BSUB -q normal
     #BSUB -o minimap_log.txt
     #BSUB -e minimap_Error.txt
     #BSUB -n 24

     minimap2 -t 24 -a meph-pri.mmi P3.3.fastq.gz >P3.3.sam
     ```
   - This maps `P3.3.fastq.gz` onto the `meph-pri.mmi` reference and outputs the SAM file `P3.3.sam`.
   - It is critical to map onto a high-quality haploid assembly. The pipeline assumes the snps relative to this assembly are another haplotype. 

### 1.2 **Convert SAM to BAM with `samtools view`**:
   - After mapping, convert the `.sam` file to a `.bam` file using `samtools view`. Example command:
     ```bash
     samtools view -@ 20 -S -b P3.3.sam > P3.3.bam
     ```
   - This converts the `P3.3.sam` file to the binary `P3.3.bam` file.
     
### 1.3 **Sort and Index the BAM file**:
   ```bash
     samtools sort -@ 20 -o P3.3_sorted.bam P3.3.bam
     samtolls index P3.3_sorted.bam
   ```

### 1.4 **Run bcftools to call variants on the HPC**:
   - Use the following `LSF` script to run `bcftools` (LSF stands for Load-Sharing Facility and is the job scheduler for our High-Performance Cluster, HPC). 
```bash
     #!/bin/bash
     #BSUB -J bcftools
     #BSUB -q normal
     #BSUB -o bcftools_Log.txt
     #BSUB -e bcftools_Error.txt
     #BSUB -n 48

     bcftools mpileup -f mephisto_alpha_renamed_polish.fasta_primary.fasta P3.3_sorted.bam | bcftools call -mv -Ov -o P3.3.vcf
 ```
   - submit this file on the HPC's LSF sytem using the bsub command:
     ```bash
        bsub<bcftools.lsf
      ```
   - This script generates the variant calls in the VCF file `P3.3.vcf`.

### 1.5 **Convert vcf file to tab-delimited snp-only file for plotting**:
- **step 1.5.1** : Use the Python script parseVCF-freq2.py to process the VCF file into a tab-delimited text file. Input is FILENAME.vcf, Output is FILENAME.vcf_columns2.txt
   ```
     ./parseVCF-freq2.py P3.3.vcf
   ```

- **Step 1.5.2**: **Filter SNPs**: Filter SNPs using `filter-text-files.py`.
   - Run `filter-text-files.py` to ensure only SNPs are kept. Input is FILENAME.vcf_columns2.txt, Output is FILENAME.vcf_columns2.txt_snps_only.txt
     ```bash
     ./filter-text-files.py P3.3.vcf_columns2.txt
     ```
### 1.6 **Plot using KaroploteR package**:
- **Step 1.6.1**: **get LOH regions**: Use get-LOH.py to process the FILENAME.vcf_columns2.txt file into FILENAME.vcf_columns2.txt_LOH.txt
   ```
     ./get-LOH.py P3.3.vcf_columns2.txt_snps_only.txt
   ```

- **Step 1.6.2**: **Plot using KaryoploteR packages**:

The text file FILENAME.vcf_columns2.txt_snps_only.txt is the input for plotting karyotypes. Follow the tutorial given at: https://bernatgel.github.io/karyoploter_tutorial/. Plot using kpPoints() command and contig, position, and alternative allele frequency data (columns 1, 2, and 7). For Rmarkdown, please go here: http://htmlpreview.github.io/?https://github.com/brachtlab/KaroploteR-heterozygosity/blob/main/Karyoplots-F1-A1.html

### 1.7 **Plot using Scaptterplot method**:
- **Step 1.7.1**: **Plot Scatterplots, part 1**: Compare the filtered files using `compare-text-files.py` **Compare Text Files**:
   - Use `compare-text-files.py` to compare filtered SNP files. Input is FILENAME1.vcf_columns2.txt_snps_only.txt and FILENAME2.vcf_columns2.txt_snps_only.txt
     ```bash
     ./compare-text-files.py FILENAME1.vcf_columns2.txt_snps_only.txt FILENAME2.vcf_columns2.txt_snps_only.txt
     ```
- **Step 1.7.2**: **Plot Scatterplots, part 2**: Using ggplot2 in R to visualize 'ALT-Fraction Alternative Variant: Parent vs. Child' figures, using the geom_hex() command within ggplot2.

  ```
  a3p <- ggplot(a3   , aes(x = V7, y = V17)) +
   geom_hex(bins = 43) + # Increased the number of bins for smaller hexagons
   scale_fill_continuous(type = "viridis") + # Use a continuous color scale
   theme_minimal(base_size = 15) +
   labs(
     title = "Fraction alternative variant parent vs child",
     x = "Parent (F1)",
     y = "Child (A3)",
     fill = "Count"
   ) +
   theme(
     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
     axis.title = element_text(size = 20, face = "bold"),
     axis.text = element_text(size = 20, face = "bold"),
     panel.grid.major = element_line(color = "grey80"),
     panel.grid.minor = element_line(color = "grey90"),
     panel.background = element_rect(fill = "white")
   ) +
   coord_fixed(ratio = 1)
  ```

  





