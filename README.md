


---
## Measuring Allele Frequency Between Parent and Progeny 

1. **Map with `minimap2`**:
   - Use the following `LSF` script for mapping reads with `minimap2`:
     ```bash
     #!/bin/bash
     #BSUB -J minimap2
     #BSUB -q long
     #BSUB -o minimap5_log.txt
     #BSUB -e minimap5_Error.txt
     #BSUB -n 24

     minimap2 -t 24 -a meph-pri.mmi P3.3-omega.fastq.gz >P3.3.sam
     ```
   - This maps `P3.3-omega.fastq.gz` onto the `meph-pri.mmi` reference and outputs the SAM file `P3.3.sam`.

2. **Convert SAM to BAM with `samtools view`**:
   - After mapping, convert the `.sam` file to a `.bam` file using `samtools view`. Example command:
     ```bash
     samtools view -@ 20 -S -b P3.3.sam > P3.3.bam
     ```
   - This converts the `P3.3.sam` file to the binary `P3.3.bam` file.
     
3. **Sort and index BAM **:
   ```bash
     samtools sort -@ 20 -0 P3.3_sorted.bam P3.3.bam
     samtolls index P3.3_sorted.bam
     ```

4. **Run `bcftools`**:
   - Use the following `LSF` script to run `bcftools`:
     ```bash
     #!/bin/bash
     #BSUB -J bcftools_JB7
     #BSUB -q normal
     #BSUB -o bcftools_Log-7.txt
     #BSUB -e bcftools_Error-7.txt
     #BSUB -n 48

     /home/jbracht/bcftools-1.20/bcftools mpileup -f mephisto_alpha_renamed_polish.fasta_primary.fasta P3.3_sorted.bam | bcftools call -mv -Ov -o P3.3.vcf
     ```
   - This script generates the variant calls in the VCF file `P3.3.vcf`.

5. **Call Variants**:
   - After generating the VCF file, process the variants as:

    - **Step 4.1**: **Variant Calling using `bcftools`**:Use parseVCF-freq2.py to process the VCF file, requiring at least two reads for both the alternate and reference alleles. 
   Use `parseVCF-freq2.py`to process the VCF file.
   ```bash
   bcftools mpileup --threads 20 -f ../mephisto_alpha_renamed_polish.fasta_primary.fasta P3.3-new2.bam | bcftools call --threads 20 -mv -Ov -o P3.3-new2-calls.vcf
   ```
     ```bash
 
     python parseVCF-freq2.py
     ```

- **Step 5.2**: **Filter SNPs**: Filter SNPs using `filter-text-files.py`.
   - Run `filter-text-files.py` to ensure only SNPs are kept:
     ```bash
     python filter-text-files.py
     ```

- **Step 5.3**: Compare the filtered files using `compare-text-files.py` **Compare Text Files**:
   - Use `compare-text-files.py` to compare filtered SNP files:
     ```bash
     python compare-text-files.py
     ```
- **Step 5.4**: Using ggplot2 in R to visualize 'ALT-Fraction Alternative Variant: Parent vs. Child' figures.
     
# Calling-Recombinant-Reads

This pipeline performs variant calling and recombination detection by mapping reads to a genome, filtering SNPs, and analyzing recombination events from forward and reverse genome alignments. Below are the steps and scripts used in the process:

### 1. Map Reads to Genome:
- Map P3, P3.1, and P3.3 reads onto the genome, generating the following SAM files:
  - `P3.sam`
  - `P3.1.sam`
  - `P3.3.sam`

### 2. Generate VCF File:
- Call P3 variants by running `parseVCF-highConfidenceSNPs.py`, ensuring at least 5 reads support both reference and alternative allele calls. This generates the P3 VCF file.

### 3. Filter SNPs:
- Use `filter-text-files.py` to filter the SNPs.

### 4. Add Genomic Context:
- Run `add-context_fixed.py` (requires the genome file as input) to add genomic context to the SNP file.

### 5. Detect Forward Recombination Events:
- Run `find-recombination.py` using the forward P3 variant calls file:  
  `P3-omega-zorro.vcf_high-conf-snps.txt_snps_only.txt_contextFIXED.txt`  
  This generates the file:  
  `omega_reads_calls2.txt`

### 6. Reverse Complement Mapping:
- Map P3, P3.1, and P3.3 onto the reverse complement of the genome (`rcomega`).
- Call P3 variants as above, filter using `filter-text-files.py`, and add context with `add-context_fixed.py` to generate:  
  `P3-rcomega-zorro.vcf_high-conf-snps.txt_snps_only.txt_contextFIXED.txt`

### 7. Detect Reverse Complement Recombination Events:
- Run `rc-find-recombination.py` using the reverse complement SNPs file. This generates:  
  `rcomega_reads_calls2.txt`

### 8. Reverse Complement Coordinates:
- Run `reverse-complement-coordinates-of-call-file.py` (requires reverse complement genome input) to convert the coordinates back to the original genome's direction. This generates:  
  `reads_calls2.txt_rc.txt`

### 9. Merge Forward and Reverse Call Files:
- Combine forward and reverse complement calls using `cat`:
  ```bash
  cat [forward, omega] [rcomega, reverse complemented] > combined-calls.txt
  ```

### 10. Sorting:
- Open `combined-calls.txt` in Excel.
- Sort by `contig` and `position`.
- Save as a tab-delimited text file.

### 11. Analyze Recombination Clusters:
- Run `find-clusters.py` to analyze the read number cutoff:
  ```bash
  find-clusters.py <read_cutoff>
  ```

### 12. Compare Cluster Analyses:
- Run `find-clusters-compare.py` to compare two cluster analyses:
  ```bash
  find-clusters-compare.py <file1> <read_cutoff_1> <file2> <read_cutoff_2>
  You can add these commands under a new section in your `README.md` file or another documentation file to clearly explain how to run them as part of the workflow. Here’s an example of how you can format it:





