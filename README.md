


---
## 1. Measuring Allele Frequency Between Parent and Progeny 

1.1 **Map with `minimap2 onto HAPLOID genome assembly`**:
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
   - It is critical to map onto a high-quality haploid assembly. The pipeline assumes the snps relative to this assembly are the other genotype. 

1.2 **Convert SAM to BAM with `samtools view`**:
   - After mapping, convert the `.sam` file to a `.bam` file using `samtools view`. Example command:
     ```bash
     samtools view -@ 20 -S -b P3.3.sam > P3.3.bam
     ```
   - This converts the `P3.3.sam` file to the binary `P3.3.bam` file.
     
1.3 **Sort and Index the BAM file**:
   ```bash
     samtools sort -@ 20 -o P3.3_sorted.bam P3.3.bam
     samtolls index P3.3_sorted.bam
   ```

1.4 **Run `bcftools`**:
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

1.5 **Call Variants**:
   - After generating the VCF file, process the variants:

 ** step 1.5.1** : Use parseVCF-freq2.py to process the VCF file into a tab-delimited text file. Input is FILENAME.vcf, Output is FILENAME.vcf_columns2.txt
   ```
     ./parseVCF-freq2.py P3.3.vcf
     ```

- **Step 1.5.2**: **Filter SNPs**: Filter SNPs using `filter-text-files.py`.
   - Run `filter-text-files.py` to ensure only SNPs are kept. Input is FILENAME.vcf_columns2.txt, Output is FILENAME.vcf_columns2.txt_snps_only.txt
     ```bash
     ./filter-text-files.py P3.3.vcf_columns2.txt
     ```

- **Step 1.5.3**: Compare the filtered files using `compare-text-files.py` **Compare Text Files**:
   - Use `compare-text-files.py` to compare filtered SNP files. Input is FILENAME1.vcf_columns2.txt_snps_only.txt and FILENAME2.vcf_columns2.txt_snps_only.txt
     ```bash
     ./compare-text-files.py FILENAME1.vcf_columns2.txt_snps_only.txt FILENAME2.vcf_columns2.txt_snps_only.txt
     ```
- **Step 1.5.4**: Using ggplot2 in R to visualize 'ALT-Fraction Alternative Variant: Parent vs. Child' figures.
     
# Calling-Recombinant-Reads

This pipeline performs variant calling and recombination detection by mapping reads to a genome, filtering SNPs, and analyzing recombination events from forward and reverse genome alignments. Below are the steps and scripts used in the process:

### 1. Map Reads to Genome:
- If you haven't already, follow step 1 above to map reads reads onto the (haploid) genome, generating the following SAM files:
  - `P3.sam`
  - `P3.1.sam`
  - `P3.3.sam`
-Then, follow steps 2-4 above, to generate .vcf files for the next step.

### 2. Identify the snps in haplotype 2:
- Call P3 variants by running `parseVCF-highConfidenceSNPs.py`, ensuring at least 5 reads support both reference and alternative allele calls. Input is FILENAME.vcf and output is FILENAME.vcf_high-conf-snps.txt.
 ```bash
     ./parseVCF-highConfidenceSNPs.py P3.3.vcf
 ```
### 3. Filter SNPs:
- Use `filter-text-files.py` to filter the SNPs. Input is FILENAME.vcf_high-conf-snps.txt and output is FILENAME.vcf_high-conf-snps.txt_snps_only.txt
 ```bash
     ./filter-text-files.py P3.3.vcf_high-conf-snps.txt
 ```

### 4. Add Genomic Context:
- Run `add-context_fixed.py` to add genomic context to the SNP file. (Context is the previous 7 bp and 12 bp relative to the snp. Both are used for identifying snps within the reads using find-recombination.py.) Input is FILENAME.vcf_high-conf-snps.txt_snps_only.txt and the genome.fasta, and output is FILENAME.vcf_high-conf-snps.txt_snps_only.txt_contextFIXED.txt
```bash
     ./add-context_fixed.py P3.3.vcf_high-conf-snps.txt_snps_only.txt genome.fasta
 ```
### 5. Detect Forward Recombination Events:
- Run `find-recombination.py` which requires the SAM file (containing the reads and their mapping positions) along with the snp file (FILENAME.vcf_high-conf-snps.txt_snps_only.txt_contextFIXED.txt). The output is SAMFILE.sam_reads_calls2.txt.    
```bash
     ./add-context_fixed.py P3.3.vcf_high-conf-snps.txt_snps_only.txt genome.fasta
 ```
The output file has 10 columns and one row per read (if the read is long enough to be scored). The recombination-flag is the key because it is set if the read is thought to be recombinant. Error-counts is the number of snps that did not match either a reference or alt basepair, yet the snp site was present in the read (based on the context search). These are presumed errors owing to the fairly high error rate of the long reads.The 'call-list' shows a list of each snp call within the read, whether reference, 'ref' or alternate, 'alt'.  

  

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





