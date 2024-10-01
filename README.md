# Recombinant-reads-analyzer
Here’s a refined version of your README file for the Calling-Recombinant-Reads project:

---

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

### 13. Merging and Calling Variants with `samtools` and `bcftools`

To merge multiple `.bam` files into one and perform variant calling, follow these steps:

1. **Merge `.bam` files using `samtools`:**
   ```bash
   samtools merge -o P0_merged.bam PASS/MAPPED/*.fastq/*.bam
   ```

   **Note**: If you encounter an error with too many files being open, use the following command to increase the file limit:
   ```bash
   ulimit -n unlimited
   ```
   Sometimes, this requires restarting the computer to take effect.

2. **Variant Calling using `bcftools`**:
   ```bash
   bcftools mpileup --threads 20 -f ../mephisto_alpha_renamed_polish.fasta_primary.fasta P3.3-new2.bam | bcftools call --threads 20 -mv -Ov -o P3.3-new2-calls.vcf
   ```

3. **Process VCF File**:
   - Use `parseVCF-freq.py` to process the VCF file, or use `parseVCF-freq3.py` if you want to only consider alternate reads greater than 3:
     ```bash
     python parseVCF-freq.py
     # or
     python parseVCF-freq3.py
     ```

4. **Filter SNPs**:
   - Run `filter-text-files.py` to ensure only SNPs are kept:
     ```bash
     python filter-text-files.py
     ```

5. **Compare Text Files**:
   - Use `compare-text-files.py` to compare filtered SNP files:
     ```bash
     python compare-text-files.py
     ```



