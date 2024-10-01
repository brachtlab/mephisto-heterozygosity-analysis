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
  ```

