


---
# 1. Measuring Allele Frequency Between Parent and Progeny 

### 1.1 **Map with `minimap2 onto HAPLOID genome assembly of haplotype 1`**:
  - Haplotype 1 can be arbitrarily determined, and snps will represent the other haplotype.
   - Use the following `LSF` script for mapping reads with `minimap2`:
     ```bash
     #!/bin/bash
     #BSUB -J minimap2
     #BSUB -q long
     #BSUB -o minimap5_log.txt
     #BSUB -e minimap5_Error.txt
     #BSUB -n 24

     minimap2 -t 24 -a meph-pri.mmi P3.3.fastq.gz >P3.3.sam
     ```
   - This maps `P3.3-omega.fastq.gz` onto the `meph-pri.mmi` reference and outputs the SAM file `P3.3.sam`.
   - It is critical to map onto a high-quality haploid assembly. The pipeline assumes the snps relative to this assembly are the other genotype. 

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

### 1.4 **Run bcftools to call variants**:
   - Use the following `LSF` script to run `bcftools`:
```bash
     #!/bin/bash
     #BSUB -J bcftools_JB7
     #BSUB -q normal
     #BSUB -o bcftools_Log-7.txt
     #BSUB -e bcftools_Error-7.txt
     #BSUB -n 48

     bcftools mpileup -f mephisto_alpha_renamed_polish.fasta_primary.fasta P3.3_sorted.bam | bcftools call -mv -Ov -o P3.3.vcf
 ```
   - This script generates the variant calls in the VCF file `P3.3.vcf`.

### 1.5 **Convert vcf file to tab-delimited snp-only file for plotting**:
- **step 1.5.1** : Use parseVCF-freq2.py to process the VCF file into a tab-delimited text file. Input is FILENAME.vcf, Output is FILENAME.vcf_columns2.txt
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
     
# 2. Calling Recombinant Reads

This pipeline performs recombination detection in reads from a sam file. Below are the steps and scripts used in the process:

### 2.1 Map Reads to Haploid Genome of Haplotype 1:
- Follow step 1.1, 1.2, 1.3, and 1.4 above to map the heterozygote (parental) reads (in our case, the reads are from parent we labeled P3, raw reads are the file P3.fastq.gz) reads onto the (haplotype 1) genome, generating a .vcf file for identifying the haplotype 2 snps. This file will be a reference for looking for identification of recombinant reads.
 

### 2.2 Identify the snps in haplotype 2:
- Call P3 variants by running `parseVCF-highConfidenceSNPs.py`, ensuring at least 5 reads support both reference and alternative allele calls. Input is FILENAME.vcf and output is FILENAME.vcf_high-conf-snps.txt.
 ```bash
     ./parseVCF-highConfidenceSNPs.py P3.vcf
 ```
### 2.3 Filter SNPs:
- Use `filter-text-files.py` to filter the SNPs. Input is FILENAME.vcf_high-conf-snps.txt and output is FILENAME.vcf_high-conf-snps.txt_snps_only.txt
 ```bash
     ./filter-text-files.py P3.vcf_high-conf-snps.txt
 ```

### 2.4 Add Genomic Context to Haplotype 2:
- Run `add-context_fixed.py` to add genomic context to the SNP file. (Context is the previous 7 bp and 12 bp relative to the snp. Both are used for identifying snps within the reads using find-recombination.py.) Input is FILENAME.vcf_high-conf-snps.txt_snps_only.txt and the genome.fasta, and output is FILENAME.vcf_high-conf-snps.txt_snps_only.txt_contextFIXED.txt
```bash
     ./add-context_fixed.py P3.vcf_high-conf-snps.txt_snps_only.txt genome.fasta
 ```
### 2.5 Detect Forward Recombination Events:
- 2.5.1 Generate the SAM File.
  - **important** The SAM file from previous must not be used. You want a sorted read set to make the script run faster. Therefore, delete the SAM file from previous and use your sorted, indexed bam file to re-generate a SAM file:
```bash
     rm P3.1.sam
     samtools view -@ 25 -o P3.1.sam P3.1_sorted.bam 
 ```
- 
- 2.5.2 Run `find-recombination.py` which requires the SAM file (containing the reads and their mapping positions) along with the snp file (FILENAME.vcf_high-conf-snps.txt_snps_only.txt_contextFIXED.txt). The output is SAMFILE.sam_reads_calls2.txt.    
```bash
     ./find-recombination.py P3.1.sam P3.vcf_high-conf-snps.txt_snps_only.txtcontextFIXED.txt 
 ```
The output file has 10 columns and one row per read (if the read is long enough to be scored). The recombination-flag is the key because it is set if the read is thought to be recombinant. Error-counts is the number of snps that did not match either a reference or alt basepair, yet the snp site was present in the read (based on the context search). These are presumed errors owing to the fairly high error rate of the long reads.The 'call-list' shows a list of each snp call within the read, whether reference, 'ref' or alternate, 'alt'.  

  

### 2.6 Reverse Complement Mapping:
-The previous steps only check for recombinant reads in the forward direction. To complete the analysis, map onto the reverse complemented genome (of haplotype 1), then call snps and then combine all calls together (steps 2.6-2.10).
-First, reverse complement the genome with python script 'reverse_complement.py'. Input is genome.fasta, output is genome.fasta_reverse_complement.fasta:
```bash
     ./reverse_complement.py genome.fasta
 ```
### 2.7 Map Parental Reads to Reverse Complement Genome to identify reverse complement snps in Haplotype 2
- Map parental read file (in this case, P3.fastq.gz,) onto the reverse complement of the genome, convert to vcf, then to text file, and add context following steps 2.1, 2.2, 2.3, and 2.4. **Important** when mapping use a prefix like 'rc_P3.sam' as the output to track that this is the reverse complement. You will end up with a file like `rc_P3.vcf_high-conf-snps.txt_snps_only.txt_contextFIXED.txt` which will be used as input in the next step.

### 2.7 Detect Reverse Complement Recombination Events:
- Run `rc-find-recombination.py` using the reverse complement SNPs file. **Important** the only difference between find-recombination.py and rc-find-recombination.py is that rc-find-recombination.py saves the position of the 3' end of the read in the calls2.txt file so that when the overal coordinates are reversed (in step 2.8) they will be correctly combined with the forward mapping file.
-  This generates:  
  `rc_reads_calls2.txt`
  ```bash
     ./rc-find-recombination.py P3.1.sam rc_P3.vcf_high-conf-snps.txt_snps_only.txt_contextFIXED.txt
 ```
- generating rc_reads_calls2.txt file.

### 2.8 Reverse Complement Coordinates:
- Run `reverse-complement-coordinates-of-call-file.py` (requires reverse complement genome input) to convert the coordinates back to the original genome's direction. This generates:  
  `rc_reads_calls2.txt_rc.txt` which now can be combined with the forward read calls (it has the positions of the reads correctly marked).
```bash
     ./recombination-omega/reverse-complement-coordinates-of-call-file.py rcP3.1.sam_reads_calls2.txt reverse_complement_genome.fasta
 ```

### 2.9 Merge Forward and Reverse Call Files:
- Combine forward and reverse complement calls using `cat`:
```bash
  cat P3.1sam_reads_calls2.txt rcP3.sam_reads_calls.txt_rc.txt > P3.1-combined-calls.txt
  ```

### 2.10 Sorting:
- Open `combined-calls.txt` in Excel.
- Sort by `contig` and `position`.
- Save as a tab-delimited text file. You can overwrite the old file name so it doesn't create yet another file. 

### 2.11 Add read lengths to the file:
- Run `add-read-lengths.py`: Input is the sorted combined-calls.txt, output is combined-calls.txt_added-lengths.txt
  ```bash
  ./add-read-lengths.py P3.1-combined-calls.txt
  ```

### 2.12 Compare Cluster Analyses:
- Run `find-clusters-compare.py` to compare two cluster analyses. Each filename is followed by the number of reads that constitute a cluster. The example below is for a 10-read threshold. 
  ```bash
  find-clusters-compare.py P3-combined-calls.txt_added-lengths 10 P3.1-combined-calls.txt_added-lengths 10
  ```
  





