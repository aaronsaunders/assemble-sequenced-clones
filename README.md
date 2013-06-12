Assemble sanger sequenced clones with Phrap
===========================================

## Requirements

 - python 2.7
 - Biopython python module
 - phrap installed and in PATH

Reads should be in a file format convertable by Biopython SeqIO (Eg. php, ab1).

## Process

 1. make a reads, samples, contigs folder
 2. list the phd reads in the input folder
 3. group them into samples in the list
 4. open the reads of each sample
    - export to samples folder as one merged file per sample in both fasta and qual format
 5. add the '.f' or '.r' to the name of the read as appropriate
 5. run phrap on the samples folder
 6. copy the contig+contig.qual into the contigs folder
 7. open each:
     - search for the pcr primers
     - slice the sequences at the primers
     - slice for bad quality
 
