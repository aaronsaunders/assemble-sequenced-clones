#!/usr/bin/python

import phrap
import os
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import shutil
from Bio import SeqIO
import sys
"""
This script is a wrapper for the functions in the phrap.py module.

"""

args = sys.argv[1:]
infilepath = args[0]
rootpath = args[1]

##########################################
# PARAMETERS

#! this should be in a parameters file for records
fpcrprimer = Seq('AGAGTTTGATCCTGGCTCA', IUPAC.IUPACAmbiguousDNA())
rpcrprimer = Seq('GGYTACCTTGTTACGACTT', IUPAC.IUPACAmbiguousDNA())

# Expects filename base delim sample delim primer '.' ext
# eg. seq-A1-m13f.fasta delim = '-'
filename_delim = '-'

##########################################################
"""
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

"""
# print parameters
print 'INPUT PARAMETERS'
print 'pwd : {0}'.format(os.getcwd())
print 'infile path : {0}'.format(infilepath)
print 'root path : {0}'.format(rootpath)
print 'fprimer: ' + fpcrprimer
print 'rprimer: ' + rpcrprimer + '\n'

# set out the folders into which the data will be stored
reads_path = os.path.join(rootpath, 'reads')
samples_path = os.path.join(rootpath, 'samples')
contigs_path = os.path.join(rootpath, 'contigs')
paths = [ rootpath, reads_path, samples_path, contigs_path ]
for path in paths:
    if not os.path.exists(path):
        os.mkdir(path)
        print 'made dir: ' + os.path.basename(path)

# convert the phd files to fasta + qual
if not os.path.exists(infilepath):
    print 'infilepath invalid!'
    sys.exit()
    
# make a list of the files(reads) in the infilepath
# (base, sample): [ read1, read2, ...readn ]

os.chdir(infilepath)
print os.getcwd()

reads = [ phrap.parse_readfilename(filename, filename_delim)  
          for filename in os.listdir(infilepath) 
          if filename.partition('.')[2].startswith('phd') ]

samples = {}
for (sample, primer, filename) in reads:
    if sample in samples:
        samples[sample].append(filename)
    if sample not in samples:
        samples[sample] = [filename]

print 'sample reads grouped as:'
for sample, reads in samples.items():
    print '{0} : {1}'.format( sample, ', '.join( reads )  )

toassemble = []
print 'writing to {0}'.format( samples_path )
for sample, reads in samples.items():
    sequences = [ SeqIO.read(open(filename, 'r'), 'phd') 
                  for filename in reads ]
    print 'writing file: {0} .fasta .qual'.format(sample)
    outfilepath = os.path.join(
        samples_path, '{0}.fasta'.format(sample))
    with open(outfilepath, 'w') as outfh:
        fcount= SeqIO.write(sequences, outfh, 'fasta')
    if fcount > 1:
        toassemble.append('{0}.fasta'.format(sample))
    with open(outfilepath + '.qual', 'w') as outfh:
        SeqIO.write(sequences, outfh, 'qual')

print ''
os.chdir(samples_path)

# run phrap
print 'assembling samples'
with open('phrap_logfile.log', 'w') as logfh:
    assembled = [ phrap.assemble(filename, logfh) 
                  for filename in toassemble ]

contigs = [ name + '.contigs' for name in toassemble  ]

# trim contig at pcr primersfor (base, sample), reads in samples.items():
for filename in contigs:
    record = phrap.trim_pcrprimers(filename, fpcrprimer, rpcrprimer)
    with open(filename+'.trim.fasta', 'w') as fh:
        SeqIO.write(record, fh, 'fasta')    
    outfilepath = os.path.join(contigs_path, filename + '.trim.fasta')
    with open(outfilepath, 'w') as fh:
        SeqIO.write(record, fh, 'fasta')    

# TODO: delete the qual file

print 'Contigs assembled:' 
for filename in os.listdir(contigs_path):
    print filename
