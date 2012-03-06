#!/usr/bin/python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import os.path
import re
import subprocess

def phd2fasta(filename, inpath, qualflag = 'TRUE'):
    '''Biopython implementation of phd2fasta '''
    if filename[-6:] == '.phd.1':
        stubname = filename[:-6]
    elif filename [-4:] == '.phd':
        stubname = filename[:-4]
    else:
        return None
    filename = os.path.join(inpath, filename)
    print filename
    fastaname = stubname + '.fasta'
    result = SeqIO.convert(filename, "phd", fastaname, "fasta")
    if result == 1:
	print 'converted' + fastaname
    if qualflag:
        qualname = stubname + '.fasta.qual'
        print qualname
        SeqIO.convert(filename, "phd", qualname, "qual")
    
    return


def parse_readfilename(filename, filename_delimiter):
    """ """
    filename, delim, ext = filename.partition('.')
    base, sample, primer = filename.split(filename_delimiter, 2)
    
    return (base + '-' + sample,  primer, filename + '.' + ext)

def assemble(filename, logfh):
    print 'asembling {0}'.format( filename )
    cmd = 'phrap -ace {fname}'.format(fname=filename)
    process = subprocess.Popen(cmd, shell = True, \
                           stdout=subprocess.PIPE, \
                           stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    process.wait()
    logfh.write(stdout)
    logfh.write('\n-----\n')

    return filename

def trim_pcrprimers(filename, fprimer, rprimer):
    """ """
    records = SeqIO.QualityIO.PairedFastaQualIterator(
        open(filename, 'rU'),
        open(filename+'.qual', 'rU'))
    outrecs = []
    print filename
    
    for record in records:
        result = trim_fprimer(record, fprimer) 
        if not result[1]:
            record, match = trim_fprimer_minus_strand(
                                        record, fprimer)
            if match:
                record = reverse_complement_SeqRecord(record)
                # TODO: error message if not match
        record, match = trim_rprimer(record, rprimer)
        outrecs.append(record)

    # TODO: add qual score screen
    
    return outrecs

def regex_from_string(seq):
    '''Makes an IUPAC DNA sequence string into a regular expression'''
    bases = { 'A' : 'A', 'T': 'T', 'G': 'G', 'C': 'C',
	      'R' : '[AG]', 'Y': '[CT]', 'S': '[GC]', 
	      'W' : '[AT]', 'K': '[GT]', 'M': '[AC]',
	      'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]',
	      'V': '[ACG]', 'N': '[AGCT]'}
    seq = seq .upper()
  
    return ''.join([ bases[character] for character in seq ])

def match_primer(record, primer_string):
    """Takes a primer and returns first match in the first 150 bp """
    regex = regex_from_string(primer_string)
    m = re.search(regex, record.seq.tostring())
    if m:
        return m.span()

    return 

def trim_fprimer(record, primer):
    """takes a list of SeqObjects or sequences and trims 5' at primer """
    primer = primer.tostring()
    m = match_primer(record, primer)
    if not m:
        return record, False
    match_start, match_end = m
    record = record[match_end:]

    return record, True

def trim_fprimer_minus_strand(record, primer):
    primer = primer.reverse_complement().tostring()
    m = match_primer(record, primer)
    if not m:
        return record, False
    match_start, match_end = m
    record = record[:match_start]

    return record, True

def trim_rprimer(record, primer):
    """
    takes a Seq or seqobj and primer (opt. qual file)
    trims the reverse primer and replaces with G (sense) and C (anti) 
    """
    primer = primer.reverse_complement().tostring()
    m = match_primer(record, primer)
    if not m:
        return record, False
    match_start, match_end = m
    record = record[:match_end]

    return record, True

def reverse_complement_SeqRecord(record):
    """Returns  a  new  SeqRecord  with  the  reverse  complement  sequence."""
    return  SeqRecord(seq = record.seq.reverse_complement(), \
                    id = record.id, description = "reverse complement" )

def trim_read_quality(record):
    #if min(record.letter_annotations["phred_quality"]) >= 20:
    
    return 
