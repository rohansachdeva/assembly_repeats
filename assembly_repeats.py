#!/usr/bin/env python3
#Rohan Sachdeva (rohansach@berkeley.edu)
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
from Bio.SeqUtils import GC
import os,sys,argparse,regex
from xopen import xopen
from tqdm import tqdm

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--in_fasta', action = 'store', dest = 'in_fasta', required = True, help = 'Fasta file containing sequences/scaffolds/contigs/assemblies')
parser.add_argument('-o', '--out_fasta', action = 'store', dest = 'out_fasta', required = True, help = 'Create a fasta with error sequences')
parser.add_argument('-t', '--table_out', action = 'store', dest = 'out_table', required = True, help = 'Output table destination/path')
parser.add_argument('--force', action = 'store_true', dest = 'force', required = False, help = 'Overwrite existing file(s)')
parser.add_argument('-c', '--circular_only', action = 'store_true', dest = 'circular_only', required = False)
parser.add_argument('-m', '--min_bp', action = 'store', required = False)

args = parser.parse_args()

fasta = os.path.expanduser(args.in_fasta)
out_fasta = os.path.expanduser(args.out_fasta)
out_table = os.path.expanduser(args.out_table)

if not args.force:
    if os.path.exists(out_fasta) and os.path.exists(out_table):
        print('Output fasta and table exist')
        sys.exit()
    elif os.path.exists(out_fasta):
        print('Output fasta exists')
        sys.exit()
    elif os.path.exists(out_table):
        print('Output table exists')
        sys.exit()
else:
    if os.path.exists(out_fasta):
        os.remove(out_fasta)
    if os.path.exists(out_table):
        os.remove(out_table)

overlap_length = 50

def write_fasta():
    fa = '>' + header + '\n' + seq
    out_fasta.write(fa + '\n')

def write_table():
    out_line = '\t'.join([id_, str(bp), str(gc), error_type, str(repeat_length)])
    out_table.write(out_line + '\n')

def check_direct_repeat(sequence):
    fraction_bp = int(0.2 * bp)
    fraction_seq = seq[:fraction_bp]
    rx = f'({fraction_seq})' + regex_subs
    fraction_count = regex.findall(rx, seq)
    fraction_count = len(fraction_count)
    if fraction_count >= 2:
        return len(fraction_seq)
    else:
        return False

def check_circular(sequence, overlap_length):
    seq_length = len(sequence)
    half_seq_length = int(seq_length / 2)
    beg_seq,end_seq = sequence[:50],sequence[-half_seq_length:]
    beg_seq_in_end_seq_index = end_seq.rfind(beg_seq)

    if beg_seq_in_end_seq_index != -1:
        end_match = end_seq[beg_seq_in_end_seq_index:]
        len_end_match = len(end_match)
        beg_match = seq[:len_end_match]

        if beg_match == end_match:
            return len_end_match
        else:
            return False

def check_palindrome(sequence):
    sequence_rc = reverse_complement(seq)
    if sequence == sequence_rc:
        return True

substitutions = '2'
regex_subs = '{s<=' + substitutions + '}'

with xopen(out_fasta, 'w') as out_fasta, xopen(out_table, 'w') as out_table:
    out_header = ['sequence_id', 'bp', 'gc', 'repeat_type', 'repeat_length']
    out_header = '\t'.join(out_header)
    out_table.write(out_header + '\n')

    for header,seq in tqdm(SimpleFastaParser(xopen(fasta)), unit=' Sequences processed'):
        header,seq = header.strip(),seq.strip()
        seq = seq.upper()
        id_ = header.split(' ')[0]
        bp,gc = len(seq),round(GC(seq),2)

        if args.circular_only:
            circular_repeat_length = check_circular(seq,overlap_length)

            if circular_repeat_length:
                error_type = 'potential_circular'
                repeat_length = circular_repeat_length

            else:
                continue

        else:

            try:
                seq_rc = reverse_complement(seq)
            except ValueError:
                continue

            if check_palindrome(seq):
                error_type = 'full_palindrome'
                repeat_length = bp

            else:
                direct_repeat_length = check_direct_repeat(seq)

                if direct_repeat_length:
                    error_type = 'direct_repeat'
                    repeat_length = direct_repeat_length

                else:
                    circular_repeat_length = check_circular(seq,overlap_length)

                    if circular_repeat_length:
                        error_type = 'potential_circular'
                        repeat_length = circular_repeat_length

                    else:
                        continue

        write_fasta()
        write_table()
