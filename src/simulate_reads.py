import os
import argparse
from CircSim import parse_arguments
import numpy as np
import pysam as ps
from Bio import SeqIO
from Bio.Seq import Seq
import random
import string
import math

def read_bed_file(input_bed):
    circle_bed = []

    with open(input_bed, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            # Convert start and end positions to integers
            fields[1] = int(fields[1])
            fields[2] = int(fields[2])
            circle_bed.append(fields)

    return circle_bed

class SimulateReads:
    def __init__(self, coverage, reads_length, insert_length, circle_bed, genome_fasta, output_fastq):
        self.coverage = coverage
        self.reads_length = reads_length
        self.insert_length = insert_length
        self.circle_bed = circle_bed
        self.genome_fasta = genome_fasta
        self.output_fastq = output_fastq

    def simulate_reads(coverage, reads_length, insert_length, circle_bed, genome_fasta):
    fasta = ps.FastaFile(genome_fasta)
        all_reads = []
        
        if sequence == 'short':
            for circle_info in circle_bed:
                chromosome, circle_start, circle_end = circle_info
                circle_length = circle_end - circle_start
                reads = math.ceil((circle_length * coverage) / (reads_length * 2)) ####### Round up!!!!!!!
                for _ in range(round(reads)): 
                    insert_start = np.random.randint(circle_start, circle_end) # CHANGE!!!
                    insert_end = insert_start + (insert_length % circle_length)
                    #n_bsj = math.ceil(insert_length / circle_length) #### insert_length // circle_length
                    left_read_start = insert_start
                    left_read_end = left_read_start + reads_length
                    right_read_end = insert_end
                    right_read_start = insert_end - reads_length
                    code = random_string = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(10))

                    # CORRECCIÓN: No teniamos en cuenta si el R2 pasaba del círculo, solo si pasaba pero era mayor que la longitud del circulo en sí
                    if right_read_end >= circle_end:
                        right_read_end -= circle_length
                        right_read_start -= circle_length             
                    if (left_read_end <= circle_end) & (right_read_start >= circle_start):
                        left_read = fasta.fetch(chromosome, left_read_start, left_read_end)
                        right_read = fasta.fetch(chromosome, right_read_start, right_read_end)
                        if right_read_start < left_read_start: # DISCORDANT
                            n_bsj += 1
                            label = 'DR'
                        else: # CONCORDANT
                            name = f'{code}|{chromosome}:{circle_start}-{circle_end}|{left_read_start}-{right_read_end}|CR '
                            label = 'CR'
                    elif (left_read_end > circle_end) & (right_read_start < circle_start): # LEFT RIGHT SPLIT READ
                        left_read = fasta.fetch(chromosome, left_read_start, circle_end) + fasta.fetch(chromosome, circle_start, circle_start + left_read_end - circle_end)
                        right_read = fasta.fetch(chromosome, circle_end - (circle_start - right_read_start), circle_end) + fasta.fetch(chromosome, circle_start, right_read_end)
                        label = 'LRSR'
                    elif left_read_end > circle_end: # LEFT SPLIT READ
                        left_read = fasta.fetch(chromosome, left_read_start, circle_end) + fasta.fetch(chromosome, circle_start, circle_start + left_read_end - circle_end)
                        right_read = fasta.fetch(chromosome, right_read_start, right_read_end)
                        label = 'LSR'
                    elif right_read_start < circle_start: # RIGHT SPLIT READ
                        left_read = fasta.fetch(chromosome, left_read_start, left_read_end)
                        right_read = fasta.fetch(chromosome, circle_end - (circle_start - right_read_start), circle_end) + fasta.fetch(chromosome, circle_start, right_read_end)
                        label = 'RSR'
                    else:
                        raise Exception(f"""CIR_START: {circle_start} | CIR_END: {circle_end}
                                            INS_START: {insert_start} | INS_END: {insert_end}
                                            CIR_LEN: {circle_end - circle_start} | INS_LEN: {insert_length}
                                            R1: {left_read_start}-{left_read_end}
                                            R2: {right_read_start}-{right_read_end}""")
                    name = f'{code}|{chromosome}:{circle_start}-{circle_end}|{left_read_start}-{right_read_end}|{label} '
                    left_name = name + '1:N:0:CGCTGTG'
                    right_name = name + '2:N:0:CGCTGTG'
                    # Create the reverse-complement read
                    right_read = str(Seq(right_read).reverse_complement())
                    left_read = left_read.upper() #########
                    right_read = right_read.upper() ########
                    all_reads.append((left_name, right_name, left_read, right_read))
                
            # Write reads to a fastq file
            write_fastq_files(output_fastq, all_reads, phred_score=40, sequence='short')

        if sequence == 'long':
            for circle_info in circle_bed:
                chromosome, circle_start, circle_end = circle_info
                circle_length = circle_end - circle_start
                reads = math.ceil((circle_length * coverage) / (reads_length * 2)) ##### Round up!!!
                for _ in range(round(reads)): 
                    read_start = np.random.randint(circle_start, circle_end) # CHANGE!!!
                    read_end = read_start + reads_length
                    while read_end > circle_end:
                        read_end -= reads_length
                    n_bsj = math.ceil((reads_length / circle_length)) #### CHECK
                    complete_bsj = n_bsj - 2
                    fasta = ps.FastaFile(genome_fasta)
                    start_sequence = fasta.fetch(chromosome, read_start, circle_end)
                    mid_sequence = fasta.fetch(chromosome, circle_start, circle_end)
                    end_sequence = fasta.fetch(chromosome, circle_start, read_end)
                    read = start_sequence + complete_bsj * mid_sequence + end_sequence 
                    code = random_string = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(10))
                    name = f'{code}|{chromosome}:{circle_start}-{circle_end}|{read_start}-{read_end}'
                    read = read.upper() #########
                    all_reads.append((name, read))
            
            # Write reads to a fastq file
            write_fastq_files(output_fastq, all_reads, phred_score=40, sequence='long')

    def generate_quality_line(phred_score, length):
        return chr(phred_score + 33) * length

    def write_fastq_files(file_path, all_reads, phred_score=40, sequence='short'):
        if sequence == 'short':
            r1_file_path = f"{file_path}_R1.fastq"
            r2_file_path = f"{file_path}_R2.fastq"

            r1_records = []
            r2_records = []
            for left_name, right_name, left_read, right_read in all_reads:
                # Create SeqRecord for left read
                left_record = SeqIO.SeqRecord(Seq(left_read), id=f"{left_name} R1", description="")
                left_record.letter_annotations['phred_quality'] = [phred_score] * len(left_read)
                r1_records.append(left_record)

                # Create SeqRecord for right read
                right_record = SeqIO.SeqRecord(Seq(right_read), id=f"{right_name} R2", description="")
                right_record.letter_annotations['phred_quality'] = [phred_score] * len(right_read)
                r2_records.append(right_record)
            
            # Write R1 records to the R1 file
            with open(r1_file_path, 'w') as r1_file:
                SeqIO.write(r1_records, r1_file, 'fastq')

            # Write R2 records to the R2 file
            with open(r2_file_path, 'w') as r2_file:
                SeqIO.write(r2_records, r2_file, 'fastq')
                    
        elif sequence == 'long':
            records = []
            for name, read in all_reads:
                # Create SeqRecord 
                record = SeqIO.SeqRecord(Seq(read), id=f"{name}", description="")
                record.letter_annotations['phred_quality'] = [phred_score] * len(read)
                records.append(record)
            # Write records to the file
            with open(file_path, 'w') as file:
                SeqIO.write(records, file, 'fastq')

if __name__ == "__main__":
    args = parse_arguments()
    circle_bed = read_bed_file(args.output_bed)
    reads_instance = SimulateReads(args.coverage, args.reads_length, args.insert_length, circle_bed, args.genome_fasta, args.output_fastq)
    reads_instance.simulate_reads()