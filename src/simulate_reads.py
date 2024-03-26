import os
import numpy as np
import pysam as ps
from Bio import SeqIO
from Bio.Seq import Seq
import random
import string
import math
from scipy.stats import beta
from utils import parse_reads_arguments 

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
    def __init__(self, args, circle_bed):
        self.genome_fasta = args.genome_fasta
        self.sequence = args.sequence
        self.coverage = args.coverage
        self.reads_length = args.reads_length
        self.insert_length = args.insert_length
        self.input_bed = args.input_bed
        self.output_fastq = args.output_fastq
        self.circle_bed = circle_bed

    def simulate_reads(self):
        fasta = ps.FastaFile(self.genome_fasta)
        all_reads = []
        alpha_value = 0.5
        beta_value = 0.5
        if self.sequence == 'short':
            for circle_info in self.circle_bed:
                chromosome, circle_start, circle_end = circle_info
                circle_length = circle_end - circle_start
                reads = math.ceil((circle_length * self.coverage) / (self.reads_length * 2))
                for _ in range(round(reads)): 
                    insert_start = int(beta.rvs(alpha_value, beta_value, loc=circle_start, scale=circle_length, size=1))
                    insert_end = insert_start + (self.insert_length % circle_length)
                    n_bsj = 0
                    n_bsj = math.ceil(self.insert_length / circle_length)  #### insert_length // circle_length
                    left_read_start = insert_start
                    left_read_end = left_read_start + self.reads_length
                    right_read_end = insert_end
                    right_read_start = insert_end - self.reads_length
                    code = random_string = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(10))
                    # CORRECCIÓN: No teníamos en cuenta si el R2 pasaba del círculo, solo si pasaba pero era mayor que la longitud del círculo en sí
                    if right_read_end >= circle_end:
                        right_read_end -= circle_length
                        right_read_start -= circle_length
                    if (left_read_end <= circle_end) & (right_read_start >= circle_start):
                        left_read = fasta.fetch(chromosome, left_read_start, left_read_end)
                        right_read = fasta.fetch(chromosome, right_read_start, right_read_end)
                        if right_read_start < left_read_start:  # DISCORDANT
                            label = 'DR'
                        else:  # CONCORDANT
                            name = f'{code}|{chromosome}:{circle_start}-{circle_end}|{left_read_start}-{right_read_end}|CR '
                            label = 'CR'
                    elif (left_read_end > circle_end) & (right_read_start < circle_start):  # LEFT RIGHT SPLIT READ
                        left_read = fasta.fetch(chromosome, left_read_start, circle_end) + fasta.fetch(chromosome, circle_start, circle_start + left_read_end - circle_end)
                        right_read = fasta.fetch(chromosome, circle_end - (circle_start - right_read_start), circle_end) + fasta.fetch(chromosome, circle_start, right_read_end)
                        label = 'LRSR'
                    elif left_read_end > circle_end:  # LEFT SPLIT READ
                        left_read = fasta.fetch(chromosome, left_read_start, circle_end) + fasta.fetch(chromosome, circle_start, circle_start + left_read_end - circle_end)
                        right_read = fasta.fetch(chromosome, right_read_start, right_read_end)
                        label = 'LSR'
                    elif right_read_start < circle_start:  # RIGHT SPLIT READ
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
                    left_read = left_read.upper()  #########
                    right_read = right_read.upper()  ########
                    all_reads.append((left_name, right_name, left_read, right_read))

            # Write reads to a fastq file
            self.write_fastq_files(self.output_fastq, all_reads, phred_score=40, sequence='short')

        if self.sequence == 'long':
            for circle_info in self.circle_bed:
                chromosome, circle_start, circle_end = circle_info
                circle_length = circle_end - circle_start
                reads = math.ceil((circle_length * self.coverage) / (self.reads_length * 2))  ##### Round up!!!
                for _ in range(round(reads)): 
                    read_start = np.random.randint(circle_start, circle_end)  # CHANGE!!!
                    read_end = read_start + self.reads_length
                    while read_end > circle_end:
                        read_end -= self.reads_length
                    n_bsj = 0
                    n_bsj = math.ceil((self.reads_length / circle_length))  #### CHECK
                    complete_bsj = n_bsj - 2
                    fasta = ps.FastaFile(self.genome_fasta)
                    start_sequence = fasta.fetch(chromosome, read_start, circle_end)
                    mid_sequence = fasta.fetch(chromosome, circle_start, circle_end)
                    end_sequence = fasta.fetch(chromosome, circle_start, read_end)
                    read = start_sequence + complete_bsj * mid_sequence + end_sequence 
                    code = random_string = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(10))
                    name = f'{code}|{chromosome}:{circle_start}-{circle_end}|{read_start}-{read_end}'
                    read = read.upper()  #########
                    all_reads.append((name, read))

            # Write reads to a fastq file
            self.write_fastq_files(self.output_fastq, all_reads, phred_score=40, sequence='long')

    def generate_quality_line(phred_score, length):
        return chr(phred_score + 33) * length

    def write_fastq_files(self, file_path, all_reads, sequence, phred_score=40):
        if sequence == 'short':
            r1_file_path = f"{file_path}R1.fastq"
            r2_file_path = f"{file_path}R2.fastq"

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

    def run(self):
        circle_bed = read_bed_file(self.input_bed)
        self.simulate_reads()

if __name__ == "__main__":
    args = parse_reads_arguments()  
    circle_bed = read_bed_file(args.input_bed)
    reads_instance = SimulateReads(circle_bed, args)
    reads_instance.simulate_reads()
