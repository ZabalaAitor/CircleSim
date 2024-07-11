import os
import numpy as np
import pysam as ps
from Bio import SeqIO
from Bio.Seq import Seq
import random
import string
import math
from scipy.stats import beta
from pyfaidx import Fasta
from utils import parse_reads_arguments

def read_bed_file(input_bed):
    circle_bed = []
    with open(input_bed, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            fields[1] = int(fields[1])
            fields[2] = int(fields[2])
            circle_bed.append(fields)
    return circle_bed

kimura_matrix = {
    'A': {'A': 0.6, 'C': 0.2, 'G': 0.1, 'T': 0.1},
    'C': {'A': 0.2, 'C': 0.6, 'G': 0.1, 'T': 0.1},
    'G': {'A': 0.1, 'C': 0.1, 'G': 0.6, 'T': 0.2},
    'T': {'A': 0.1, 'C': 0.1, 'G': 0.2, 'T': 0.6}
}

def introduce_substitution(seq, position, substitution_matrix):
    if position >= 0 and position < len(seq):
        original_base = seq[position]
        if original_base in substitution_matrix:
            bases, probabilities = zip(*substitution_matrix[original_base].items())
            new_base = np.random.choice(bases, p=probabilities)
            seq[position] = new_base
    else:
        raise IndexError("Invalid position for substitution: position out of range.")

def introduce_insertion(seq, position):
    new_base = random.choice("ACGT")
    seq.insert(position, new_base)

def introduce_deletion(seq, position):
    if position >= 0 and position < len(seq):
        del seq[position]
    else:
        raise IndexError("Invalid position for deletion: position out of range.")

def introduce_sequencing_errors(seq, error_rate):
    error_bases = ['A', 'C', 'G', 'T']
    for i in range(len(seq)):
        if random.random() < error_rate:
            original_base = seq[i]
            new_base = random.choice([base for base in error_bases if base != original_base])
            seq[i] = new_base
    return seq

def introduce_mutations(seq, mutation_rate, substitution_matrix, error_rate):
    mutated_seq = list(seq)
    seq_length = len(seq)
    
    for i in range(seq_length):
        if random.random() < mutation_rate:
            mutation_type = random.choices(
                ['substitution', 'insertion', 'deletion'],
                [0.8, 0.1, 0.1],
                k=1
            )[0]
            if mutation_type == 'substitution':
                if i < len(mutated_seq):
                    introduce_substitution(mutated_seq, i, substitution_matrix)
            elif mutation_type == 'insertion':
                if i < len(mutated_seq):
                    introduce_insertion(mutated_seq, i)
            elif mutation_type == 'deletion':
                if i < len(mutated_seq):
                    introduce_deletion(mutated_seq, i)

    mutated_seq = introduce_sequencing_errors(mutated_seq, error_rate)
    return ''.join(mutated_seq)

class SimulateReads:
    def __init__(self, args, circle_bed):
        self.type = args.type
        self.molecule = args.molecule
        self.genome_fasta = args.genome_fasta
        self.sequence = args.sequence
        self.coverage = args.coverage
        self.reads_length = args.reads_length
        self.insert_length = args.insert_length
        self.alpha_value = args.alpha
        self.beta_value = args.beta
        self.mutation = args.mutation
        self.mutation_rate = args.mutation_rate
        self.error_rate = args.error_rate
        self.save_unmutated = args.save_unmutated
        self.input_bed = args.input_bed
        self.output_fastq = args.output_fastq
        self.circle_bed = circle_bed

    def simulate_reads(self):
        fasta = ps.FastaFile(self.genome_fasta)
        all_reads = []
        if self.save_unmutated:
            all_unmutated_reads = []
        alpha_value = self.alpha_value
        beta_value = self.beta_value
        if self.sequence == 'short':
            for circle_info in self.circle_bed:
                chromosome, circle_start, circle_end = circle_info
                circle_length = circle_end - circle_start
                reads = math.ceil((circle_length * self.coverage) / (self.reads_length * 2))
                for _ in range(round(reads)):
                    if self.molecule == 'linear':
                        insert_start = np.random.randint(circle_start, circle_end - self.insert_length)
                    else:
                        insert_start = int(beta.rvs(alpha_value, beta_value, loc=circle_start, scale=circle_length, size=1))
                    insert_end = insert_start + (self.insert_length % circle_length)
                    left_read_start = insert_start
                    left_read_end = left_read_start + self.reads_length
                    right_read_end = insert_end
                    right_read_start = insert_end - self.reads_length
                    code = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(10))
                    if right_read_end >= circle_end:
                        right_read_end -= circle_length
                        right_read_start -= circle_length
                    if (left_read_end <= circle_end) & (right_read_start >= circle_start):
                        left_read = fasta.fetch(chromosome, left_read_start, left_read_end)
                        right_read = fasta.fetch(chromosome, right_read_start, right_read_end)
                        if right_read_start < left_read_start:
                            label = 'DR'
                        else:
                            name = f'{code}|{chromosome}:{circle_start}-{circle_end}|{left_read_start}-{right_read_end}|CR '
                            label = 'CR'
                    elif (left_read_end > circle_end) & (right_read_start < circle_start):
                        left_read = fasta.fetch(chromosome, left_read_start, circle_end) + fasta.fetch(chromosome, circle_start, circle_start + left_read_end - circle_end)
                        right_read = fasta.fetch(chromosome, circle_end - (circle_start - right_read_start), circle_end) + fasta.fetch(chromosome, circle_start, right_read_end)
                        label = 'LRSR'
                    elif left_read_end > circle_end:
                        left_read = fasta.fetch(chromosome, left_read_start, circle_end) + fasta.fetch(chromosome, circle_start, circle_start + left_read_end - circle_end)
                        right_read = fasta.fetch(chromosome, right_read_start, right_read_end)
                        label = 'LSR'
                    elif right_read_start < circle_start:
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
                    right_read = str(Seq(right_read).reverse_complement())
                    left_read = left_read.upper()
                    right_read = right_read.upper()

                    if self.mutation:
                        mutated_left_read = introduce_mutations(left_read, self.mutation_rate, kimura_matrix, self.error_rate)
                        mutated_right_read = introduce_mutations(right_read, self.mutation_rate, kimura_matrix, self.error_rate)
                        all_reads.append((left_name, right_name, mutated_left_read, mutated_right_read))
                        if self.save_unmutated:
                            all_unmutated_reads.append((left_name, right_name, left_read, right_read))
                    else:
                        all_reads.append((left_name, right_name, left_read, right_read))

            self.write_fastq_files(self.output_fastq, all_reads, phred_score=40, sequence='short')
            if self.save_unmutated:
                output_path = self.output_fastq + 'unmutated_'
                self.write_fastq_files(output_path, all_unmutated_reads, phred_score=40, sequence='short')

        if self.sequence == 'long':
            for circle_info in self.circle_bed:
                chromosome, circle_start, circle_end = circle_info
                circle_length = circle_end - circle_start
                reads = math.ceil((circle_length * self.coverage) / (self.reads_length * 2))
                for _ in range(round(reads)):
                    read_start = int(beta.rvs(alpha_value, beta_value, loc=circle_start, scale=circle_length, size=1))
                    read_end = read_start + self.reads_length
                    while read_end > circle_end:
                        read_end -= self.reads_length
                    n_bsj = 0
                    n_bsj = math.ceil((self.reads_length / circle_length))
                    complete_bsj = n_bsj - 2
                    fasta = ps.FastaFile(self.genome_fasta)
                    start_sequence = fasta.fetch(chromosome, read_start, circle_end)
                    mid_sequence = fasta.fetch(chromosome, circle_start, circle_end)
                    end_sequence = fasta.fetch(chromosome, circle_start, read_end)
                    read = start_sequence + complete_bsj * mid_sequence + end_sequence
                    code = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(10))
                    name = f'{code}|{chromosome}:{circle_start}-{circle_end}|{read_start}-{read_end}'
                    read = read.upper()
                    all_reads.append((name, read))

            self.write_fastq_files(self.output_fastq, all_reads, phred_score=40, sequence='long')

    def generate_quality_line(self, phred_score, length):
        return chr(phred_score + 33) * length

    def write_fastq_files(self, file_path, all_reads, sequence, phred_score=40):
        if sequence == 'short':
            r1_file_path = f"{file_path}R1.fastq"
            r2_file_path = f"{file_path}R2.fastq"

            r1_records = []
            r2_records = []
            for left_name, right_name, left_read, right_read in all_reads:
                left_record = SeqIO.SeqRecord(Seq(left_read), id=f"{left_name} R1", description="")
                left_record.letter_annotations['phred_quality'] = [phred_score] * len(left_read)
                r1_records.append(left_record)

                right_record = SeqIO.SeqRecord(Seq(right_read), id=f"{right_name} R2", description="")
                right_record.letter_annotations['phred_quality'] = [phred_score] * len(right_read)
                r2_records.append(right_record)

            with open(r1_file_path, 'w') as r1_file:
                SeqIO.write(r1_records, r1_file, 'fastq')

            with open(r2_file_path, 'w') as r2_file:
                SeqIO.write(r2_records, r2_file, 'fastq')

        elif sequence == 'long':
            records = []
            for name, read in all_reads:
                record = SeqIO.SeqRecord(Seq(read), id=f"{name}", description="")
                record.letter_annotations['phred_quality'] = [phred_score] * len(read)
                records.append(record)
            with open(file_path, 'w') as file:
                SeqIO.write(records, file, 'fastq')

    def run(self):
        circle_bed = read_bed_file(self.input_bed)
        self.simulate_reads()

if __name__ == "__main__":
    args = parse_reads_arguments()
    circle_bed = read_bed_file(args.input_bed)
    reads_instance = SimulateReads(args, circle_bed)
    reads_instance.simulate_reads()
