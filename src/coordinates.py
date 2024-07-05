import os
import numpy as np
import pandas as pd
import pysam as ps
import random
from pyfaidx import Fasta
from scipy.stats import lognorm
from utils import parse_coordinates_arguments

class SimulateCoordinates:

    def __init__(self, args):
        self.type = args.type
        self.molecule = args.molecule
        self.number = args.number
        self.distribution = args.distribution
        self.length_min = args.length_min
        self.length_max = args.length_max
        self.mean = args.mean
        self.sd = args.sd
        self.genome_fasta = args.genome_fasta
        self.transcript_bed = args.transcript_bed
        self.output_bed = args.output_bed
        self.split = args.split

    def select_chromosomes(self):
        chromosomes = {}
        whole_genome_length = 0
        with Fasta(self.genome_fasta) as fasta:
            for name, record in fasta.items():
                if all(keyword not in name.lower() for keyword in ['alt', 'random', 'un', 'm', 'ki', 'gl']):  # Filter unusual chromosomes
                    chromosomes[name] = {'length': len(record), 'weight': 0}
        whole_genome_length = sum(chromosomes[contig]['length'] for contig in chromosomes)
        for contig, info in chromosomes.items():
            info['weight'] = info['length'] / whole_genome_length
        return chromosomes

    def find_next_splice_site(self, fasta, chromosome, start, end, target_sequence):
        circle_sequence = fasta.fetch(chromosome, start, end).upper()
        pos = circle_sequence.find(target_sequence)
        if pos == -1:
            return -1
        return start + pos

    def find_next_last_splice_site(self, fasta, chromosome, start, end, target_sequence):
        circle_sequence = fasta.fetch(chromosome, start, end).upper()
        pos = circle_sequence.rfind(target_sequence)
        if pos == -1:
            return -1
        return start + pos

    def simulate_coordinates(self, chromosomes):
        circle_bed = []
        n_circles = 0
        fasta = ps.FastaFile(self.genome_fasta)
        while n_circles < self.number:
            if self.distribution == 'uniform':
                circle_length = np.random.randint(self.length_min, self.length_max)
            elif self.distribution == 'lognormal':
                while True:
                    circle_length = int(lognorm.rvs(0.5, loc=0, scale=self.mean, size=1))
                    if self.length_min <= circle_length <= self.length_max:
                        break
            chromosome = np.random.choice(list(chromosomes.keys()), p=[chromosomes[contig]['weight'] for contig in chromosomes])
            if self.type == 'DNA':
                circle_start = np.random.randint(0, (chromosomes[chromosome]['length'] - circle_length))
            elif self.type == 'RNA':
                #transcript_bed = pd.read_csv(self.transcript_bed, sep='\t', header=None)
                # Example assuming the first column is integer (adjust dtype accordingly)
                transcript_bed = pd.read_csv(self.transcript_bed, sep='\t', header=None, dtype={1: str})
                filtered_transcripts = transcript_bed[(transcript_bed[1] == chromosome) & (transcript_bed[4] > circle_length)]
                if filtered_transcripts.empty:
                    continue
                transcript = filtered_transcripts.iloc[random.randint(0, len(filtered_transcripts)-1)]
                circle_start = np.random.randint(int(transcript[2]), int(transcript[3]) - circle_length)

            circle_end = circle_start + circle_length

            if self.split:
                # Start split sequence
                start_sequence = fasta.fetch(chromosome, circle_start-2, circle_start).upper()
                split_seq_start = self.split[:2]  # First half of the split sequence
                new_start = self.find_next_last_splice_site(fasta, chromosome, circle_start-100, circle_start-1, split_seq_start)
                if new_start == -1:
                    continue
                circle_start = new_start + 3

                # Check if circle_start is within the transcript
                if self.type == 'RNA':
                    if not (int(transcript[2]) <= circle_start < int(transcript[3])):
                        continue

                # End split sequence
                end_sequence = fasta.fetch(chromosome, circle_end, circle_end+2).upper()
                split_seq_end = self.split[2:]  # Second half of the split sequence
                if end_sequence != split_seq_end:
                    new_end = self.find_next_last_splice_site(fasta, chromosome, circle_start, circle_end, split_seq_end)
                    if new_end == -1:
                        continue
                    circle_end = new_end

                # Check if circle_end is within the transcript
                if self.type == 'RNA':
                    if not (int(transcript[2]) < circle_end <= int(transcript[3])):
                        continue

            circle_sequence = fasta.fetch(chromosome, circle_start, circle_end)

            circle_bed.append((chromosome, circle_start, circle_end))
            n_circles += 1

        return circle_bed


    def run(self):
        chromosomes = self.select_chromosomes()
        circle_bed = self.simulate_coordinates(chromosomes)
        try:
            with open(self.output_bed, 'w') as file:
                for line in circle_bed:
                    file.write('\t'.join(map(str, line)) + '\n')
        except (FileNotFoundError, PermissionError) as e:
            raise RuntimeError(f"Error writing to BED file '{self.output_bed}': {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulate coordinates of circular and linear molecules')
    parse_coordinates_arguments(parser)
    args = parser.parse_args()
    circles_instance = SimulateCoordinates(args)
    circles_instance.run()
