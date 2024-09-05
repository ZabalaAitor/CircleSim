import os
import argparse
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

    # Method to select chromosomes from the genome and calculate their weights based on length
    def select_chromosomes(self):
        chromosomes = {}
        whole_genome_length = 0
        with Fasta(self.genome_fasta) as fasta:
            # Filter out unwanted chromosome sequences (e.g., alt, random, unplaced sequences)
            for name, record in fasta.items():
                if all(keyword not in name.lower() for keyword in ['alt', 'random', 'un', 'm', 'ki', 'gl']):  # Filter unusual chromosomes
                    chromosomes[name] = {'length': len(record), 'weight': 0}
        # Calculate the total genome length
        whole_genome_length = sum(chromosomes[contig]['length'] for contig in chromosomes)
        # Assign a weight to each chromosome proportional to its length
        for contig, info in chromosomes.items():
            info['weight'] = info['length'] / whole_genome_length
        return chromosomes

    # Method to find the next occurrence of a target sequence within a specified range
    def find_next_splice_site(self, fasta, chromosome, start, end, target_sequence):
        circle_sequence = fasta.fetch(chromosome, start, end).upper()
        pos = circle_sequence.find(target_sequence)
        if pos == -1:
            return -1
        return start + pos

    # Method to find the previous occurrence of a target sequence within a specified range
    def find_next_last_splice_site(self, fasta, chromosome, start, end, target_sequence):
        circle_sequence = fasta.fetch(chromosome, start, end).upper()
        pos = circle_sequence.rfind(target_sequence)
        if pos == -1:
            return -1
        return start + pos

    # Main method to simulate the circle coordinates based on the genome and input arguments
    def simulate_coordinates(self, chromosomes):
        circle_bed = []
        n_circles = 0
        fasta = ps.FastaFile(self.genome_fasta)
        while n_circles < self.number:
            # Determine the length of the circle based on the distribution type
            if self.distribution == 'uniform':
                circle_length = np.random.randint(self.length_min, self.length_max)
            elif self.distribution == 'lognormal':
                while True:
                    circle_length = int(lognorm.rvs(s=self.sd, loc=0, scale=self.mean, size=1))
                    if self.length_min <= circle_length <= self.length_max:
                        break
            # Randomly choose a chromosome based on their weights
            chromosome = np.random.choice(list(chromosomes.keys()), p=[chromosomes[contig]['weight'] for contig in chromosomes])
            # Simulate DNA or RNA coordinates
            if self.type == 'DNA':
                circle_start = np.random.randint(0, (chromosomes[chromosome]['length'] - circle_length))
            elif self.type == 'RNA':
                # Load the transcript BED file and filter by the selected chromosome and circle length
                transcript_bed = pd.read_csv(self.transcript_bed, sep='\t', header=None, dtype={1: str})
                filtered_transcripts = transcript_bed[(transcript_bed[1] == chromosome) & (transcript_bed[4] > circle_length)]
                if filtered_transcripts.empty: # Skip if no transcripts match
                    continue
                # Select a random transcript and simulate the circle start within its range
                transcript = filtered_transcripts.iloc[random.randint(0, len(filtered_transcripts)-1)]
                circle_start = np.random.randint(int(transcript[2]), int(transcript[3]) - circle_length)
            # Define the circle end position
            circle_end = circle_start + circle_length
            # Handle the case where sequences are split based on a specific target sequence
            if self.split:
                # Check for a valid split sequence at the start
                start_sequence = fasta.fetch(chromosome, circle_start-2, circle_start).upper()
                split_seq_start = self.split[:2]  # First half of the split sequence
                new_start = self.find_next_last_splice_site(fasta, chromosome, circle_start-100, circle_start-1, split_seq_start)
                if new_start == -1: # If no valid start is found, skip this circle
                    continue
                circle_start = new_start + 3 # Adjust circle start

                # Ensure that the circle start is within the transcript for RNA molecules
                if self.type == 'RNA':
                    if not (int(transcript[2]) <= circle_start < int(transcript[3])):
                        continue

                # Ensure that the circle start is within the transcript for RNA molecules
                end_sequence = fasta.fetch(chromosome, circle_end, circle_end+2).upper()
                split_seq_end = self.split[2:]  # Second half of the split sequence
                if end_sequence != split_seq_end:
                    new_end = self.find_next_last_splice_site(fasta, chromosome, circle_start, circle_end, split_seq_end)
                    if new_end == -1:
                        continue
                    circle_end = new_end # Adjust circle end

                # Ensure that the circle end is within the transcript for RNA molecules
                if self.type == 'RNA':
                    if not (int(transcript[2]) < circle_end <= int(transcript[3])):
                        continue

            # Fetch the sequence for the circle
            circle_sequence = fasta.fetch(chromosome, circle_start, circle_end)

            # Append the simulated circle coordinates to the BED list
            circle_bed.append((chromosome, circle_start, circle_end))
            n_circles += 1

        return circle_bed

    # Main method to run the simulation and save the output to a BED file
    def run(self):
        chromosomes = self.select_chromosomes()  # Select valid chromosomes for simulation
        circle_bed = self.simulate_coordinates(chromosomes)  # Simulate the circle coordinates
        try:
            # Write the resulting coordinates to the output BED file
            with open(self.output_bed, 'w') as file:
                for line in circle_bed:
                    file.write('\t'.join(map(str, line)) + '\n')
        except (FileNotFoundError, PermissionError) as e:
            # Raise an error if there is an issue writing to the file
            raise RuntimeError(f"Error writing to BED file '{self.output_bed}': {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulate coordinates of circular and linear molecules')
    parse_coordinates_arguments(parser)
    args = parser.parse_args()
    circles_instance = SimulateCoordinates(args)
    circles_instance.run()
