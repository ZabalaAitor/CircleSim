import os
from utils import parse_coordinates_arguments
import numpy as np
import pandas as pd
import pysam as ps
import random
from pyfaidx import Fasta
from scipy.stats import lognorm

class SimulateCoordinates:

    def __init__(self, args):
        self.type = args.type
        self.number = args.number
        self.distribution = args.distribution
        self.length_min = args.length_min
        self.length_max = args.length_max
        self.mean = args.mean   
        self.sd = args.sd
        self.genome_fasta = args.genome_fasta
        self.transcript_bed = args.transcript_bed
        self.output_bed = args.output_bed

    def select_chromosomes(self):
        chromosomes = {}
        whole_genome_length = 0
        # Filter unusual chromosomes
        with Fasta(self.genome_fasta) as fasta:
            for name, record in fasta.items():
                if all(keyword not in name.lower() for keyword in ['alt', 'random', 'un', 'm', 'ki', 'gl']): # Specify unusual chromosomes
                    chromosomes[name] = {'length': len(record), 'weight': 0}              
        # Calculate the total genome length
        whole_genome_length = sum(chromosomes[contig]['length'] for contig in chromosomes)       
        # Calculate and add weights to the chromosomes dictionary for probability
        for contig, info in chromosomes.items():
            info['weight'] = info['length'] / whole_genome_length
        return chromosomes

    def simulate_coordinates(self, chromosomes):
        circle_bed = []
        n_circles = 0
        fasta = ps.FastaFile(self.genome_fasta)     
        # Create circles 
        while n_circles < self.number:
            # Select distribution
            if self.distribution == 'uniform':  
                # Select a circle length between a max and min value by a uniformal distribution
                circle_length = np.random.randint(self.length_min, self.length_max) 
            elif self.distribution == 'lognormal':
                while True:
                    # Generate circle length from lognormal distribution
                    circle_length = int(lognorm.rvs(0.5, loc=0, scale=self.mean, size=1))
                    # Check if circle length falls within specified range
                    if self.length_min <= circle_length <= self.length_max:
                        break  # Exit loop if valid length is generated
            # Select chromosome base on weight
            chromosome = np.random.choice(list(chromosomes.keys()), p=[chromosomes[contig]['weight'] for contig in chromosomes])
            # Specify molecule type
            if self.type == 'DNA': # All genome
                # Randomly select a position within the selected transcript
                circle_start = np.random.randint(0, (chromosomes[chromosome]['length'] - circle_length))
            elif self.type == 'RNA': # Transcript
                transcript_bed = pd.read_csv(self.transcript_bed, sep='\t', header=None)
                # Filter transcripts that are not in the selected chromosome and have a length greater than circle length
                filtered_transcripts = transcript_bed[(transcript_bed[1] == chromosome) & (transcript_bed[4] > circle_length)]
                if filtered_transcripts.empty:
                    continue  # Skip if no valid transcripts found for this chromosome
                # Select a random transcript from filtered 
                transcript = random.choice(filtered_transcripts.values)
                # Randomly select a position within the selected transcript
                circle_start = np.random.randint(int(transcript[2]), int(transcript[3]) - circle_length)
            circle_end = circle_start + circle_length
            line = [chromosome, circle_start, circle_end]
            ## Check if there are N in sequence
            circle_sequence = fasta.fetch(chromosome, circle_start, circle_end)

            if circle_sequence.count('n') > 20 or circle_sequence.count('N') > 20:
                continue

            n_circles += 1
            circle_bed.append(line)
        return circle_bed  

    def run(self):
        chromosomes = self.select_chromosomes()
        circle_bed = self.simulate_coordinates(chromosomes)
        try:
            # Write circles to a bed file 
            with open(self.output_bed, 'w') as file:
                for line in circle_bed:
                    file.write('\t'.join(map(str, line)) + '\n')
        except (FileNotFoundError, PermissionError) as e:
            raise RuntimeError(f"Error writing to BED file '{self.output_bed}': {e}")

if __name__ == "__main__":
    args = parse_coordinates_arguments()
    circles_instance = SimulateCoordinates(args)
    circles_instance.run()