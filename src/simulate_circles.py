import os
from utils import parse_circles_arguments
import numpy as np
from pyfaidx import Fasta

class SimulateCircles:

    def __init__(self, args):
        self.total = args.total
        self.length_min = args.length_min
        self.length_max = args.length_max
        self.genome_fasta = args.genome_fasta
        self.output_bed = args.output_bed

    def select_chromosomes(self):
        chromosomes = {}
        whole_genome_length = 0

        # Filter chromosomes
        with Fasta(self.genome_fasta) as fasta:
            for name, record in fasta.items():
                if all(keyword not in name.lower() for keyword in ['alt', 'random', 'un', 'm']):
                    chromosomes[name] = {'length': len(record), 'weight': 0}
                    
        # Calculate the total genome length
        whole_genome_length = sum(chromosomes[contig]['length'] for contig in chromosomes)
        
        # Calculate and add weights to the chromosomes dictionary
        for contig, info in chromosomes.items():
            info['weight'] = info['length'] / whole_genome_length

        return chromosomes

    def simulate_circles(self, chromosomes):
        circle_bed = []
        n_circles = 0
        
        # Create circles
        while n_circles < self.total:
            n_circles += 1

            chromosome = np.random.choice(list(chromosomes.keys()), p=[chromosomes[contig]['weight'] for contig in chromosomes])
            circle_length = np.random.randint(self.length_min, self.length_max) # CHANGE!!!
            circle_start = np.random.randint(0, (chromosomes[chromosome]['length'] - circle_length))
            circle_end = circle_start + circle_length

            line = [chromosome, circle_start, circle_end]
            circle_bed.append(line)

        try:
            # Write circles to a bed file 
            with open(self.output_bed, 'w') as file:
                for line in circle_bed:
                    file.write('\t'.join(map(str, line)) + '\n')
        except Exception as e:
            raise RuntimeError(f"Error writing to BED file '{self.output_bed}': {e}")

    def run(self):
        chromosomes = self.select_chromosomes()
        self.simulate_circles(chromosomes)

if __name__ == "__main__":
    args = parse_circles_arguments()
    circles_instance = SimulateCircles(args)
    circles_instance.run()
