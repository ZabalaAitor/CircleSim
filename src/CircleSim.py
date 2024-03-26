import argparse
import os
import sys
from simulate_circles import SimulateCircles
from simulate_reads import SimulateReads
from simulate_reads import read_bed_file
from utils import parse_circles_arguments, parse_reads_arguments


class CircleSim:

    def __init__(self):
        self.cs_version = "1.0"

        self.parser = argparse.ArgumentParser(
            description='CircSim',
            usage=f'''CircSim <subprograms> [options]

version = {self.cs_version}

CircleSim

Commands:

    circles  Simulate circles
    reads   Simulate circle reads 
'''
        )
        subparsers = self.parser.add_subparsers(dest='subcommand')
        self.circles = subparsers.add_parser(
            name='circles',
            description='Simulate circles',
            prog='CircSim circles',
            usage='''CircSim circles [options]'''
        )
        self.circles.add_argument('-t', '--type', type=str, required=True, choices=['DNA', 'RNA'], help='Type of circle (DNA or RNA)')
        self.circles.add_argument('-n', '--number', type=int, default=100, help='Number of circles to simulate')
        self.circles.add_argument('-d', '--distribution', type=str, choices=['uniform', 'lognormal'], default='uniform', help='Distribution of circle length (uniform or lognormal)')
        self.circles.add_argument('-l', '--length_min', type=int, default=300, help='Minimum length of eccDNA')
        self.circles.add_argument('-L', '--length_max', type=int, default=2000, help='Maximum length of eccDNA')
        self.circles.add_argument('-m', '--mean', type=int, default=1000, help='Circle length mean')
        self.circles.add_argument('-sd', '--sd', type=int, default=100, help='Circle length standard deviation')
        self.circles.add_argument('-g', '--genome_fasta', type=str, help='Path to the genome FASTA file')
        self.circles.add_argument('-o', '--output_bed', type=str, default=os.getcwd(), help='Path to the output BED file for eccDNA')

        self.reads = subparsers.add_parser(
            name='reads',
            description='Simulate circle reads',
            prog='CircSim reads',
            usage='''CircSim reads [options]'''
        )
        self.reads.add_argument('-t', '--type', type=str, required=True, choices=['DNA', 'RNA'], help='Type of circle (DNA or RNA)')
        self.reads.add_argument('-c', '--coverage', type=float, default=30, help='Coverage for reads simulation')
        self.reads.add_argument('-r', '--reads_length', type=int, default=150, help='Length of simulated reads')
        self.reads.add_argument('-i', '--insert_length', type=int, default=500, help='Insert length for reads simulation')
        self.reads.add_argument('-g', '--genome_fasta', type=str, required=True, help='Path to the genome FASTA file')
        self.reads.add_argument('-b', '--input_bed', type=str, help='Path to the input BED file for eccDNA')
        self.reads.add_argument('-o', '--output_fastq', type=str, default=os.getcwd(), help='Path to the output FASTQ file for reads')

    def args_circles(self):
        return self.parser.parse_args()

    def args_reads(self):
        return self.parser.parse_args()
        
    def main(self):
        args = self.parser.parse_args()

        if args.subcommand == 'circles':
            # Check if genome FASTA file is not provided
            if not args.genome_fasta:
                # Set default genome FASTA file based on the circle type
                if args.circle == 'DNA':
                    args.genome_fasta = '' ##########################
                elif args.circle == 'RNA':
                    args.genome_fasta = '' ##########################
                # Check if default genome FASTA file exists locally
                if os.path.exists(default_genome_fasta):
                    args.genome_fasta = default_genome_fasta
                else:
                    print("Default genome FASTA file not found locally. Downloading from GitHub...")
                    # Download default genome FASTA file
                    subprocess.run(['wget', default_genome_fasta])
                    args.genome_fasta = default_genome_fasta

            circles_instance = SimulateCircles(self.args_circles())
            circles_instance.run()
        elif args.subcommand == 'reads':
            # Check if genome FASTA file is not provided
            if not args.genome_fasta:
                # Set default genome FASTA file based on the circle type
                if args.circle == 'DNA':
                    args.genome_fasta = '' ##########################
                elif args.circle == 'RNA':
                    args.genome_fasta = '' ##########################
                # Check if default genome FASTA file exists locally
                if os.path.exists(default_genome_fasta):
                    args.genome_fasta = default_genome_fasta
                else:
                    print("Default genome FASTA file not found locally. Downloading from GitHub...")
                    # Download default genome FASTA file
                    subprocess.run(['wget', default_genome_fasta])
                    args.genome_fasta = default_genome_fasta

            circle_bed = read_bed_file(args.input_bed)
            reads_instance = SimulateReads(args, circle_bed)
            reads_instance.simulate_reads()


if __name__ == "__main__":
    circsim = CircSim()
    circsim.main()

