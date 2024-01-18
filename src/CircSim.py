import argparse
import os
import sys
from simulate_eccDNA import SimulateEccDNA
from simulate_reads import SimulateReads
from simulate_reads import read_bed_file
from utils import parse_eccDNA_arguments, parse_reads_arguments


class CircSim:

    def __init__(self):
        self.cs_version = "1.0"

        self.parser = argparse.ArgumentParser(
            description='CircSim',
            usage=f'''CircSim <subprograms> [options]

version = {self.cs_version}

CircSim

Commands:

    eccDNA  Simulate eccDNA
    reads   Simulate eccDNA reads
'''
        )
        subparsers = self.parser.add_subparsers(dest='subcommand')

        self.eccDNA = subparsers.add_parser(
            name='eccDNA',
            description='Simulate eccDNA',
            prog='CircSim eccDNA',
            usage='''CircSim eccDNA [options]'''
        )
        self.eccDNA.add_argument('-t', '--total_eccDNA', type=int, default=100, help='Total number of eccDNA to simulate')
        self.eccDNA.add_argument('-l', '--length_min', type=int, default=300, help='Minimum length of eccDNA')
        self.eccDNA.add_argument('-L', '--length_max', type=int, default=2000, help='Maximum length of eccDNA')
        self.eccDNA.add_argument('-g', '--genome_fasta', type=str, required=True, help='Path to the genome FASTA file')
        self.eccDNA.add_argument('-o', '--output_bed', type=str, default=os.getcwd(), help='Path to the output BED file for eccDNA')

        self.reads = subparsers.add_parser(
            name='reads',
            description='Simulate eccDNA reads',
            prog='CircSim reads',
            usage='''CircSim reads [options]'''
        )
        self.reads.add_argument('-s', '--sequence', type=str, default='short', help='Coverage for reads simulation')
        self.reads.add_argument('-c', '--coverage', type=float, default=30, help='Coverage for reads simulation')
        self.reads.add_argument('-r', '--reads_length', type=int, default=150, help='Length of simulated reads')
        self.reads.add_argument('-i', '--insert_length', type=int, default=500, help='Insert length for reads simulation')
        self.reads.add_argument('-g', '--genome_fasta', type=str, required=True, help='Path to the genome FASTA file')
        self.reads.add_argument('-b', '--input_bed', type=str, required=True, help='Path to the input BED file for eccDNA')
        self.reads.add_argument('-o', '--output_fastq', type=str, default=os.getcwd(), help='Path to the output FASTQ file for reads')

    def args_eccDNA(self):
        return self.parser.parse_args()

    def args_reads(self):
        return self.parser.parse_args()
        
    def main(self):
        args = self.parser.parse_args()

        if args.subcommand == 'eccDNA':
            eccdna_instance = SimulateEccDNA(self.args_eccDNA())
            eccdna_instance.run()
        elif args.subcommand == 'reads':
            circle_bed = read_bed_file(args.input_bed)
            reads_instance = SimulateReads(args, circle_bed)
            reads_instance.simulate_reads()


if __name__ == "__main__":
    circsim = CircSim()
    circsim.main()

