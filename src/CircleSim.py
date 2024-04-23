import argparse
import os
from reads import SimulateReads, read_bed_file
from coordinates import SimulateCoordinates
from join import Join
from utils import parse_coordinates_arguments, parse_reads_arguments, parse_join_arguments, set_default_genome_fasta


class CircleSim:

    def __init__(self):
        self.cs_version = "1.0"

        self.parser = argparse.ArgumentParser(
            description='CircleSim',
            usage=f'''CircleSim <subprograms> [options]

version = {self.cs_version}

CircleSim

Commands:

    coordinates  Simulate coordinates
    reads   Simulate reads 
    join   Join fastq files 
'''
        )
        subparsers = self.parser.add_subparsers(dest='subcommand')

        self.coordinates = subparsers.add_parser(
            name='coordinates',
            description='Simulate coordinates',
            prog='CircleSim coordinates',
            usage='''CircleSim coordinates [options]'''
        )
        parse_coordinates_arguments(self.coordinates)

        self.reads = subparsers.add_parser(
            name='reads',
            description='Simulate circle reads',
            prog='CircleSim reads',
            usage='''CircleSim reads [options]'''
        )
        parse_reads_arguments(self.reads)

        self.join = subparsers.add_parser(
            name='join',
            description='Join fastq files',
            prog='CircleSim join',
            usage='''CircleSim join [options]'''
        )
        parse_join_arguments(self.join)

    def main(self):
        args = self.parser.parse_args()

        if args.subcommand == 'coordinates':
            set_default_genome_fasta(args)
            circles_instance = SimulateCoordinates(args)
            circles_instance.run()
        elif args.subcommand == 'reads':
            set_default_genome_fasta(args)
            circle_bed = read_bed_file(args.input_bed)
            reads_instance = SimulateReads(args, circle_bed)
            reads_instance.simulate_reads()
        elif args.subcommand == 'join':
            join_instance = Join(args)
            join_instance.join()

if __name__ == "__main__":
    circlesim = CircleSim()
    circlesim.main()
