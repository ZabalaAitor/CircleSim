# CircleSim

**CircleSim** is a command-line tool designed to simulate circular DNA and RNA sequences and generate corresponding sequencing reads. This tool can be helpful for researchers working with circular DNA (like eccDNA) or RNA, providing an easy way to simulate the coordinates of these molecules and generate sequencing files for further analysis.

## Features

- **Simulate Coordinates**: Generates genomic coordinates for circular DNA or RNA molecules with customizable parameters.
- **Simulate Reads**: Creates simulated sequencing reads for the generated circular sequences, including options to introduce mutations and sequencing errors.
- **Merge FASTQ Files**: Merges FASTQ files from different sequencing simulations into a single output.

## Installation

To install CircleSim in a Conda environment, follow these steps:

```bash
git clone https://github.com/yourusername/CircleSim.git
cd CircleSim
conda create -n circlesim python=3.8
conda activate circlesim
pip install -r requirements.txt
```

## Using CircleSim

**Simulate Coordinates:**

```bash
python CircleSim.py coordinates -t DNA -T circular -n 100 -l 500 -L 1500 -o circular_coordinates.bed
python CircleSim.py coordinates -t DNA -T linear -n 100 -l 500 -L 1500 -o linear_coordinates.bed
```

**Simulate reads:**

```bash
python CircleSim.py reads -b circular_coordinates.bed -o circular_reads.fastq
python CircleSim.py reads -b linear_coordinates.bed -o linear_reads.fastq
```

**Merge FASTQ Files:**

```bash
python CircleSim.py  -c circular_reads.fastq -l linear_reads.fastq -o reads.fastq
```

## License

CircleSim is freely available under the [MIT license](./LICENSE).

## Acknowledgements

CircleSim is developed by Aitor Zabala, Iñigo Prada-Luengo, Alex Martínez Ascensión and David Otaegui at Biogipuzkoa Health Research Institute.