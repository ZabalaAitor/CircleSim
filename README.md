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

## Usage

CircleSim provides three main commands:

### 1. Simulate Coordinates

```bash
python circlesim.py coordinates [options]
```

### 2. Simulate Reads

```bash
python circlesim.py reads [options]
```

### 2. Merge FASTQ Files

```bash
python circlesim.py join[options]
```

## Examples

**Simulate Coordinates:**

```bash
python circlesim.py coordinates -t DNA -T circular -n 100 -d uniform -l 500 -L 1500 -o output_coordinates.bed
```

**Simulate reads:**

```bash
python circlesim.py reads -l 150 -e 0.01 -m 0.001 -b output_coordinates.bed -o simulated_reads.fastq
```

**Merge FASTQ Files:**

```bash
python circlesim.py merge_fastq -i read1.fastq,read2.fastq -o merged_reads.fastq
```