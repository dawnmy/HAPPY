# HAPPY
Hi-C data Analysis and Processing PIpeline (HAPPY). This software outputs the contact matrix with different bin sizes (resolutions): 1k, 5k, 10k, 50k, 100k.

## Requirement

1. conda
2. snakemake
3. BWA-MEM
4. pairtools
5. cooler
6. SAMtools

`conda` and `snakemake` need to be installed manually and the rest will be automatically installed during the first launch of the program.

The output pair file can be easily filered based on provided condition.

## Installation

```shell
git clone git@github.com:dawnmy/HAPPY.git
```

## Usage

1. Adapt the config file for the pipeline
Modify the `config/config.yaml` file in the program folder to adapt to your data location.

2. Download the reference genome and create BWA index and fasta index (`.fai`). 
For instance GRCh38 for homo sapiens:
```shell
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O GRCh38.fasta.gz

pigz -d GRCh38.fasta.gz
bwa index GRCh38.fasta
samtools fai GRCh38.fasta
```

3. Make the chromosome sizes file based on the `.fai`
```shell
cut -f1,2 GRCh38.fasta.fai > chrom.all.sizes
```

4. Launch the pipeline
With 20 threads
```shell
snakemake -s runHiC.smk --use-conda -j 20
```

If you use SGE submission system:
```shell
snakemake -s runHiC.smk --use-conda -c "qsub -cwd -pe multislot {threads} -i /dev/null -v PATH" -j 2
```

[![Analytics](https://ga-beacon.appspot.com/UA-143780277-1/HAPPY/readme?pixel)](https://github.com/igrigorik/ga-beacon)
