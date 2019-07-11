from os import path

configfile: "config/config.yaml"
condaenv = "config/conda_hic.yaml"

# Define paths for input and output
fq_dir = config["fq_dir"]
out_dir = config["out_dir"]

ref_genome = config["ref_genome"]

# Create the size file with fai file
chr_sizes = config["chr_sizes"]
restrict_enzyme = config["restrict_enzyme"]

qc_dir = path.join(out_dir, "qc_fq")
bam_dir = path.join(out_dir, "bam")
pair_dir = path.join(out_dir, "pair")
reports_dir = path.join(out_dir, "reports")
results_dir = path.join(out_dir, "results")
cd = path.dirname(path.abspath(__file__))

samples, = glob_wildcards(fq_dir + "/{sample}_1.fastq.gz")

# Define the wildcard for sample name
wildcard_constraints:
    sample = "[^\.\/]+"


# Per default, the reads files end with _1.fastq.gz and _2.fastq.gz
def get_fq(wc):
    return [path.join(fq_dir, wc.sample + end + ".fastq.gz") for end in ["_1", "_2"]]


rule all:
    input:
        pair = expand(
            pair_dir + "/{sample}.dedup.frag.noseq.pair.gz", sample=samples)

# pairtools (v0.3.0) should be installed manualy instead of using conda (v0.2.2)
rule install_pairtools:
    input:
        fq_dir + samples[0] + "_1.fastq.gz"
    output:
        touch("pairtools.installed")
    conda: condaenv
    shell:
        """
        conda uninstall pairtools -y
        wget https://github.com/mirnylab/pairtools/archive/v0.3.0.tar.gz -O pairtools-0.3.0.tar.gz
        tar zvxf pairtools-0.3.0.tar.gz
        pip uninstall pairtools -y
        python pairtools-0.3.0/setup.py install
        """

# Quality control
rule fastp:
    input:
        get_fq
    output:
        or1 = qc_dir + "/{sample}.qc_1.fq",
        or2 = qc_dir + "/{sample}.qc_2.fq",
        html = reports_dir + "/fastp/{sample}.qc.report.html"
    threads: 8
    conda: condaenv
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output.or1} -O {output.or2} \
            -5 20 -3 20 -l 50 -h {output.html} -w {threads}
        """

# Alignment, map two pairs independently and keep 5' when there is a cliping
rule bwa:
    input:
        r1 = rules.fastp.output.or1,
        r2 = rules.fastp.output.or2,
        ref = ref_genome
    output: bam_dir + "/{sample}.bwa.bam",
    threads: 20
    conda: condaenv
    benchmark: reports_dir + "/benchmarks/{sample}.bwa.txt"
    shell:
        """
        bwa mem -SP5M -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -Shb -@ {threads} \
            -q 10 -o {output}
        """

# Convert bam to pairs format
rule pair_sort:
    input:
        bam = rules.bwa.output,
        chr_sizes = chr_sizes,
        pairtools_installed = rules.install_pairtools.output
    output: pair_dir + "/{sample}.noseq.pair.gz"
    threads: 20
    conda: condaenv
    benchmark: reports_dir + "/benchmarks/{sample}.pairsort.txt"
    shell:
        """
        pairtools parse -c {input.chr_sizes} --assembly hg38 --nproc-in {threads} --drop-seq --drop-sam \
            --add-columns mapq {input.bam} | pairtools sort --nproc {threads} --nproc-out {threads} -o {output}
        """


# rule pair_merge:

# Generate the restriction fragment location
rule make_frag:
    input:
        ref = ref_genome,
        chr_sizes = chr_sizes
    output: ref_genome + ".digested"
    params: enzyme = restrict_enzyme
    conda: condaenv
    shell:
        """
        cooler digest {input.chr_sizes} {input.ref} {params.enzyme} -o {output}
        """

# Remove pairs from PCR duplicates
rule pair_dedup:
    input:
        pair = rules.pair_sort.output
    output:
        dedup_pair = pair_dir + "/{sample}.dedup.noseq.pair.gz",
        dedup_stats = pair_dir + "/{sample}.dedup.noseq.pair.stats"
    threads: 20
    conda: condaenv
    benchmark: reports_dir + "/benchmarks/{sample}.pairdedup.txt"
    shell:
        """
        pairtools dedup --mark-dups --output {output.dedup_pair} --output-stats {output.dedup_stats} \
            {input.pair}
        """


# Add the restriction fragment information in the pair file
rule pair_restrict:
    input:
        pair = rules.pair_dedup.output.dedup_pair,
        frag = rules.make_frag.output
    output:
        pair = pair_dir + "/{sample}.dedup.frag.noseq.pair.gz",
        idx = pair_dir + "/{sample}.dedup.frag.noseq.pair.gz.px2"
    threads: 20
    conda: condaenv
    params: restrict_py = path.join(cd, "bin/add_restrict.py")
    benchmark: reports_dir + "/benchmarks/{sample}.pair.restrcit.txt"
    shell:
        """
        python {params.restrict_py} -f {input.frag} -t {threads} {input.pair} | pbgzip -n {threads} -c > {output.pair}
        pairix -p pairs {output.pair}
        """

# TODO

# Filter the pairs by pair_type, and the nearest restriction sites
# rule pair_filter:

# Convert the pair file into .cool or .mcool files with different resolutions and
# visualize in HiGlass
# rule cool:
