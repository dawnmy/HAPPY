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
            pair_dir + "/{sample}.dedup.frag.filtered.pair.gz", sample=samples),
        mcool = expand(
            results_dir + "/{sample}.dedup.frag.filtered.mcool", sample=samples)

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
        pair = pair_dir + "/{sample}.dedup.frag.pair.gz",
    threads: 20
    conda: condaenv
    params: restrict_py = path.join(cd, "bin/add_restrict.py")
    benchmark: reports_dir + "/benchmarks/{sample}.pair.restrcit.txt"
    shell:
        """
        python {params.restrict_py} -f {input.frag} -t {threads} {input.pair} | pbgzip -n {threads} -c > {output.pair}
        """


# Filter the pairs by pair_type, and the nearest restriction sites
rule pair_filter:
    input: rules.pair_restrict.output.pair
    output: pair_dir + "/{sample}.dedup.frag.filtered.pair.gz"
    threads: 20
    params: restrict_py = path.join(cd, "bin/add_restrict.py")
    benchmark: reports_dir + "/benchmarks/{sample}.pair.filter.txt"
    run:
        import pipes
        import shutil
        import sys
        print(input, str(input).endswith('.gz'))
        print(output, str(output).endswith('.gz'))

        def read_pairs(in_pairs, nproc=1):
            if in_pairs.endswith('.gz'):
                if shutil.which('pbgzip') is None:
                    raise ValueError(
                        'pbgzip is not found, cannot decompress input')
                else:
                    t = pipes.Template()
                    t.append('pbgzip -dc -n {}'.format(nproc), "--")
                    f = t.open(in_pairs, 'r')
                return f
            else:
                return open(in_pairs, 'r')

        def write_pairs(out_pairs, nproc=1):
            if out_pairs.endswith('.gz'):
                if shutil.which('pbgzip') is None:
                    raise ValueError(
                        'pbgzip is not found, cannot decompress input')
                else:
                    t = pipes.Template()
                    t.append('pbgzip -c -n {}'.format(nproc), "--")
                    f = t.open(out_pairs, 'w')
                return f
            else:
                return open(out_pairs, 'w')

        def split_header_body(instream, comment_char='#'):
            header = []
            comment_byte = comment_char.encode()
            # get peekable buffer for the instream
            inbuffer = instream.buffer
            current_peek = inbuffer.peek()
            while current_peek.startswith(comment_byte):
                line = inbuffer.readline()
                header.append(line.decode().strip())
                # peek into the remainder of the instream
                current_peek = inbuffer.peek()

            return header, instream

        def filer_pair(input, output, **kwargs):
            instream = read_pairs(input,
                                  nproc=kwargs.get('threads', 3))

            outstream = write_pairs(output,
                                    nproc=kwargs.get('threads', 3))

            header, body = split_header_body(instream)

            header.append(
                "# The pairs were filtered by 1. not on the same frag; 2. nearst restriction site <= 5kb; \
                    3. MapQ < 10; 4. not uniq map and not one end containing ligation site"
            )
            outstream.writelines((line + '\n' for line in header))
            for line in body:
                cols = line.rstrip().split("\t")
                chrom1, pos1, chrom2, pos2, pair_type = cols[1], int(
                    cols[2]), cols[3], int(cols[4]), cols[7]

                (mapq1, mapq2, idx1, start1, end1, idx2, start2, end2) = map(
                    int, cols[-8:])
                if chrom1 == chrom2 and idx1 == idx2:
                    continue
                elif pos1 - start1 > 5000 and end1 - pos1 > 5000:
                    continue
                elif pos2 - start2 > 5000 and end2 - pos2 > 5000:
                    continue
                elif mapq1 < 10 or mapq2 < 10:
                    continue
                elif pair_type not in ["UU", "UR", "RU"]:
                    continue
                else:
                    outstream.write(line)
            instream.close()
            outstream.close()
        # ! input and output are specific object not str thus need to be converted to str
        filer_pair(str(input), str(output), threads=threads)


# Convert the pair file into .cool, and create index to accelarate the processing
rule cooler:
    input:
        pair = pair_dir + "/{sample}.dedup.frag.filtered.pair.gz",
        chr_sizes = chr_sizes
    output:
        idx = pair_dir + "/{sample}.dedup.frag.filtered.pair.gz.px2",
        cool = results_dir + "/{sample}.dedup.frag.filtered.cool"
    threads: 20
    conda: condaenv
    benchmark: reports_dir + "/benchmarks/{sample}.pair.cooler.txt"
    shell:
        """
        pairix -p pairs {input.pair}
        #cooler cload pairs --assembly hg38 -c1 2 -p1 3 -c2 4 -p2 5 {input.chr_sizes}:1 {input.pair} {output.cool}
        cooler cload pairix --assembly hg38 -p {threads} -s 20 {input.chr_sizes}:1 {input.pair} {output.cool}
        """

# bin the contact into different resolutions and generate .mcool for visualizing
#  in HiGlass
rule bin_cool:
    input:
        rules.cooler.output.cool
    output:
        results_dir + "/{sample}.dedup.frag.filtered.mcool"
    threads: 20
    conda: condaenv
    benchmark: reports_dir + "/benchmarks/{sample}.pair.bin_cool.txt"
    shell:
        """
        cooler zoomify -n {threads} --balance {input} -r 1000,5000,10000,50000,100000 -o {output}
        """
