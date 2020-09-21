import glob
import os


bacteria_name = glob.glob("*.fasta")[0][:-6]
storage_prefix = "/scratch/mlebedev"

rule all:
    input: "unaligned.fa",
            expand("hists/{sample}.png", sample=bacteria_name),


rule download_db:
    output: targz = expand("{storage_prefix}/downloads/bacteria.tar.gz", storage_prefix=storage_prefix)
    params:
        dir = expand("{storage_prefix}/downloads/", storage_prefix=storage_prefix)
    shell: "wget https://zenodo.org/record/3725706/files/bacteria.tar.gz -P {params.dir}"


rule unzip:
    input: "{storage_prefix}/downloads/bacteria.tar.gz"
    output: "{storage_prefix}/bacteria/"
    shell: "mkdir {output} | tar -zxf {input} --directory {output}"


rule referenceseeker:
    input:
         genome = "{bacteria}.fasta",
         bacteria_db = expand("{storage_prefix}/bacteria/", storage_prefix=storage_prefix)
    output:
        "{bacteria}.refseq"
    conda:
        "envs/referenceseeker.yml"
    shell:
         "referenceseeker --threads 10 {input.bacteria_db} {input.genome} > {output}"


rule path_to_relatives:
    input:
         reference = expand("{bacteria}.refseq", bacteria=bacteria_name),
         bacteria_db = expand("{storage_prefix}/bacteria/", storage_prefix=storage_prefix)
    output:
         "path_to_refs.txt"
    run:
        with open(f"{input.reference}") as ref:
            genomes = [line.strip().split()[0] for line in ref.readlines()[1::]]
        with open(f"{output}", "w") as out:
            for i in genomes:
                out.write(f"{input.bacteria_db}{i}.fna\n")


rule cp_relatives:
    input: "path_to_refs.txt"
    output: "references/"
    shell:
        'mkdir -p references && '
        'cp `cat {input}` references && '
        'mv references/* {output}'


IDS = glob_wildcards("references/{id}.fna")

## TODO: messed up the wildcards
rule mashmap:
    input:
        reference = expand("references/{ref_genomes}.fna", ref_genomes=IDS),
        genome = expand("{bacteria}.fasta", bacteria=bacteria_name)
    conda:
        "envs/mashmap.yml"
    threads: 10
    output:
        expand("mashmap/{bacteria}_{ref_genomes}_out.mashmap", ref_genomes=IDS, bacteria=bacteria_name)
    shell:
        'mkdir -p mashmap &&'
        'mashmap -r {input.reference} -q {input.genome} -s 500b -t {threads} --perc_identity 95 -o mashmap/{output}'


rule get_nonaligned:
    input:
        mashmaps = expand("mashmap/{bacteria}_{ref_genomes}_out.mashmap", ref_genomes=IDS, bacteria=bacteria_name)
    params:
        mashmaps_dir = "mashmap/"
    output: "instruction_to_cut.bed"
    run:
        for root, dirs, files in os.walk(f"{params.mashmaps_dir}"):
            for file in files:
                ranges = []
                name = None
                with open(file) as m:
                    for line in m.readlines():
                        ranges.append([int(i) for i in line.split()[2:4]])
                        name = line.split()[0]

                ranges.sort(key=lambda x: x[0])

                for num, i in enumerate(ranges[1::]):
                    num_to_append = ranges[num][1]
                    if i[0] - ranges[num][1] > 1:
                        print(ranges[num][1], i[0])
                        with open(f"{output}", "a") as out:
                            out.write("\t".join((name, str(ranges[num][1]), str(i[0]))) +  "\n")


rule bedtools:
    input: cut = "instruction_to_cut.bed",
           genome = expand("{bacteria}.fasta", bacteria=bacteria_name)
    output: "unaligned.fa"
    shell: "bedtools getfasta -fi {input.genome} -bed {input.cut} -fo {output}"




rule prodigal:
	input: expand("{bacteria}.fasta", bacteria=bacteria_name)
	output: "proteins/{bacteria}.faa"
	shell: "prodigal -a {output} -q -i {input}"


rule diamond:
	input: "proteins/{bacteria}.faa"
	output: "lengths/{bacteria}.data"
	threads: 10
	params:
		db="/scratch/mlebedev/uniprot_trembl.diamond.dmnd",
		of="6 qlen slen"
	shell: "diamond blastp --threads {threads} --max-target-seqs 1 --db {params.db} --query {input} --outfmt {params.of} --out {output}"


rule hist:
	input: "lengths/{bacteria}.data"
	output: "hists/{bacteria}.png"
	shell: "Rscript ../ideel/scripts/hist.R {input} {output}"

## TODO Change to make it work with multiple genomes
rule prokka:
    input: expand("{bacteria}.fasta", bacteria=bacteria_name)
    output: "prokka/{bacteria_name}.gff"
    threads: 4
    shell: "prokka --outdir {output} --prefix {bacteria_name} {input} -cpus {threads}"


rule roary:
    input: expand("prokka/{bacteria}.gff", bacteria=bacteria_name)
    output: "roary/"
    threads: 4
    shell: "roary -p 4 -f {output}"
