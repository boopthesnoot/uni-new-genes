import glob


bacteria_name = glob.glob("*.fasta")[0][:-6]
storage_prefix = "/scratch/mlebedev"

rule all:
    input: "unaligned.fa",
           "viz.done",
           f"hists/{bacteria_name}.png",


rule download_db:
    output: targz = expand("{storage_prefix}/downloads/bacteria.tar.gz", storage_prefix=storage_prefix)
    params:
        dir = storage_prefix + "/downloads/"
    shell: "wget https://zenodo.org/record/3725706/files/bacteria.tar.gz -P {params.dir}"


rule unzip:
    input: "{storage_prefix}/downloads/bacteria.tar.gz"
    output: "{storage_prefix}/bacteria/"
    shell: "mkdir -p {output} && tar -zxf {input} --directory {output}"


rule referenceseeker:
    input:
         genome = "{bacteria}.fasta",
         bacteria_db = storage_prefix + "/bacteria/"
    output:
        "{bacteria}.refseq"
    conda:
        "envs/referenceseeker.yml"
    shell:
         "referenceseeker --threads 10 {input.bacteria_db} {input.genome} > {output}"


rule path_to_relatives:
    input:
         reference = bacteria_name + ".refseq",
         bacteria_db = storage_prefix + "/bacteria/"
    output:
         "path_to_refs.txt"
    run:
        with open(f"{input.reference}") as ref:
            genomes = [line.strip().split()[0] for line in ref.readlines()[1::]]
        with open(f"{output}", "w") as out:
            for i in genomes:
                out.write(f"{input.bacteria_db}{i}.fna\n")


rule mashmap:
    input:
        reference = "path_to_refs.txt",
        genome = "{bacteria_name}.fasta"
    conda:
        "envs/mashmap.yml"
    threads: 10
    output:
        "mashmap/{bacteria_name}_out.mashmap"
    shell:
        '''
        mashmap --rl {input.reference} -q {input.genome} -s 500 -t {threads} --perc_identity 95 -o {output}
        '''


rule get_nonaligned:
    input: f"mashmap/{bacteria_name}_out.mashmap"
    output: "instruction_to_cut.bed"
    run:
        ranges = []
        name = None
        with open(f"{input}") as m:
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
		db = "/scratch/mlebedev/uniprot_trembl.diamond.dmnd",
		of = "6 qlen slen"
	shell:
            '''
            diamond blastp --threads {threads} --max-target-seqs 1 --db {params.db} \
            --query {input} --outfmt {params.of} --out {output}
            '''


rule hist:
	input: "lengths/{bacteria}.data"
	output: "hists/{bacteria}.png"
	shell: "Rscript ../ideel/scripts/hist.R {input} {output}"


rule make_prokka_input:
    input: refs="path_to_refs.txt"
    params: bacteria_name
    output: directory("prokka_inp/")
    shell:
        '''
        mkdir -p {output} && 
        cp `cat {input.refs}` {output} && 
        cp {params}.fasta {output}
        mv {output}{params}.fasta {output}{params}.fna
        '''


rule prokka:
    input: "prokka_inp/"
    output: directory("prokka/")
    threads: 10
    shell: '''
            for file in $(ls {input} | rev | cut -c5- | rev)
            do prokka --outdir prokka --prefix $file {input}$file.fna -cpus {threads} --force
            done
           '''

IDS_roary_in = glob_wildcards("prokka/{id}.gff")

rule roary:
    input: expand("prokka/{id}.gff", id=IDS_roary_in)
    output: directory("roary_out/")
    params: "prokka/*.gff"
    threads: 10
    shell: "roary -p {threads} -e -n -v -f {output} {params}"

IDS_roary_out = glob_wildcards("roary_out/{dirname}/core_gene_alignment.aln")[0]

rule fast_tree:
    input:
         f"roary_out/{IDS_roary_out[0]}/core_gene_alignment.aln"
    output:
         f"roary_out/{IDS_roary_out[0]}/core_gene_alignment.newick"
    shell:
         "FastTree -nt -gtr {input} > {output}"

rule roary_grahs:
    input:
         tree=f"roary_out/{IDS_roary_out[0]}/core_gene_alignment.newick",
         csv=f"roary_out/{IDS_roary_out[0]}/gene_presence_absence.csv"
    output: "viz.done"
    shell: "python roary_plots/roary_plots.py {input.tree} {input.csv} && touch {output}"
