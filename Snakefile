import glob


bacteria_name = glob.glob("*.fasta")[0][:-6]
storage_prefix = "."


rule all:
    input: expand("{bacteria}.refseq", bacteria=bacteria_name),


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
         "referenceseeker {input.bacteria_db} {input.genome} > {output}"
