import os

configfile:
    "configs/config.yaml"

def get_expected_output(reads_dir):
    return(list(map(lambda file: file.split(".")[0], os.listdir(reads_dir))))

genome = config["reference"]
index_path = genome.rsplit(".", maxsplit=1)[0]
iter_num = config["output_prefix"]
start, stop = config["extraction_size"][0], config["extraction_size"][1]+1
reads_prefixes = get_expected_output(config["trimmed_samples"])

rule output:
    input:
        expand("results/"+iter_num+"/original/{dataset_id}_mapped_{n}.pdf",
                    dataset_id=reads_prefixes,
                    n=[i for i in range(start, stop)]),
        expand("results/"+iter_num+"/fastq/{dataset_id}_{type}.fq", dataset_id=reads_prefixes, type=["mapped", "unmapped"])

rule bowtie_build:
    input:
        genome
    output:
        expand(index_path+".{index}.ebwt", index=range(1,5)),
        expand(index_path+".rev.{index}.ebwt", index=range(1,3))
    params:
        index=index_path
    shell:
        "bowtie-build {input} {params.index}"

rule bowtie_map_raw:
    input:
        reads=config["trimmed_samples"]+"{dataset_id}.fq",
        indexes=rules.bowtie_build.output
    output:
        "results/"+iter_num+"/alignments/{dataset_id}.bam"
    params:
        index=index_path
    threads:
        32
    shell:
        "bowtie -p{threads} {params.index} {input.reads} -S -n 0 | samtools view - -Sb | samtools sort - -@{threads} -o {output}"

rule filer_mapped_bam:
    input:
        rules.bowtie_map_raw.output
    output:
        "results/"+iter_num+"/alignments/{dataset_id}_mapped.bam"
    shell:
        "samtools view -b -F 4 {input} > {output}"

rule filter_unmapped_bam:
    input:
        rules.bowtie_map_raw.output
    output:
        "results/"+iter_num+"/alignments/{dataset_id}_unmapped.bam"
    shell:
        "samtools view -b -f 4 {input} > {output}"

rule extract_fastq_from_bam:
    input:
        "results/"+iter_num+"/alignments/{dataset_id}_{type}.bam"
    output:
        "results/"+iter_num+"/fastq/{dataset_id}_{type}.fq"
    shell:
        "samtools fastq {input} > {output}"

rule extract_original_fastq:
    input:
        original=config["original_samples"]+"{dataset_id}.fq",
        trimmed="results/"+iter_num+"/fastq/{dataset_id}_mapped.fq"
    output:
        "results/"+iter_num+"/original/{dataset_id}_mapped.fq"
    shell:
        "seqtk subseq {input.original} <(awk 'NR % 4 == 1' {input.trimmed} | sed 's/@//') > {output}"

rule split_original_reads_by_length:
    input:
        rules.extract_original_fastq.output
    output:
        "results/"+iter_num+"/original/{dataset_id}_mapped_{n}.fa"
    shell:
        "reformat.sh in={input} out={output} minlength={wildcards.n} maxlength={wildcards.n}"

rule create_logo:
    input:
        rules.split_original_reads_by_length.output
    output:
        "results/"+iter_num+"/original/{dataset_id}_mapped_{n}.pdf"
    shell:
        """
        weblogo -S 1.0 \
            -s large \
            --ticmarks 0.5 \
            --errorbar-gray 0.1 \
            -y bits \
            --number-fontsize 17 \
            --fontsize 21 \
            --small-fontsize 1 \
            -c classic \
            --text-font Helvetica \
            --logo-font Helvetica \
            --title-font Helvetica \
            --number-interval 5 \
            --stack-width 12 \
            --resolution 600 \
            --aspect-ratio 5 \
            --format pdf \
            --fin {input} \
            --fout {output}
        """
