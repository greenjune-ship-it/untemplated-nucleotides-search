import os

configfile:
    "config/config.yaml"

def get_expected_output(reads_dir):
    return(list(map(lambda file: file.split(".")[0], os.listdir(reads_dir))))


def get_bowtie_map_input(include_cutadapt):
    if include_cutadapt:
        return(rules.cut_adapt.output)
    else:
        return(path_to_samples + "{dataset_id}.fq")


genome = config["reference"]
index_path = genome.rsplit(".", maxsplit=1)[0]
iter_num = config["output_prefix"]
start, stop = config["extraction_size"][0], config["extraction_size"][1] + 1
samples = get_expected_output(config["trimmed_samples"])
path_to_samples = config["trimmed_samples"]
include_cutadapt = config["include_cutadapt"]
if include_cutadapt:
    path_to_samples = "results/" + iter_num + "/trimmed_reads/"

rule output:
    input:
        expand(
            "results/" + iter_num + "/original/mapped/{dataset_id}_{n}.pdf",
            dataset_id=samples,
            n=[i for i in range(start, stop)],
            mapping_status=["mapped", "unmapped"]
        ),
        expand(
            "results/" + iter_num + "/fastq_from_alignments/{mapping_status}/{dataset_id}.fq",
            dataset_id=samples,
            mapping_status=["mapped", "unmapped"]
        )

rule cut_adapt:
    input:
        config["trimmed_samples"] + "{dataset_id}.fq"
    output:
        path_to_samples + "{dataset_id}_trimmed.fq"
    run:
        if include_cutadapt:
            shell("cutadapt --cores=32 -u -1 -o {output} {input}")

rule bowtie_build:
    input:
        genome
    output:
        expand(index_path + ".{index}.ebwt", index=range(1,5)),
        expand(index_path + ".rev.{index}.ebwt", index=range(1,3))
    params:
        index=index_path
    shell:
        "bowtie-build {input} {params.index}"

rule bowtie_map_raw:
    input:
        reads=get_bowtie_map_input(include_cutadapt),
        indexes=rules.bowtie_build.output
    output:
        "results/" + iter_num + "/alignments/{dataset_id}.bam"
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
        "results/" + iter_num + "/alignments/mapped/{dataset_id}.bam"
    shell:
        "samtools view -b -F 4 {input} > {output}"

rule filter_unmapped_bam:
    input:
        rules.bowtie_map_raw.output
    output:
        "results/" + iter_num + "/alignments/unmapped/{dataset_id}.bam"
    shell:
        "samtools view -b -f 4 {input} > {output}"

rule extract_fastq_from_bam:
    input:
        "results/" + iter_num + "/alignments/{mapping_status}/{dataset_id}.bam"
    output:
        "results/" + iter_num + "/fastq_from_alignments/{mapping_status}/{dataset_id}.fq"
    shell:
        "samtools fastq {input} > {output}"

rule extract_original_fastq:
    input:
        original=config["original_samples"] + "{dataset_id}.fq",
        trimmed="results/" + iter_num + "/fastq_from_alignments/{mapping_status}/{dataset_id}.fq"
    output:
        "results/" + iter_num + "/original/{mapping_status}/{dataset_id}.fq"
    shell:
        "seqtk subseq {input.original} <(awk 'NR % 4 == 1' {input.trimmed} | sed 's/@//') > {output}"

rule split_original_reads_by_length:
    input:
        rules.extract_original_fastq.output
    output:
        "results/" + iter_num + "/original/{mapping_status}/{dataset_id}_{n}.fa"
    shell:
        "reformat.sh in={input} out={output} minlength={wildcards.n} maxlength={wildcards.n}"

rule create_logo:
    input:
        rules.split_original_reads_by_length.output
    output:
        "results/" + iter_num + "/original/{mapping_status}/{dataset_id}_{n}.pdf"
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
