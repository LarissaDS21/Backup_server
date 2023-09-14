#Alinhamento para análise direcionada ao transcrito
bams = open('barcodes_samples.txt').read().strip().split('\n')

rule all:
    input:
        expand("Alinhamento-Star/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai", samples=bams)


rule bedtools:
    input:
        BAM_ion="IonXpress_{samples}_rawlib.bam"
    output:
        FASTQ="Fastqs/{samples}.fastq"
    log:
        "log/Fastqs/{samples}.log"
    threads: 12
    resources:
        mem_gb=80
    params:
        jobname = "Bedtools_{samples}"
    shell:
        '/conda/bin/bedtools bamtofastq -i {input} -fq {output}'

rule gzip:
    input:
        FASTQ="Fastqs/{samples}.fastq"
    output:
        FASTQ_gz="Fastqs/{samples}.fastq.gz"
    log:
        "log/Fastqs/{samples}_gz.log"
    threads: 12
    resources:
        mem_gb=80
    params:
        jobname = "gzip_{samples}"
    shell:
        'gzip -c {input} > {output}'

rule STAR_align:
    input:
        FASTQ_gz="Fastqs/{samples}.fastq.gz"
    output:
        BAM="Alinhamento-Star/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    log:
        "log/StarAlign/{samples}.log"
    threads: 12
    resources:
        mem_gb=80
    params:
        out_dir = "Alinhamento-Star/{samples}/{samples}_",
        jobname = "Star_Align_{samples}"
    shell:
        'STAR --runThreadN 12 --genomeDir /mnt/gluster01/lgbm/Larissa.Souza/STAR/Genome-index-GRCh38/GenomeDir/ --readFilesIn {input} ' 
        '--readFilesCommand zcat --outFileNamePrefix {params.out_dir} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif '
        '--outFilterMismatchNmax 3 --alignIntronMax 299999 --alignSJDBoverhangMin 6 --alignEndsType EndToEnd --chimSegmentMin 2 --sjdbGTFfile /mnt/gluster01/lgbm/Larissa.Souza/anotacoes_homo_sapiens/Homo_sapiens.GRCh38.104.gtf'

rule samtools:
    input:
        BAM_Star="Alinhamento-Star/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    output:
        BAM_BAI="Alinhamento-Star/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai"
    log:
        "log/samtools/{samples}.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "samtools_{samples}"
    shell:
        '/conda/bin/samtools index -b {input} {output}'