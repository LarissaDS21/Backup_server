fastqs = open('/mnt/data4/lgbm/Larissa.Souza/PSI/rMATS/fastqs/SRA/SRR-samples.txt').read().strip().split('\n')


rule all:
    input:
        expand("bams/StarAlign/CTRLs/{samples}/{samples}_Aligned.sortedByCoord.out.bam", samples=fastqs)
        

rule STAR_align:
    input:
        FASTQ_R1="fastqs/SRA/{samples}_1.fastq.gz",
        FASTQ_R2="fastqs/SRA/{samples}_2.fastq.gz"
    output:
        BAM="bams/StarAlign/CTRLs/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    log:
        "bams/StarAlign/tmp_output/{samples}.log"
    threads: 12
    resources:
        mem_gb=80
    params:
        out_dir = "bams/StarAlign/CTRLs/{samples}/{samples}_",
        jobname = "Star_Align_{samples}"
    shell:
        'STAR --runThreadN 12 --genomeDir /mnt/data4/lgbm/Larissa.Souza/STAR/Genome-index-GRCh38/GenomeDir/ --readFilesIn {input.FASTQ_R1} {input.FASTQ_R2} ' 
        '--readFilesCommand zcat --outFileNamePrefix {params.out_dir} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif '
        '--outFilterMismatchNmax 3 --alignIntronMax 299999 --alignSJDBoverhangMin 6 --alignEndsType EndToEnd --chimSegmentMin 2 --sjdbGTFfile Homo_sapiens.GRCh38.104.gtf'
        
