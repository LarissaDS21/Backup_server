configfile: "rMATS_length.yaml"
BAMS = open('bams/StarAlign/input-pipe.txt').read().strip().split('\n')
localrules: echo
EVENTS=('SE',"MXE",...)

rule all:
    input:
        expand("results/{samples}/A3SS.MATS.JC.txt", samples=BAMS)
        
rule echo:
    input:
        path_bam="bams/StarAlign/Casos/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    output: 
        b2_txt="arquivos-Casos-TXT/{samples}_b2.txt"
    shell:
        'echo {input} > {output}'

rule rMATS:
    input:
        b1_txt="controle_b1.txt",
        b2_txt="arquivos-Casos-TXT/{samples}_b2.txt"
    output:
        JC_txt="results/{samples}/A3SS.MATS.JC.txt"
    log:
        "results/tmp_output/{samples}/{samples}_snakemake.log"
    threads: 4
    resources:
        mem_gb=1
    params:
        out_dir = "results/{samples}/",
        tmp_dir = "results/tmp_output/{samples}/",
        length = lambda wc: config['length'][wc.samples],
        jobname = "rMATS_BAM_{samples}"
    shell:
        'rmats.py --b1 {input.b1_txt} --b2 {input.b2_txt} --gtf Homo_sapiens.GRCh38.104.gtf -t paired --variable-read-length --readLength {params.length} --novelSS '
        '--nthread 4 --od {params.out_dir} --tmp {params.tmp_dir} > {log} 2>&1'