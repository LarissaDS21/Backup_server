VCF = open('VCFs/input-pipe.txt').read().strip().split('\n')
CAMINHO= "VCFs/Todos"

rule all:
    input:
         "VCFs/Todos/CountVariants/VCFs_cheios.txt"

rule sed:
    input:
        vcf="{CAMINHO}/{samples}.vcf.gz"
    output:
        vcf_c="{CAMINHO}/{samples}_corrigido.vcf.gz"
    log:
        "{CAMINHO}/logs/sed_vcf_corrigidos/{samples}.log"
    threads: 1
    resources:
        mem_gb=1
    params:
        jobname = "sed_all_samples_{samples}"
    shell:
        "zcat {input} | sed 's/ /_/g' | bgzip -c > {output}"

rule tabix:
    input:
        vcf_c="{CAMINHO}/{samples}_corrigido.vcf.gz"
    output:
        vcf_c_tbi="{CAMINHO}/{samples}_corrigido.vcf.gz.tbi"
    log:
        "{CAMINHO}/logs/tabix_vcfs/{samples}.log"
    threads: 1
    resources:
        mem_gb=1
    params:
        jobname = "tabix_all_samples_{samples}"
    shell:
        "tabix -p vcf {input}"

rule gatk4:
    input:
        vcf_c="{CAMINHO}/{samples}_corrigido.vcf.gz",
        vcf_c_tbi="{CAMINHO}/{samples}_corrigido.vcf.gz.tbi",
    output:
        variants_txt="{CAMINHO}/CountVariants/{samples}.txt"
    conda:
        "gatk4"
    log:
        "{CAMINHO}/logs/CountVariants/{samples}.log"
    threads: 1
    resources:
        mem_gb=1
    params:
        jobname = "gatk4_all_samples_{samples}"
    shell:
        "gatk CountVariants -V {input.vcf_c} > {output}"


rule vcf_cheio:
    input:
        variants_txt=expand("VCFs/Todos/CountVariants/{samples}.txt", samples=VCF)
    output:
        variants_txt="{CAMINHO}/CountVariants/VCFs_cheios.txt"
    log:
        "{CAMINHO}/logs/vcfs_cheios.log"
    threads: 1
    resources:
        mem_gb=1
    params:
        jobname = "vcf_cheio_all_samples_joincounts"
    shell:
        "for arquivo in {input.variants_txt}; do "
        "nvariants=$(tac $arquivo | head -n1); "
        "if [ $nvariants -gt 0 ]; then "
        "nome_arquivo=$(basename $arquivo | sed 's/.txt//'); "
        "echo $nome_arquivo >> {output}; "
        "elif [ $nvariants -eq 0 ]; then "  
        "nome_arquivo2=$(basename $arquivo | sed 's/.txt//'); "
        "echo $nome_arquivo2 >> VCFs/Todos/CountVariants/VCFs_vazios.txt; fi; "
        "done"
