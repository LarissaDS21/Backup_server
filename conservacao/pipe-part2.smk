VCF = open('VCFs/input_cheios_todos_unicos.txt').read().strip().split('\n')

rule all:
    input:
        "Arquivos_finais/TableFrom-VCFs/joined/all.samples.tab"

rule bcftools_annotate:
    input:
        vcf_c="VCFs/Todos/{samples}_corrigido.vcf.gz",
        vcf_c_tbi="VCFs/Todos/{samples}_corrigido.vcf.gz.tbi",
    output:
        vcf_info_selected="Arquivos_finais/All_Vcfs/{samples}.vcf.gz"
    log:
        "Arquivos_finais/All_Vcfs/logs/bcftools_annotate/{samples}.log"
    threads: 1
    resources:
        mem_gb=1
    params:
        jobname = "bcftools_annotate_{samples}"
    shell:
        "/conda/bin/bcftools annotate -x ^INFO/GeneNames,INFO/SequenceOntologyCombined,INFO/TranscriptNameClinicallyRelevant,INFO/HGVScClinicallyRelevant,INFO/HGVSpClinicallyRelevant,"
        "INFO/dbNSFP_gnomAD_exomes_AF,INFO/dbNSFP_1000Gp3_AF,INFO/dbNSFP_ExAC_AF,INFO/dbNSFP_GERP___NR,INFO/dbNSFP_GERP___RS,INFO/dbNSFP_GERP___RS_rankscore,"
        "INFO/PhyloP7wayVertebrate,INFO/PhyloP7wayVertebrateRankscore,INFO/PhyloP20wayMammalian,INFO/PhyloP20wayMammalianRankscore,INFO/MetaSVMScore,INFO/MetaSVMPred,INFO/MetaLRScore {input.vcf_c} -Oz -o {output}"

rule tabix:
    input:
        vcf_gz="Arquivos_finais/All_Vcfs/{samples}.vcf.gz"
    output:
        vcf_gz_tbi="Arquivos_finais/All_Vcfs/{samples}.vcf.gz.tbi"
    log:
        "Arquivos_finais/All_Vcfs/logs/tabix-vcfs/{samples}_vep-vcfs.log"
    threads: 1
    resources:
        mem_gb=1
    params:
        jobname = "tabix_all-vcfs_{samples}"
    shell:
        "tabix -p vcf {input}"

rule bcftools_query:
    input:
        vcf_gz="Arquivos_finais/All_Vcfs/{samples}.vcf.gz",
        vcf_gz_tbi="Arquivos_finais/All_Vcfs/{samples}.vcf.gz.tbi",
    output:
        table="Arquivos_finais/TableFrom-VCFs/{samples}.tab"
    log:
        "Arquivos_finais/All_Vcfs/logs/bcftools_query/{samples}.log"
    threads: 1
    resources:
        mem_gb=1
    params:
        jobname = "bcftools_query_{samples}"
    shell:
        "/conda/bin/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/GeneNames\t%INFO/SequenceOntologyCombined\t%INFO/TranscriptNameClinicallyRelevant\t%INFO/HGVScClinicallyRelevant"
        "\t%INFO/HGVSpClinicallyRelevant\t%INFO/dbNSFP_gnomAD_exomes_AF\t%INFO/dbNSFP_1000Gp3_AF\t%INFO/dbNSFP_ExAC_AF\t%INFO/dbNSFP_GERP___NR"
        "\t%INFO/dbNSFP_GERP___RS\t%INFO/dbNSFP_GERP___RS_rankscore\t%INFO/PhyloP7wayVertebrate\t%INFO/PhyloP7wayVertebrateRankscore\t%INFO/PhyloP20wayMammalian"
	    "\t%INFO/PhyloP20wayMammalianRankscore\t%INFO/MetaSVMScore\t%INFO/MetaSVMPred\t%INFO/MetaLRScore\t[\t%SAMPLE]\n' {input.vcf_gz} > {output}"

rule bcftools_join:
    input:
        table=expand("Arquivos_finais/TableFrom-VCFs/{samples}.tab",samples=VCF)
    output:
        table="Arquivos_finais/TableFrom-VCFs/joined/all.samples.tab"
    log:
        "Arquivos_finais/All_Vcfs/logs/bcftools_join-vcfs/join.log"
    threads: 1
    resources:
        mem_gb=1
    params:
        jobname = "bcftools_join_all_tab_samples"
    shell:
        "echo 'CHROM\tPOS\tID\tREF\tALT\tGeneNames\tSequenceOntologyCombined\tTranscriptNameClinicallyRelevant\tHGVScClinicallyRelevant\tHGVSpClinicallyRelevant\tgnomAD"
        "\t1000G\tExAC\tGERP_NR\tGERP_RS\tGERP_RS_rankscore\tPhyloP7wayVertebrate\tPhyloP7wayVertebrateRankscore\tPhyloP20wayMammalian\tPhyloP20wayMammalianRankscore"
        "\tMetaSVMScore\tMetaSVMPred\tMetaLRScore\t\tSAMPLE' > {output}; "
        "cat {input.table} >> {output}"
