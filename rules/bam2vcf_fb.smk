
rule intervals:
    """
    Because snakemake input MUST be a file, this hack is necessary to iterate through intervals in the next step,
    inputing one interval in the command at a time. This information is stored in intervals_fb
    """
    output:
        interval = temp(touch(intervalDir + "{i}"))

rule bam2vcf:
    """
    Call freebayes.
    """
    input:
        ref = config['ref'],
        bams = expand(bamDir + "{sample}_dedup.bam", sample=SAMPLES) # input all bams at once
    output: 
        vcf = temp(vcfDir_fb + "{i}.vcf") # output one vcf per interval
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "freebayes -f {input.ref} "
        "-r {wildcards.i} "
        "{input.bams} "
        "> {output.vcf}\n"


rule gatherVcfs:
    """
    This rule gathers all of the VCFs, one per interval, into one final VCF
    """
    input:
        intervals = expand(intervalDir + "{i}", i=intervals_fb),
        vcfs = expand(vcfDir_fb + "{i}.vcf", i=intervals_fb),
        ref = config['ref']
    output: 
        vcf = temp(config["fbDir"] + "Combined.vcf"),
        vcfidx = temp(config["fbDir"] + "Combined.vcf.idx"),
        vcfFiltered = config["fbDir"] + "Combined_hardFiltered.vcf"
    params:
        gatherCommand = myfunc(vcfDir_fb, intervals_fb)
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk GatherVcfs "
        "{params.gatherCommand} "
        "-O {output.vcf}\n"

        "sleep 10\n" # the variant filtration step was failing in an unreproducible way, so added this in case
        # Hard filter Combined.vcf
        "gatk VariantFiltration "
        "-R {input.ref} "
        "-V {output.vcf} "
        "--output {output.vcfFiltered} "
        "--filter-expression \"QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || ExcessHet > 30.0\" "
        "--filter-name \"filteredOut\" "

rule vcftools:
    input:
        vcf = config["fbDir"] + "Combined_hardFiltered.vcf",
        int = "intervals.bed"
    output: 
        missing = config["fbDir"] + "missing_data_per_ind.txt",
        SNPsPerInt = config["fbDir"] + "SNP_per_interval.txt"
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "vcftools --vcf {input.vcf} --remove-filtered-all --minDP 1 --stdout --missing-indv > {output.missing}\n"
        "bedtools intersect -a {input.int} -b {input.vcf} -c > {output.SNPsPerInt}"
