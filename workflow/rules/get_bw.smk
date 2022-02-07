## convert the bam file to bigwig format (signal normalized to RPM) for viewing the data on Genome Browse
## get rpm_factor
factors = pd.read_table('report/rpm_factor.txt', header=None)
sample_factor_dict = dict(zip(factors.iloc[:,0], factors.iloc[:,1]))

rule get_bw:
    input: 
        "STAR_align/{sample}.bam"
    output: 
        "bw_rpm/{sample}.bw"
    threads: 1
    params:
        mem = '10G',
        rpmFactor = lambda wildcards: sample_factor_dict[wildcards.sample],
        jobName = "get_bw.{sample}"
    conda: "../envs/rnaseq_env.yaml"
    shell:
        #"echo {params.rpmFactor} && "
        "genomeCoverageBed -split -bg -ibam {input} -scale {params.rpmFactor} > bw_rpm/{wildcards.sample}.bg && "
        "bedtools sort -i bw_rpm/{wildcards.sample}.bg > bw_rpm/{wildcards.sample}.sort.bg && "
        "bedGraphToBigWig bw_rpm/{wildcards.sample}.sort.bg STAR_mm10_index/chrNameLength.txt {output} && "
        "rm bw_rpm/{wildcards.sample}.bg bw_rpm/{wildcards.sample}.sort.bg"

