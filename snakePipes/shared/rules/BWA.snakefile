### BWA ####################################################################
if paired:
    rule BWA:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "BWA/{sample}.BWA_summary.txt",
            bam = temp("BWA/{sample}.sorted.bam")# removing since we keep the sambamba output (dupmarked)
        params:
            bwa_opts = str(bwa_opts or ''),
            mate_orientation = mate_orientation,
            insert_size_max = insert_size_max
        benchmark:
            "BWA/.benchmark/BWA.{sample}.benchmark"
        threads: 24  # 1G per core
        conda: CONDA_DNA_MAPPING_ENV
        shell:
            "bwa mem "
            "-M "+bwa_index+" {input.r1} {input.r2} "
            "{params.bwa_opts} "
            "-t {threads} "
            "2> {output.align_summary} | "
            "samtools view -Sb - | "
            "samtools sort -m 2G -T ${{TMPDIR}}{wildcards.sample} -@ 2 -O bam - > {output.bam}"
else:
    rule BWA:
        input:
            fastq_dir+"/{sample}"+".fastq.gz"
        output:
            align_summary = "BWA/{sample}.BWA_summary.txt",
            bam = temp("BWA/{sample}.sorted.bam")
        params:
            bwa_opts = str(bwa_opts or '')
        benchmark:
            "BWA/.benchmark/BWA.{sample}.benchmark"
        threads: 24  # 1G per core
        conda: CONDA_DNA_MAPPING_ENV
        shell:
            "bwa mem "
            "-M "
            "-t {threads} "
            "{params.bwa_opts} "
            "-p "+bwa_index+" {input} "
            "2> {output.align_summary} | "
            "samtools view -Sbu - | "
            "samtools sort -m 2G -T ${{TMPDIR}}{wildcards.sample} -@ 2 -O bam - > {output.bam}"
