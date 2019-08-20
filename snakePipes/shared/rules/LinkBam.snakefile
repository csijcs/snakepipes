rule link_bam:
    input:
        indir+"/{sample}"+bam_ext
    output:
        mapping_prg+"/{sample}.bam"
    shell:
        "( [ -f {output} ] || ( ln -s -r {input} {output} ) )"

rule sorting_bam:
    input:
        mapping_prg+"/{sample}.bam"
    log:
        mapping_prg+"/logs/{sample}.log"
    output:
        temp(mapping_prg+"/{sample}.sorted.bam")
    conda: CONDA_SHARED_ENV
    shell:
        """
            samtools sort -@ 4 -m 2G {input} -o {output} 2> {log}
        """

rule picard_mark_duplicates:
    input:
        mapping_prg+"/{sample}.sorted.bam"
    output:
        bam=temp(mapping_prg+"/{sample}.marked.bam"), # duplicate marked
        metrics=mapping_prg + "/logs/{sample}.dup.metrics"
    log:
        out=mapping_prg + "/logs/{sample}.picard_markdup.out",
        err=mapping_prg + "/logs/{sample}.picard_markdup.err"
    benchmark: mapping_prg + "/.benchmark/picard_markdup.{sample}.benchmark"
    conda: CONDA_CHIPSEQ_ENV
    threads: 10
    shell:
        """
            picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} VALIDATION_STRINGENCY=LENIENT -XX:ParallelGCThreads={threads} 2> {log.err} > {log.out}
           """

rule samtools_filter:
    input:
        bam = mapping_prg+"/{sample}.marked.bam"
    output:
        bam_out = "filtered_bam/{sample}.filtered.bam"
    conda: CONDA_SHARED_ENV
    shell: "samtools view -bq 20 {input} > {output}"

#rule link_bam_bai_external:
#    input:
#        bam = mapping_prg+"/{sample}.bam",
#        bai = mapping_prg+"/{sample}.bam.bai"
#    output:
#        bam_out = "filtered_bam/{sample}.filtered.bam",
#        bai_out = "filtered_bam/{sample}.filtered.bam.bai",
#    shell:
#        "( [ -f {output.bam_out} ] || ( ln -s -r {input.bam} {output.bam_out} && ln -s -r {input.bai} {output.bai_out} ) )"

rule samtools_index_filtered:
     input:
         "filtered_bam/{sample}.filtered.bam"
     output:
         "filtered_bam/{sample}.filtered.bam.bai"
     conda: CONDA_SHARED_ENV
     shell: "samtools index {input}"

rule sambamba_flagstat:
       input:
           mapping_prg+"/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       conda: CONDA_SAMBAMBA_ENV
       shell: """
           sambamba flagstat -p {input} > {output}
           """
