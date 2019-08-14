rule phantom:
    input:
        test_bam='filtered_bam/{sample}.filtered.bam',
        control_bam=lambda wildcards: "filtered_bam/"+get_control(wildcards.sample)+".filtered.bam" if get_control(wildcards.sample)
                else []
    output:
        'phantom/{sample}.phantom'
    log:
        'phantom/logs/{sample}.phantom.log'
    conda: CONDA_PHANTOM_ENV
    params:
        script=os.path.join(maindir, "shared", "rscripts", "run_spp.R")
    shell:
        """
            Rscript {params.script} -rf -c={input.test_bam} -i={input.control_bam} -savp -out={output} &> {log}
        """

#def get_ext(sample):
#        phantom = "phantom/{sample}.phantom"
#        command = "awk '{print $3 }' < " + phantom + """ | tr ","  "\t" | awk '{if($1!=0) print $1; else print $2}' """
#        p = sp.Popen([command], stdout=sp.PIPE, shell=True)
#        return chr(float(p.stdout.read().decode('utf-8').strip('\n'))/2)

rule fragment_length:
    input: "phantom/{sample}.phantom"
    output: temp("phantom/{sample}.fragment_length")
    shell: """
            awk '{{print $3}}' < {input} | tr ',' '\\t' | awk '{{if($1!=0) print $1; else print $2}}' > {output}
            """
