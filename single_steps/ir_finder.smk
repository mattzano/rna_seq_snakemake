import os
# a top level folder where the bams reside
project_dir = "/SAN/vyplab/alb_projects/data/bdnf_4su_i3lmn/trial_run/"
out_spot = "IRFinder/"
bam_spot = "STAR_aligned/"
bam_suffix = ".Aligned.sorted.out.bam"
IRfinder_path = "/SAN/vyplab/alb_projects/tools/IRFinder/bin/IRFinder"
IRfinder_reference = "/SAN/vyplab/alb_projects/tools/IRFinder/REF/human_gencode_v34/irfinder/"
#IRfinder_reference = "/SAN/vyplab/alb_projects/tools/IRFinder/REF/mouse_gencode_v25/irfinder/"
# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------

output_dir = os.path.join(project_dir,out_spot)
bam_dir = os.path.join(project_dir,bam_spot)

SAMPLES, = glob_wildcards(bam_dir + "{sample}" + bam_suffix)
print(SAMPLES)

rule all:
  input:
    expand(output_dir + "{sample}_namesorted.bam", sample = SAMPLES),
    expand(output_dir + "{sample}/IRFinder-IR-nondir.txt", sample = SAMPLES)

rule name_sort:
    input:
        aligned_bam = bam_dir + "{sample}" + bam_suffix
    output:
       out_name = temp(output_dir + "{sample}_namesorted.bam")
    shell:
        """
        mkdir -p {output_dir}
        t=/scratch0/$USER/$RANDOM
        mkdir -p $t
        samtools sort -T $t -n -@ 2 {input.aligned_bam} -o {output.out_name}
        """
rule run_ir_finder:
    input:
        output_dir + "{sample}_namesorted.bam"
    wildcard_constraints:
        sample="|".join(SAMPLES)
    output:
        output_dir + "{sample}/IRFinder-IR-nondir.txt"
    params:
        output_folder = output_dir + "{sample}"
    shell:
        """
        rm -r {params.output_folder}
        {IRfinder_path} -m BAM -r {IRfinder_reference} -d {params.output_folder} {input}
        """
