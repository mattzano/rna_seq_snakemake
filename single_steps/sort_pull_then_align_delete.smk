import os
# a top level folder where the bams reside
my_new_salmon = "/SAN/vyplab/alb_projects/tools/salmon-1.5.1_linux_x86_64/bin/salmon"
project_dir = "/SAN/vyplab/alb_projects/data/klim_tdp/"
out_spot = "name_sortbams/"
bam_spot = "STAR_aligned/"
fastq_dir = "pulled_fastq/"
bam_suffix = ".bam"
end_type = "pe"
out_spot = "salmon_full/"
salmon_strand_info = "A"
salmon_index = "/SAN/vyplab/vyplab_reference_genomes/salmon/human/GRCh38/full/gencode_v34.kmer_31/"
gtf = "/SAN/vyplab/vyplab_reference_genomes/salmon/human/gencode_v34geneMap.txt"
# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------

output_dir = os.path.join(project_dir,out_spot)
bam_dir = os.path.join(project_dir,bam_spot)
fastq_dir = os.path.join(project_dir, fastq_dir)

SAMPLES, = glob_wildcards(bam_dir + "{sample}" + bam_suffix)
print(SAMPLES)

rule all:
  input:
    expand(output_dir + "{sample}/" + "quant.sf", sample = SAMPLES)

rule name_sort:
    input:
        aligned_bam = bam_dir + "{sample}" + bam_suffix
    output:
       temp(output_dir + "{sample}_namesorted.bam")
    shell:
        """
        mkdir -p {output_dir}
        samtools sort -n -@ 2 {input.aligned_bam} -o {output}
        """
if end_type == "pe":
  rule bam_to_fastq:
      input:
          name_sort_bam = output_dir + "{sample}_namesorted.bam"
      output:
          one = temp(fastq_dir + "{sample}_1.merged.fastq"),
          two = temp(fastq_dir + "{sample}_2.merged.fastq")
      shell:
          """
          bedtools bamtofastq -i {input} \
                        -fq {output.one} \
                        -fq2 {output.two}
          """
  rule gunzip_fastq:
      input:
          one = fastq_dir + "{sample}_1.merged.fastq",
          two = fastq_dir + "{sample}_2.merged.fastq"
      output:
          one_out = temp(fastq_dir + "{sample}_1.merged.fastq.gz"),
          two_out = temp(fastq_dir + "{sample}_2.merged.fastq.gz")
      shell:
          """
          gzip {input.one}
          gzip {input.two}
          """
else:
  rule bam_to_fastq:
      input:
          name_sort_bam = output_dir + "{sample}_namesorted.bam"
      output:
          one = temp(fastq_dir + "{sample}_1.merged.fastq")
      shell:
          """
          bedtools bamtofastq -i {input} \
                        -fq {output.one}
          """
  rule gunzip_fastq:
      input:
          one = fastq_dir + "{sample}_1.merged.fastq"
      output:
          one_out = temp(fastq_dir + "{sample}_1.merged.fastq.gz")
      shell:
          """
          gzip {input.one}
          """

rule salmon_quant:
    input:
        fast1 = fastq_dir  + "{sample}_1.merged.fastq.gz",
        fast2 = fastq_dir  + "{sample}_2.merged.fastq.gz",
    output:
        output_dir + "{sample}/" + "quant.sf"
    params:
        out = output_dir + "{sample}",
        index_dir = salmon_index,
        threads = 4
    shell:
        """
        {my_new_salmon} quant \
        --index {params.index_dir} \
        --libType {salmon_strand_info} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --geneMap {gtf} \
        --threads {params.threads} \
        --gcBias \
        --seqBias \
        --numBootstraps 50 \
        -o {params.out} 
        """
