
VNTR_INPUT="/oak/stanford/groups/euan/projects/tannerj/ALNG_DB/VNTR.vamos.cohort_combined.ALNG_DB.AF_annotated.sorted.txt"
k = [10]
chroms = ["chr%d" % i for i in range(1,23)] + ['chrX', 'chrY']
rule all:
    input:
      expand("VNTR.mean_neighbor_distance.kneighbors_{k}.extreme_outliers.csv", k=k)

rule VNTR_outlier:
    threads: 16
    resources:
        mem=128,
        time=12
    input:
      VNTR_INPUT
    params: 
        script = "call_VNTR_outliers.R",
        tmpdir = "/tmp/tannerj/VNTR_outliers.chrom_{chr}.k_{k}"
    output:
        temp("VNTR.mean_neighbor_distance.kneighbors_{k}.extreme_outliers.chrom_{chr}.csv")
    conda: 'r'
    shell: """
        mkdir -p {params.tmpdir}
        
        head -1 {input} > {params.tmpdir}/input_VNTR.tsv
        awk "\$3==\\"{wildcards.chr}\\"" {input} >> {params.tmpdir}/input_VNTR.tsv
        Rscript {params.script} {params.tmpdir}/input_VNTR.tsv {output} {wildcards.k}
        rm -rf {params.tmpdir}
    """

rule concat_results:
  threads: 1 
  resources:
    mem=24,
    time=4
  input:
    expand("VNTR.mean_neighbor_distance.kneighbors_{{k}}.extreme_outliers.chrom_{chr}.csv", chr=chroms),
  output:
    "VNTR.mean_neighbor_distance.kneighbors_{k}.extreme_outliers.csv"
  shell: """
    head -1 {input[0]} > {output}
    cat {input} | grep -E -v "^VNTR_ID" >> {output}
  """

#rule convert_to_vcf:
#  ""
#
