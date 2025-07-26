process sort_bam {

  tag "$accession"
  publishDir "work/sorted_bam", mode: 'copy'

  input:
  tuple val(accession), path(unsorted_bam), path(star_log)

  output:
  tuple val(accession), path("${accession}.sorted.bam"), path("${accession}.sorted.bam.bai")

  script:
  """
  sambamba sort \
    --memory-limit 100GB \
    --nthreads 6 \
    --tmpdir ./tmp_sort_${accession} \
    --out ${accession}.sorted.bam \
    ${unsorted_bam}

  samtools index ${accession}.sorted.bam
  """
}