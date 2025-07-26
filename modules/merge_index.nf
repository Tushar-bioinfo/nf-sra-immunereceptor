process merge_index_bam {

  tag "$accession"
  publishDir "results/final_bams", mode: 'copy'

  input:
  tuple val(accession), path(bam_parts)

  output:
  tuple val(accession), path("${accession}.bam"), path("${accession}.bam.bai")

  script:
  """
  sambamba merge ${accession}.bam ${bam_parts.join(' ')}
  sambamba index ${accession}.bam
  """
}
