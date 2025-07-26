process star_align {

  tag "$accession"
  publishDir "work/star", mode: 'copy'

  input:
  tuple val(accession), path(read1), path(read2)
  val ref_dir
  val threads

  output:
  tuple val(accession), path("${accession}_Aligned.out.bam"), path("${accession}_Log.final.out")

  script:
  """
  STAR --runThreadN ${threads} \
       --genomeDir ${ref_dir} \
       --readFilesIn ${read1} ${read2} \
       --outFileNamePrefix ${accession}_ \
       --outSAMtype BAM Unsorted \
       --outFilterScoreMinOverLread 0.2 \
       --outFilterMatchNminOverLread 0.2 \
       --outFilterMismatchNoverReadLmax 0.04 \
       --outSAMunmapped Within \
       --outSAMattributes NH HI AS nM MD
  """
}
