process slice_receptors {

  tag "$accession"
  publishDir "work/sliced_bams", mode: 'copy'

  input:
  tuple val(accession), path(sorted_bam), path(bam_index)

  output:
  tuple val(accession), path("${accession}_merged.bam"), path("${accession}_merged.bam.bai")

  script:
  """
  # Extract immune receptor regions
  samtools view -b ${sorted_bam} chr14:105586437-106879844 -o ${accession}_IGH.bam
  samtools view -b ${sorted_bam} chr2:88857361-90235368   -o ${accession}_IGK.bam
  samtools view -b ${sorted_bam} chr22:22026076-22922913  -o ${accession}_IGL.bam
  samtools view -b ${sorted_bam} chr14:21621904-22552132  -o ${accession}_TRA.bam
  samtools view -b ${sorted_bam} chr7:142299011-142813287 -o ${accession}_TRB.bam
  samtools view -b ${sorted_bam} chr7:38240024-38368055   -o ${accession}_TRG.bam
  samtools view -b -f 4 ${sorted_bam}                     -o ${accession}_unmapped.bam

  # Merge all extracted BAMs
  sambamba merge ${accession}_merged.bam \
    ${accession}_IGH.bam \
    ${accession}_IGK.bam \
    ${accession}_IGL.bam \
    ${accession}_TRA.bam \
    ${accession}_TRB.bam \
    ${accession}_TRG.bam \
    ${accession}_unmapped.bam

  # Index final merged BAM
  sambamba index ${accession}_merged.bam

  # Clean intermediate receptor BAMs
  rm -f ${accession}_IGH.bam ${accession}_IGK.bam ${accession}_IGL.bam \\
        ${accession}_TRA.bam ${accession}_TRB.bam ${accession}_TRG.bam \\
        ${accession}_unmapped.bam
  """
}
