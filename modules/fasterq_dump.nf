process fasterq_dump {

  tag "$accession"
  publishDir "work/fastq", mode: 'copy'

  input:
  path sra_file
  val accession
  val ngc_key

  output:
  tuple val(accession), path("${accession}_1.fastq"), path("${accession}_2.fastq")

  script:
  """
  mkdir -p scratch output
  timeout 4h fasterq-dump ${sra_file} -t scratch -e 1 -x --progress -O output --ngc ${ngc_key}
  mv output/${accession}_1.fastq . 
  mv output/${accession}_2.fastq . 
  rm -rf scratch output
  """
}