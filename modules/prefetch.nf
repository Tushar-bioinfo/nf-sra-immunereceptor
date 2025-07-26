process prefetch {

  tag "$accession"
  publishDir "work/prefetch", mode: 'copy'

  input:
  val accession
  val ngc_key

  output:
  path "${accession}.sra", emit: sra_file

  script:
  """
  mkdir -p ${accession}
  prefetch --ngc ${ngc_key} -O . -p ${accession}
  mv ${accession}/${accession}.sra .
  rm -rf ${accession}
  """
}