// Import modules
include { prefetch }        from './modules/prefetch.nf'
include { fasterq_dump }    from './modules/fasterq_dump.nf'
include { star_align }      from './modules/star_align.nf'
include { sort_bam }        from './modules/sort_bam.nf'
include { slice_receptors } from './modules/slice_receptors.nf'

workflow {

  // Read manifest and get accessions
  Channel
    .fromPath(params.manifest)
    .splitCsv(header:true)
    .map { row -> row['Run'] }
    .set { accessions }

  // processing
  sra_files    = prefetch(accessions, params.ngc_key)
  fastq_files  = fasterq_dump(sra_files, params.ngc_key)
  aligned_bam  = star_align(fastq_files, params.ref_dir)
  sorted_bam   = sort_bam(aligned_bam)
  final_bams   = slice_receptors(sorted_bam)
}
