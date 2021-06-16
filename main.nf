#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { build_snpeff_db } from './modules/tb-genomic-epidemiology.nf'
include { index_ref } from './modules/tb-genomic-epidemiology.nf'
include { align } from './modules/tb-genomic-epidemiology.nf'
include { call_variants } from './modules/tb-genomic-epidemiology.nf'
include { vcf_extract_subs } from './modules/tb-genomic-epidemiology.nf'
include { vcf_to_tab } from './modules/tb-genomic-epidemiology.nf'
include { annotate_variants } from './modules/tb-genomic-epidemiology.nf'
include { make_consensus } from './modules/tb-genomic-epidemiology.nf'
include { make_consensus as make_consensus_subs } from './modules/tb-genomic-epidemiology.nf'

workflow {

  ch_fastq_input = Channel.fromFilePairs( "${params.fastq_input}/*{1,2}*.fastq.gz" )

  build_snpeff_db(Channel.fromPath(params.ref).combine(Channel.fromPath(params.gff)))
  index_ref(build_snpeff_db.out.ref_data)
  align(ch_fastq_input.combine(index_ref.out))
  call_variants(align.out.combine(index_ref.out))
  vcf_extract_subs(call_variants.out.map{ it -> [it[0], it[2]] })
  annotate_variants(call_variants.out.map{ it -> [it[0], it[2]] }.combine(build_snpeff_db.out.ref_data).combine(build_snpeff_db.out.snpeff_config))
  vcf_to_tab(annotate_variants.out.combine(index_ref.out))
  make_consensus(call_variants.out.map{ it -> [it[0], it[2]] }.combine(index_ref.out).combine(Channel.of("")))
  make_consensus_subs(vcf_extract_subs.out.combine(index_ref.out).combine(Channel.of(".subs")))

  

}
