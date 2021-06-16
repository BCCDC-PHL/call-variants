process build_snpeff_db {

  tag { ref_id }

  input:
  tuple path(ref_fa), path(ref_gff)

  output:
  tuple val(ref_id), path("data"), emit: ref_data
  path("snpEff.config"), emit: snpeff_config

  script:
  ref_header = ref_fa.toRealPath().toFile().newReader().readLine()
  ref_id = ref_header.replaceAll(">", "").split()[0]
  """
  mkdir -p data/ref
  cp ${ref_fa} data/ref/sequences.fa
  cp ${ref_gff} data/ref/genes.gff
  echo '##FASTA' >> data/ref/genes.gff
  cat data/ref/sequences.fa >> data/ref/genes.gff
  echo "codon.Bacterial_and_Plant_Plastid : TTT/F, TTC/F, TTA/L, TTG/L+, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V+, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G" >> snpEff.config
  echo "ref.genome : Reference Genome' >> snpEff.config
  echo "\tref.chromosomes : ${ref_id}' >> snpEff.config
  echo "\tref.${ref_id}.codonTable : Bacterial_and_Plant_Plastid' >> snpEff.config
  
  snpEff build -c snpEff.config -gff3 ref
  """
}

process index_ref {

  tag { ref_id }

  input:
  tuple val(ref_id), path(ref_data)

  output:
  tuple val(ref_id), path("data")

  script:
  """
  mv ${ref_data} tmp_data
  cp -R tmp_data data
  cd data
  bwa index ref/sequences.fa
  samtools faidx ref/sequences.fa
  """
  
}

process align {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}.bam*"

  input:
  tuple val(grouping_key), path(reads), val(ref_id), path(ref_data)

  output:
  tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

  script:
  if (grouping_key =~ '_S[0-9]+_') {
    sample_id = grouping_key.split("_S[0-9]+_")[0]
  } else if (grouping_key =~ '_') {
    sample_id = grouping_key.split("_")[0]
  } else {
    sample_id = grouping_key
  }
  rg_string = "@RG\\tID:${sample_id}\\tSM:${sample_id}"
  """
  bwa mem -t ${task.cpus} -Y -M -R "${rg_string}" ${ref_data}/ref/sequences.fa ${reads} | \
  samclip --max ${params.maxsoft} --ref ${ref_data}/ref/sequences.fa.fai | \
  samtools sort -@ ${task.cpus} -n -l 0 | \
  samtools fixmate -@ ${task.cpus} -m - - | \
  samtools sort -@ ${task.cpus} -l 0 | \
  samtools markdup -@ ${task.cpus} -r - - \
  > ${sample_id}.bam
  samtools index ${sample_id}.bam
  """
}

process call_variants {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}*.vcf"

  input:
  tuple val(sample_id), path(alignment), path(alignment_index), val(ref_id), path(ref_data)

  output:
  tuple val(sample_id), path("${sample_id}.raw.vcf"), path("${sample_id}.filt.vcf")

  script:
  bcf_filter = "FMT/GT=\"1/1\" && QUAL>=${params.minqual} && FMT/DP>=${params.mincov} && (FMT/AO)/(FMT/DP)>=${params.minfrac}"
  keep_vcf_tags = "^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL"
  """
  freebayes \
    --ploidy 2 \
    --pvar 0.0 \
    --min-alternate-count 2 \
    --min-alternate-fraction 0.05 \
    --min-repeat-entropy 1.0 \
    --min-base-quality ${params.basequal} \
    --min-mapping-quality ${params.mapqual} \
    --strict-vcf \
    -f ${ref_data}/ref/sequences.fa \
    ${alignment} > ${sample_id}.raw.vcf
  bcftools view --include '${bcf_filter}' ${sample_id}.raw.vcf | \
  vt normalize -r ${ref_data}/ref/sequences.fa - | \
  bcftools annotate --remove ${keep_vcf_tags} > ${sample_id}.filt.vcf
  """
}

process make_consensus {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}.consensus${consensus_filename_suffix}.fa"

  input:
  tuple val(sample_id), path(vcf), val(ref_id), path(ref_data), val(consensus_filename_suffix)

  output:
  tuple val(sample_id), path("${sample_id}.consensus${consensus_filename_suffix}.fa")

  script:
  """
  bcftools convert -Oz -o ${sample_id}.vcf.gz ${vcf}
  bcftools index -f ${sample_id}.vcf.gz
  bcftools consensus --sample ${sample_id} -f ${ref_data}/ref/sequences.fa -o consensus.fa ${sample_id}.vcf.gz
  sed 's/>.*/>${sample_id}/' consensus.fa > ${sample_id}.consensus${consensus_filename_suffix}.fa
  """
}

process annotate_variants {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}.vcf"

  input:
  tuple val(sample_id), path(filtered_vcf), val(ref_id), path(ref_data), path(snpeff_config)

  output:
  tuple val(sample_id), path("${sample_id}.vcf")

  script:
  """
  snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr \
    -c ${snpeff_config} -dataDir ${ref_data} ref ${filtered_vcf} > ${sample_id}.vcf
  """
}

process vcf_extract_subs {
  tag { sample_id }

  input:
  tuple val(sample_id), path(input_vcf)

  output:
  tuple val(sample_id), path("${sample_id}.subs.vcf")
  
  script:
  """
  vcf_extract_subs ${input_vcf} > ${sample_id}.subs.vcf
  """
}

process vcf_to_tab {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}.snps.tsv"

  input:
  tuple val(sample_id), path(input_vcf), val(ref_id), path(ref_data)

  output:
  tuple val(sample_id), path("${sample_id}.snps.tsv")
  
  script:
  """
  vcf_to_tab --ref ${ref_data}/ref/sequences.fa --gff ${ref_data}/ref/genes.gff --vcf ${input_vcf} > ${sample_id}.snps.tsv
  """
}

