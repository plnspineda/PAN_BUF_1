nextflow.enable.dsl=2

// --- Parameters ---
params.samples    = "samples.tsv"
params.reference  = "/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/reference_fa/renamed/ARVUO_withX_renamed.fa"
params.chroms     = (1..24).collect { it.toString() } + "X"
params.outdir     = "results"

// --- Workflow ---
workflow {
    // load samples
    samples_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sampleID, file(row.read1), file(row.read2)) }
    
    def ref_p = file(params.reference)
    ref_files = Channel.fromPath("${ref_p.getParent()}/${ref_p.baseName}*").collect()

    // QC and alignment
    SEQKIT_STATS(samples_ch)
    BWA_MEM(samples_ch, ref_files)
    MARK_DUPLICATES(BWA_MEM.out.bam)

    // variant calling
    hc_input_ch = MARK_DUPLICATES.out.bam_bai
        .combine(Channel.fromList(params.chroms))

    HAPLOTYPE_CALLER(hc_input_ch, ref_files)

    // genotyping
    db_import_ch = HAPLOTYPE_CALLER.out.gvcf
        .map { gvcf, chrom -> tuple(chrom, gvcf) }
        .groupTuple()

    GENOMICSDB_IMPORT(db_import_ch)
    GENOTYPE_GVCFS(GENOMICSDB_IMPORT.out.db, ref_files)

    // filtering
    SPLIT_VARIANTS(GENOTYPE_GVCFS.out.vcf)

    // FIX: FILTER_SNPS now passes a tuple of 3 elements to FILTER_SNPS_DP10
    FILTER_SNPS(SPLIT_VARIANTS.out.snps)
    FILTER_SNPS_DP10(FILTER_SNPS.out)
    
    FILTER_INDELS(SPLIT_VARIANTS.out.indels)
}

// --- Processes ---

process BWA_MEM {
    tag "$sampleID"
    cpus 8
    memory '64 GB'
    publishDir "${params.outdir}/bwa_mem", mode: 'copy'

    input:
    tuple val(sampleID), path(r1), path(r2)
    path ref

    output:
    tuple val(sampleID), path("${sampleID}.sorted.bam"), emit: bam

    script:
    def fasta = ref.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') }
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sampleID}\\tSM:${sampleID}" ${fasta} ${r1} ${r2} | \
    samtools sort -@ ${task.cpus} -o ${sampleID}.sorted.bam -
    """
}

process MARK_DUPLICATES {
    tag "$sampleID"
    publishDir "${params.outdir}/mark_dups", mode: 'copy'
    
    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("${sampleID}.md.bam"), path("${sampleID}.md*.bai"), emit: bam_bai

    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sampleID}.md.bam \
        -M ${sampleID}.metrics.txt \
        -CREATE_INDEX true
    """
}

process HAPLOTYPE_CALLER {
    tag "$sampleID / $chrom"
    publishDir "${params.outdir}/hap_call", mode: 'copy'

    input:
    tuple val(sampleID), path(bam), path(bai), val(chrom)
    path ref

    output:
    tuple path("${sampleID}.${chrom}.g.vcf.gz"), val(chrom), emit: gvcf

    script:
    def fasta = ref.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') }
    """
    gatk HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -L ${chrom} \
        -O ${sampleID}.${chrom}.g.vcf.gz \
        -ERC GVCF
    """
}

process GENOMICSDB_IMPORT {
    tag "$chrom"
    cpus 4
    memory '32 GB'
    
    input:
    tuple val(chrom), path(gvcfs)

    output:
    tuple val(chrom), path("genomicsdb_${chrom}"), emit: db

    script:
    """
    for f in ${gvcfs}; do
        sample=\$(basename \$f | cut -d. -f1)
        abs_path=\$(readlink -f \$f)
        echo -e "\${sample}\\t\${abs_path}" >> sample_map.tsv
    done

    gatk --java-options "-Xmx28g" GenomicsDBImport \
        --genomicsdb-workspace-path genomicsdb_${chrom} \
        --sample-name-map sample_map.tsv \
        --intervals ${chrom} \
        --overwrite-existing-genomicsdb-workspace true
    """
}

process GENOTYPE_GVCFS {
    tag "$chrom"
    cpus 8
    memory '64 GB'
    publishDir "${params.outdir}/vcf_raw", mode: 'copy'

    input:
    tuple val(chrom), path(genomicsdb)
    path ref

    output:
    tuple val(chrom), path("combined.${chrom}.vcf.gz"), path("combined.${chrom}.vcf.gz.tbi"), emit: vcf

    script:
    def fasta = ref.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') }
    """
    DB_PATH=\$(readlink -f ${genomicsdb})
    
    gatk --java-options "-Xmx56g" GenotypeGVCFs \
        -R ${fasta} \
        -V "gendb://\${DB_PATH}" \
        -O combined.${chrom}.vcf.gz \
        -L ${chrom}

    gatk IndexFeatureFile -I combined.${chrom}.vcf.gz
    """
}

process SPLIT_VARIANTS {
    tag "$chrom"

    input:
    tuple val(chrom), path(vcf), path(index)

    output:
    tuple val(chrom), path("raw_snps.${chrom}.vcf.gz"), path("raw_snps.${chrom}.vcf.gz.tbi"),   emit: snps
    tuple val(chrom), path("raw_indels.${chrom}.vcf.gz"), path("raw_indels.${chrom}.vcf.gz.tbi"), emit: indels

    script:
    """
    gatk SelectVariants -V ${vcf} -select-type SNP -O raw_snps.${chrom}.vcf.gz
    gatk IndexFeatureFile -I raw_snps.${chrom}.vcf.gz

    gatk SelectVariants -V ${vcf} -select-type INDEL -O raw_indels.${chrom}.vcf.gz
    gatk IndexFeatureFile -I raw_indels.${chrom}.vcf.gz
    """
}

process FILTER_SNPS {
    tag "$chrom"
    publishDir "${params.outdir}/vcf_filtered", mode: 'copy'

    input:
    tuple val(chrom), path(vcf), path(index)

    output:
    // Now returns a tuple to match the input expectation of FILTER_SNPS_DP10
    tuple val(chrom), path("filtered_snps.${chrom}.vcf.gz"), path("filtered_snps.${chrom}.vcf.gz.tbi")

    script:
    """
    gatk VariantFiltration \\
        -V ${vcf} \\
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
        --filter-name "filter1" \\
        --cluster-size 3 \\
        --cluster-window-size 10 \\
        -O filtered_snps.${chrom}.vcf.gz

    gatk IndexFeatureFile -I filtered_snps.${chrom}.vcf.gz
    """
}

process FILTER_SNPS_DP10 {
    tag "$chrom"
    publishDir "${params.outdir}/vcf_depth_filtered", mode: 'copy'

    input:
    tuple val(chrom), path(vcf), path(index)

    output:
    path "depth_filtered_snps.${chrom}.vcf.gz"

    script:
    """
    gatk VariantFiltration \\
        -V ${vcf} \\
        --filter-expression "DP < 10 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
        --filter-name "filter1" \\
        --genotype-filter-expression "GQ < 20.0" \\
        --genotype-filter-name "lowGQ" \\
        --cluster-size 3 \\
        --cluster-window-size 10 \\
        -O depth_filtered_snps.${chrom}.vcf.gz
    """
}

process FILTER_INDELS {
    tag "$chrom"
    publishDir "${params.outdir}/vcf_filtered", mode: 'copy'

    input:
    tuple val(chrom), path(vcf), path(index)

    output:
    path "filtered_indels.${chrom}.vcf.gz"

    script:
    """
    gatk VariantFiltration \\
        -V ${vcf} \\
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || ReadPosRankSum < -8.0" \\
        --filter-name "filter1" \\
        --genotype-filter-expression "GQ < 20.0" \\
        --genotype-filter-name "lowGQ" \\
        --cluster-size 3 \\
        --cluster-window-size 10 \\
        -O filtered_indels.${chrom}.vcf.gz
    """
}

process SEQKIT_STATS {
    tag "$sampleID"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sampleID), path(r1), path(r2)

    output:
    path "${sampleID}_stats.tsv"

    script:
    """
    seqkit stats ${r1} ${r2} -T -a > ${sampleID}_stats.tsv
    """
}
