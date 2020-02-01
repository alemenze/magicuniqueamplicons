
/*
 * Create a channel for input read files
 */

Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\n" }
        .into { read_files_fastqc; read_files_trimming }


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}

 if (params.FW_primer && params.RV_primer){
/*
 * STEP 2 - Trimming
 */
 process trimming{
     tag "$name"
     publishDir "${params.outdir}/trimmed", mode:'copy'

     input:
     set val(name), file(reads) from read_files_trimming

     output:
     set val(name), file("trimmed/*.*") into ch_fastq_trimmed

     script:
     """
     mkdir -p trimmed
     cutadapt --pair-filter=any --discard-untrimmed \
     -g ${params.FW_primer} -G ${params.RV_primer} \
     -o trimmed/${reads[0]} -p trimmed/${reads[1]} \
     ${reads[0]} ${reads[1]} > cutadapt_log_${name}.txt
     """
    }
 }
 else {
     println "No Trimming performed"
 }

/*
 * STEP 3 - Panda Pairing
 */

if (params.FW_primer && params.RV_primer){
 process panda_pair_trimmed {
    publishDir "${params.outdir}/pairs", mode: 'copy'
    tag "$name"
    input:
    set val(name), file(reads) from ch_fastq_trimmed
    
    output:
    set val(name),file("*_paired.fastq") into panda_results
    
    script:
    """
    pandaseq -f ${reads[0]} -r ${reads[1]} -w ${name}_paired.fastq
    """
    }
 }
 else {
    process panda_pair_untrimmed {
    publishDir "${params.outdir}/pairs", mode: 'copy'
    tag "$name"
    input:
    set val(name), file(reads) from read_files_trimming
    
    output:
    set val(name),file("*_paired.fastq") into panda_results
    
    script:
    """
    pandaseq -f ${reads[0]} -r ${reads[1]} -w ${name}_paired.fastq
    """
    }
}

 /*
 * STEP 4 - Unique Counting
 */
  process unique_count {
     publishDir "${params.outdir}/txts", mode: 'copy'
    tag "$name"
    input:
    set val(name), file(reads) from panda_results

    output:
    set val(name), file("*.txt") into count_results

    script:
    """
    zgrep -A1 ">" ${reads} | grep -v "&--\$" | grep -v ">" | sort | uniq -c | sort -k1,1rn | tee > ${name}.txt
    """
 }
 
/*
 * STEP 5 - CSV conversion
 */
 process csv_convert {
     publishDir "${params.outdir}/output", mode: 'copy'
tag "$name"
    input:
    set val(name), file(counts) from count_results

    output:
    set val(name), file("*.csv") into csv_results

    script:
    """
    cat ${counts} | tr -s '[:blank:]' ',' > ${name}.csv
    """
 }
