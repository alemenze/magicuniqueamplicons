
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
     publishDir "${params.outdir}/trimmed", mode:'copy',
     saveAs 

     input set val(pair_id), file(reads) from read_files_trimming

     output:
     set val(name), file "trimmed/*.*" into (ch_fastq_trimmed)

     script:
     """
     mkdir -p trimmed
     cutadapt --pair-filter=any --discard-untrimmed -g ${params.FW_primer} -G ${params.RV_primer} -o trimmed/$folder${params.split}${reads[0]} -p trimmed/$folder${params.split}${reads[1]} ${reads[0]} ${reads[2]} > cutadapt_log_${pair_id}.txt
     """
 }
 }
 else {
     println "No Trimming performed"
 }

/*
 * STEP 3 - Panda Pairing
 */
 process panda_pair {
    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    
    if (params.FW_primer && params.RV_primer){
        set val(name), file(reads) from ch_fastq_trimmed
    }
    else {
        set val(name), file(reads) from read_files_trimming
    }

    output:
    set val(name),file("*_paired.fastq") into panda_results
    

    script:
    """
    pandaseq -f ${reads[0]} -r ${reads[1]} -w ${name}_paired.fastq
    """

 }

 /*
 * STEP 4 - Unique Counting
 */
  process unique_count {
     publishDir "${params.outdir}/txts", mode: 'copy'

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
 * STEP 4 - CSV conversion
 */
 process csv_convert {
     publishDir "${params.outdir}/output", mode: 'copy'

    input:
    set val(name), file(counts) from count_results

    output:
    set val(name), file("*.csv") into csv_results

    script:
    """
    cat ${counts} | tr -s '[:blank:]' ',' > ${name}.csv
    """
 }
