nextflow.enable.dsl=2

process samtools_mpileup {
    input:
    path bam, type: 'file'
    optional path reference, type: 'file'
    optional path intervals, type: 'file'

    output:
    file "${bam.simpleName}.mpileup" into mpileup_files

    script:
    def fasta_arg = reference ? "-f $reference" : ""
    def intervals_arg = intervals ? "-l $intervals" : ""
    
    """
    samtools mpileup $fasta_arg $intervals_arg $bam > ${bam.simpleName}.mpileup
    """
}

workflow {
    // Define parameters
    reference = params.reference ? file(params.reference) : null
    intervals = params.intervals ? file(params.intervals) : null
    outdir = params.outdir // Output directory on S3

    // Read BAM file paths from a text file on S3
    bam_list = file(params.input_list) // This is the text file containing BAM paths

    // Load BAM files from the text file
    Channel
        .fromPath(bam_list)
        .set { bam_files }

    // Pass Bam file, interval file (optional), and reference genome(optional) into mpileup command
    bam_files
        .combine(reference, intervals) { bam, ref, intv -> tuple(bam, ref, intv) }
        | samtools_mpileup

    // Output files are saved to the specified outdir on S3
    mpileup_files
        .map { file -> file.copy(outdir + "/" + file.name) } // Copy files to S3 output directory
}
