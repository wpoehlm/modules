nextflow.enable.dsl=2

process samtools_mpileup {
    input:
    path bam, type: 'file'
    path reference, type: 'file'

    output:
    file "${bam.simpleName}.mpileup" into mpileup_files

    script:
    """
    samtools mpileup -f $reference $bam > ${bam.simpleName}.mpileup
    """
}

workflow {
    // Define parameters
    reference = file(params.reference)
    outdir = params.outdir // Output directory on S3

    // Read BAM file paths from a text file on S3
    bam_list = file(params.input_list) // This is the text file containing BAM paths

    // Load BAM files from the text file
    Channel
        .fromPath(bam_list)
        .set { bam_files }

    // Pass Bam file and reference genome into mpileup command
    bam_files
        .combine(reference) { bam, ref -> tuple(bam, ref) }
        | samtools_mpileup

    // Output files are saved to the specified outdir on S3
    mpileup_files
        .map { file -> file.copy(outdir + "/" + file.name) } // Copy files to S3 output directory
}
