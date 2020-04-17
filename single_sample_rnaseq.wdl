import "star_align.wdl" as star_align
import "samtools.wdl" as samtools
import "feature_counts.wdl" as feature_counts

workflow SingleSampleRnaSeqWorkflow{

    # This workflow operates on a single "sample" and
    # assumes that all the sequence reads are contained in 
    # 1 or 2 fastq files, for single- and paired-end sequencing,
    # respectively.
    #
    # Outputs are two count files, created from primary-filtered
    # and primary-filtered + deduplicated BAM files.

    File r1_fastq
    File? r2_fastq
    File star_index_path
    File gtf
    String tag

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1_fastq, "_R1.fastq.gz") + "." + tag

    # Perform the alignment, which outputs a sorted BAM
    call star_align.perform_align as alignment{
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            gtf = gtf,
            star_index_path = star_index_path,
            sample_name = sample_name
    }

    # Index that sorted BAM
    call samtools.samtools_index as index1 {
        input:
            input_bam = alignment.sorted_bam
    }

    # Filter for primary reads only
    call samtools.samtools_primary_filter as primary_filter{
        input:
            input_bam = alignment.sorted_bam,
            input_bam_index = index1.bam_index,
            sample_name = sample_name
    }

    # Index the primary-filtered BAM
    call samtools.samtools_index as index2 {
        input:
            input_bam = primary_filter.output_bam
    }

    # Quantify the primary-filtered BAM
    call feature_counts.count_reads as quantify_primary {
        input:
            input_bam = primary_filter.output_bam,
            gtf = gtf,
            sample_name = sample_name,
            tag = "primary"
    }

    output {
        File unfiltered_bam = alignment.sorted_bam
        File unfiltered_bam_index = index1.bam_index
        File primary_bam = primary_filter.output_bam
        File primary_bam_index = index2.bam_index 
        File primary_filter_feature_counts_file = quantify_primary.count_output
        File primary_filter_feature_counts_summary = quantify_primary.count_output_summary
        File star_log = alignment.final_log
    }

}
