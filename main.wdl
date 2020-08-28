import "single_sample_rnaseq.wdl" as single_sample_rnaseq
import "feature_counts.wdl" as feature_counts
import "multiqc.wdl" as multiqc
import "fastqc.wdl" as fastqc
import "disambiguate.wdl" as disambiguate

workflow StarDisambiguateWorkflow{
    # This workflow is a 'super' workflow that parallelizes
    # the analysis over multiple samples

    Array[File] r1_files
    Array[File] r2_files
    Int trim_length_bp
    File human_star_index_path
    File human_gtf
    File mouse_star_index_path
    File mouse_gtf
    String git_repo_url
    String git_commit_hash
    String human_tag = "human"
    String mouse_tag = "mouse"

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)


    scatter(item in fastq_pairs){

        call trim_reads {
            input:
            r1 = item.left,
            r2 = item.right,
            trim_length_bp = trim_length_bp
        }

        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = trim_reads.trimmed_r1
        }

        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = trim_reads.trimmed_r2
        }

        call single_sample_rnaseq.SingleSampleRnaSeqWorkflow as human_single_sample_process{
            input:
                r1_fastq = trim_reads.trimmed_r1,
                r2_fastq = trim_reads.trimmed_r2,
                star_index_path = human_star_index_path,
                gtf = human_gtf,
                tag = human_tag
        }

        call single_sample_rnaseq.SingleSampleRnaSeqWorkflow as mouse_single_sample_process{
            input:
                r1_fastq = trim_reads.trimmed_r1,
                r2_fastq = trim_reads.trimmed_r2,
                star_index_path = mouse_star_index_path,
                gtf = mouse_gtf,
                tag = mouse_tag
        }

        call disambiguate.DisambiguateBam as bam_disambiguate {
            input:
                human_bam = human_single_sample_process.primary_bam,
                mouse_bam = mouse_single_sample_process.primary_bam,
                human_bam_idx = human_single_sample_process.primary_bam_index,
                mouse_bam_idx = mouse_single_sample_process.primary_bam_index,
                human_gtf = human_gtf,
                mouse_gtf = mouse_gtf,
                human_tag = human_tag,
                mouse_tag = mouse_tag
        }
    }

    # merge the counts from the "full" human BAM
    call feature_counts.concatenate as merge_human_primary_counts {
        input:
            count_files = human_single_sample_process.primary_filter_feature_counts_file,
            output_filename = "raw_primary_counts.human.tsv"
    }

    # merge the counts from the "full" mouse BAM
    call feature_counts.concatenate as merge_mouse_primary_counts {
        input:
            count_files = mouse_single_sample_process.primary_filter_feature_counts_file,
            output_filename = "raw_primary_counts.mouse.tsv"
    }

    # now merge the counts from the ambig/disambig counts:
    call feature_counts.concatenate as merge_human_ambig_counts {
        input:
            count_files = bam_disambiguate.human_ambig_counts,
            output_filename = "ambig_counts.human.tsv"
    }

    call feature_counts.concatenate as merge_mouse_ambig_counts {
        input:
            count_files = bam_disambiguate.mouse_ambig_counts,
            output_filename = "ambig_counts.mouse.tsv"
    }

    call feature_counts.concatenate as merge_human_disambig_counts {
        input:
            count_files = bam_disambiguate.human_disambig_counts,
            output_filename = "disambig_counts.human.tsv"
    }

    call feature_counts.concatenate as merge_mouse_disambig_counts {
        input:
            count_files = bam_disambiguate.mouse_disambig_counts,
            output_filename = "disambig_counts.mouse.tsv"
    }

    output {
        File raw_human_counts = merge_human_primary_counts.count_matrix
        File raw_mouse_counts = merge_mouse_primary_counts.count_matrix

        File human_ambig_counts = merge_human_ambig_counts.count_matrix
        File mouse_ambig_counts = merge_mouse_ambig_counts.count_matrix
        File human_disambig_counts = merge_human_disambig_counts.count_matrix
        File mouse_disambig_counts = merge_mouse_disambig_counts.count_matrix

        Array[File] human_primary_bam = human_single_sample_process.primary_bam
        Array[File] human_primary_bam_idx = human_single_sample_process.primary_bam_index
        Array[File] mouse_primary_bam = mouse_single_sample_process.primary_bam
        Array[File] mouse_primary_bam_idx = mouse_single_sample_process.primary_bam_index

        Array[File] human_ambig_bam = bam_disambiguate.human_ambig_bam
        Array[File] human_ambig_bam_idx = bam_disambiguate.human_ambig_bam_index
        Array[File] mouse_ambig_bam = bam_disambiguate.mouse_ambig_bam
        Array[File] mouse_ambig_bam_idx = bam_disambiguate.mouse_ambig_bam_index
        Array[File] human_disambig_bam = bam_disambiguate.human_disambig_bam
        Array[File] human_disambig_bam_idx = bam_disambiguate.human_disambig_bam_index
        Array[File] mouse_disambig_bam = bam_disambiguate.mouse_disambig_bam
        Array[File] mouse_disambig_bam_idx = bam_disambiguate.mouse_disambig_bam_index
    }

    meta {
        workflow_title : "STAR + Disambiguate"
        workflow_short_description : "Disambiguating PDx samples using STAR and AZ Disambiguate process"
        workflow_long_description : "Use this workflow for aligning with STAR and disambiguating human and mouse reads with AZ's Disambiguate process."
    }
}

task trim_reads {

    File r1
    File r2
    Int trim_length_bp

    String suffix="_R1.fastq.gz"

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1, suffix)

    Int disk_size = 200


    command {
        java -jar /opt/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
            -trimlog ${sample_name}.trim.log \
            -summary ${sample_name}.trim_summary.log \
            ${r1} ${r2} \
            -baseout ${sample_name}.trimmed.fastq.gz \
            CROP:{trim_length_bp}

        mv ${sample_name}.trimmed_1P.fastq.gz ${sample_name}_R1.fastq.gz
        mv ${sample_name}.trimmed_2P.fastq.gz ${sample_name}_R2.fastq.gz
    }

    output {
        File trimmed_r1 = "${sample_name}_R1.fastq.gz"
        File trimmed_r2 = "${sample_name}_R2.fastq.gz"
    }

    runtime {
        docker: "docker.io/hsphqbrc/star_disambig:v0.2"
        cpu: 4
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}