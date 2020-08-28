import "samtools.wdl" as samtools
import "feature_counts.wdl" as feature_counts

workflow DisambiguateBam {

    File human_bam
    File mouse_bam
    File human_bam_idx
    File mouse_bam_idx
    File human_gtf
    File mouse_gtf
    String human_tag
    String mouse_tag

    # the human-aligned + primary-filtered BAM should be <sample>.<human_tag>.primary_filtered.bam
    # the suffix is then ".<human_tag>.primary_filtered.bam"
    String suffix = "." + human_tag + ".primary_filtered.bam"

    String sample_name = basename(human_bam, suffix)

    call disambiguate_bams {
        input:
            human_bam = human_bam,
            mouse_bam = mouse_bam,
            human_bam_idx = human_bam_idx,
            mouse_bam_idx = mouse_bam_idx,
            sample_name = sample_name
    }

    call samtools.samtools_index as index_human_ambig {
        input:
            input_bam = disambiguate_bams.human_ambig_bam
    }

    call samtools.samtools_index as index_mouse_ambig {
        input:
            input_bam = disambiguate_bams.mouse_ambig_bam
    }

    call samtools.samtools_index as index_human_disambig {
        input:
            input_bam = disambiguate_bams.human_disambig_bam
    }

    call samtools.samtools_index as index_mouse_disambig {
        input:
            input_bam = disambiguate_bams.mouse_disambig_bam
    }

    # Quantify them
    call feature_counts.count_reads as quantify_human_ambig {
        input:
            input_bam = disambiguate_bams.human_ambig_bam,
            gtf = human_gtf,
            sample_name = sample_name,
            tag = "human_ambig"
    }

    call feature_counts.count_reads as quantify_mouse_ambig {
        input:
            input_bam = disambiguate_bams.mouse_ambig_bam,
            gtf = mouse_gtf,
            sample_name = sample_name,
            tag = "mouse_ambig"
    }

    call feature_counts.count_reads as quantify_human_disambig {
        input:
            input_bam = disambiguate_bams.human_disambig_bam,
            gtf = human_gtf,
            sample_name = sample_name,
            tag = "human_disambig"
    }
    
    call feature_counts.count_reads as quantify_mouse_disambig {
        input:
            input_bam = disambiguate_bams.mouse_disambig_bam,
            gtf = mouse_gtf,
            sample_name = sample_name,
            tag = "mouse_disambig"
    }

    output {
        File human_ambig_bam = disambiguate_bams.human_ambig_bam
        File mouse_ambig_bam = disambiguate_bams.mouse_ambig_bam
        File human_disambig_bam = disambiguate_bams.human_disambig_bam
        File mouse_disambig_bam = disambiguate_bams.mouse_disambig_bam

        File human_ambig_bam_index = index_human_ambig.bam_index
        File mouse_ambig_bam_index = index_mouse_ambig.bam_index
        File human_disambig_bam_index = index_human_disambig.bam_index
        File mouse_disambig_bam_index = index_mouse_disambig.bam_index

        File human_ambig_counts = quantify_human_ambig.count_output
        File human_ambig_counts_summary = quantify_human_ambig.count_output_summary
        File mouse_ambig_counts = quantify_mouse_ambig.count_output
        File mouse_ambig_counts_summary = quantify_mouse_ambig.count_output_summary
        File human_disambig_counts = quantify_human_disambig.count_output
        File human_disambig_counts_summary = quantify_human_disambig.count_output_summary
        File mouse_disambig_counts = quantify_mouse_disambig.count_output
        File mouse_disambig_counts_summary = quantify_mouse_disambig.count_output_summary
    }
}

task disambiguate_bams {

    File human_bam
    File mouse_bam
    File human_bam_idx
    File mouse_bam_idx
    String sample_name

    Int disk_size = 100

    command {
        ngs_disambiguate -o disambiguate_output -s ${sample_name} -a star ${human_bam} ${mouse_bam}
        mv disambiguate_output/${sample_name}.ambiguousSpeciesA.bam ${sample_name}.ambiguous.human.bam
        mv disambiguate_output/${sample_name}.disambiguatedSpeciesA.bam ${sample_name}.disambiguated.human.bam
        mv disambiguate_output/${sample_name}.ambiguousSpeciesB.bam ${sample_name}.ambiguous.mouse.bam
        mv disambiguate_output/${sample_name}.disambiguatedSpeciesB.bam ${sample_name}.disambiguated.mouse.bam 
    }

    output {
        File human_ambig_bam = "${sample_name}.ambiguous.human.bam"
        File human_disambig_bam = "${sample_name}.disambiguated.human.bam"
        File mouse_ambig_bam = "${sample_name}.ambiguous.mouse.bam"
        File mouse_disambig_bam = "${sample_name}.disambiguated.mouse.bam"
    }

    runtime {
        docker: "docker.io/hsphqbrc/star_disambig:v0.2"
        cpu: 2
        memory: "8 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

}
