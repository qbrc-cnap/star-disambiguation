workflow feature_counts_test {
    File input_bam
    File gtf
    String sample_name
    String tag
    Array[File] countfile_array

    call count_reads {
        input:
            input_bam = input_bam,
            gtf = gtf,
            sample_name = sample_name,
            tag = tag
    }

    call concatenate {
        input:
            count_files = countfile_array,
            output_filename = "test_output_matrix.tsv"
    }

    output {
        File countfile = count_reads.count_output
        File countmatrix = concatenate.count_matrix
    }
}


task count_reads {

    File input_bam
    File gtf
    String sample_name
    String tag

    String output_counts_name = sample_name + "." + tag + ".feature_counts.tsv"

    String strand_option = "0"

    Int disk_size = 200

    command {
        featureCounts \
            -s${strand_option} \
            -p \
            -t exon \
            -g gene_name \
            -a ${gtf} \
            -o ${output_counts_name} \
            ${input_bam}
    }

    output {
        File count_output = "${output_counts_name}"
        File count_output_summary = "${output_counts_name}.summary"
    }

    runtime {
        docker: "docker.io/hsphqbrc/star_disambig:v0.1"
        cpu: 8
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task concatenate {
    # This concatenates the featureCounts count files into a 
    # raw count matrix.

    Array[File] count_files
    String output_filename

    Int disk_size = 30

    command {
        concatenate_featurecounts.py -o ${output_filename} ${sep=" " count_files}
    }

    output {
        File count_matrix = "${output_filename}"
    }

    runtime {
        docker: "docker.io/hsphqbrc/star_disambig:v0.1"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0 
    }

}
