workflow samtools_test {
    File input_bam
    String sample_name

    call samtools_index {
        input:
            input_bam = input_bam
    }

    call samtools_primary_filter {
        input:
            input_bam = input_bam,
            input_bam_index = samtools_index.bam_index,
            sample_name = sample_name
    }

    output {
        File primary_filtered_bam = samtools_primary_filter.output_bam
    }
}


task samtools_index {

    File input_bam

    String bam_index_name = basename(input_bam) + ".bai"


    Int disk_size = 200

    command {
        samtools index ${input_bam} "${bam_index_name}"
    }

    output {
        File bam_index = "${bam_index_name}"
    }

    runtime {
        docker: "docker.io/hsphqbrc/star_disambig:v0.2"
        cpu: 8
        memory: "16 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

}

task samtools_primary_filter {

    File input_bam
    File input_bam_index
    String sample_name

    String output_bam_name = sample_name + ".primary_filtered.bam"

    Int disk_size = 200

    command {
        samtools view -b -F 0x0100 ${input_bam} > "${output_bam_name}"
    }

    output {
        File output_bam = "${output_bam_name}"
    }

    runtime {
        docker: "docker.io/hsphqbrc/star_disambig:v0.2"
        cpu: 8
        memory: "16 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
