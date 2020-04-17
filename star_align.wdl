workflow align_test {
    File r1_fastq
    File? r2_fastq
    File star_index_path
    File gtf
    String sample_name

    call perform_align {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            star_index_path = star_index_path,
            gtf = gtf,
            sample_name = sample_name
    }

    output {
        File bam = perform_align.sorted_bam
    }
}

task perform_align{
    # align to the reference genome using the STAR aligner
    # The STAR alignment produces a position-sorted BAM file

    # Input params passed by a parent Workflow:
    # We require that there is a single R1 fastq and possibly a single
    # R2 fastq for paired sequencing.
    # Genome is a string that matches one of the keys in the Map below
    # sample_name helps with naming files for eventual concatenation.
    File r1_fastq
    File? r2_fastq
    File star_index_path
    File gtf
    String sample_name

    # Default disk size in GB
    Int disk_size = 500

    command {
        set -euxo pipefail
        mkdir -p workspace/index
        tar -xf ${star_index_path} -C workspace/index
        STAR \
            --readFilesIn ${r1_fastq} ${r2_fastq} \
            --genomeDir workspace/index \
            --outFileNamePrefix "${sample_name}." \
            --twopassMode Basic \
            --runThreadN 16 \
            --readFilesCommand zcat \
            --sjdbGTFfile ${gtf} \
            --outFilterType BySJout \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx
    }

    output {
        File sorted_bam = "${sample_name}.Aligned.sortedByCoord.out.bam"
        File run_log = "${sample_name}.Log.out"
        File final_log = "${sample_name}.Log.final.out"
        File unmapped_mate1= "${sample_name}.Unmapped.out.mate1"
        File? unmapped_mate2 = "${sample_name}.Unmapped.out.mate2"
    }

    runtime {
        docker: "docker.io/hsphqbrc/star_disambig:v0.1"
        cpu: 8
        memory: "40 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

}
