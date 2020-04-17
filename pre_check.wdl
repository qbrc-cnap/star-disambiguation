workflow PairedRnaSeqAndQuantWorkflow {
    
    Array[File] r1_files
    Array[File] r2_files
    String genome
    File star_index_path
    File gtf
    String output_zip_name
    String git_repo_url
    String git_commit_hash


    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)
    scatter(item in fastq_pairs){

        call assert_valid_fastq {
            input:
                r1_file = item.left,
                r2_file = item.right
        }
    }
}

task assert_valid_fastq {

    File r1_file
    File r2_file
    Int disk_size = 100

    command <<<
        /usr/bin/python3 /opt/software/precheck/check_fastq.py -r1 ${r1_file} -r2 ${r2_file}
    >>>

    runtime {
        docker: "docker.io/hsphqbrc/star_disambig:v0.1"
        cpu: 4
        memory: "50 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
