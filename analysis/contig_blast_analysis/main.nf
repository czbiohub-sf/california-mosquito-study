
def helpMessage() {
    log.info """
    ==============================================================
      _                            _       _ _         _ _
     | |_ ___ _____ ___ ___    ___|_|_____|_| |___ ___|_| |_ _ _
     | '_|___|     | -_|  _|  |_ -| |     | | | .'|  _| |  _| | |
     |_,_|   |_|_|_|___|_|    |___|_|_|_|_|_|_|__,|_| |_|_| |_  |
                                                            |___|
    ==============================================================

    Usage:

    The typical command for running the pipeline is as follows.

    With a samples.csv file containing the columns sample_id,read1,read2:

      nextflow run czbiohub/nf-kmer-similarity \
        --outdir s3://olgabot-maca/nf-kmer-similarity/ --samples samples.csv


    With read pairs in one or more semicolon-separated s3 directories:

      nextflow run czbiohub/nf-kmer-similarity \
        --outdir s3://olgabot-maca/nf-kmer-similarity/ \
        --read_pairs s3://olgabot-maca/sra/homo_sapiens/smartseq2_quartzseq/*{R1,R2}*.fastq.gz;s3://olgabot-maca/sra/danio_rerio/smart-seq/whole_kidney_marrow_prjna393431/*{R1,R2}*.fastq.gz


    With plain ole fastas in one or more semicolon-separated s3 directories:

      nextflow run czbiohub/nf-kmer-similarity \
        --outdir s3://olgabot-maca/nf-kmer-similarity/choanoflagellates_richter2018/ \
        --fastas /home/olga/data/figshare/choanoflagellates_richter2018/1_choanoflagellate_transcriptomes/*.fasta


    With SRA ids (requires nextflow v19.03-edge or greater):

      nextflow run czbiohub/nf-kmer-similarity \
        --outdir s3://olgabot-maca/nf-kmer-similarity/ --sra SRP016501


    Mandatory Arguments:
      --outdir                      Local or S3 directory to output the comparison matrix to

    Sample Arguments -- One or more of:
      --samples                     CSV file with columns id, read1, read2 for each sample
      --fastas
      --read_pairs                 Local or s3 directories containing *R{1,2}*.fastq.gz
                                    files, separated by commas
      --sra                         SRR, ERR, SRP IDs representing a project. Only compatible with
                                    Nextflow 19.03-edge or greater


    Options:
      --ksizes                      Which nucleotide k-mer sizes to use. Multiple are
                                    separated by commas. Default is '21,27,33,51'
      --molecules                   Which molecule to compare on. Default is both DNA
                                    and protein, i.e. 'dna,protein'
      --log2_sketch_sizes           Which log2 sketch sizes to use. Multiple are separated
                                    by commas. Default is '10,12,14,16'
      --one_signature_per_record    Make a k-mer signature for each record in the FASTQ/FASTA files.
                                    Useful for comparing e.g. assembled transcriptomes or metagenomes.
                                    (Not typically used for raw sequencing data as this would create
                                    a k-mer signature for each read!)
    """.stripIndent()
}



// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


 // R1, R2 pairs from a samples.csv file
 samples_ch = Channel.empty()
 
 // Provided a samples.csv file
 Channel
  .fromPath("s3://czbiohub-mosquito/contigs/*", type:"dir")
  .map{ f -> tuple(f.name, file(f))}
  //.println()
  .set{ samples_ch }

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

process lca_analysis {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'lucymli/lca_analysis'

    // If job fails, try again with more memory
    // memory { 8.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3

    input:
    set sample_id, file(sample_dir) from samples_ch

    script: // ADD PYTHON FILE TO /usr/local/bin IN DOCKER FILE
    tsv = sample_dir
    tsv_name = sample_id
    filtered_tsv = tsv.getSimpleName()+"_filtered.m9"
    lca_tsv = tsv.getName().replaceAll("blast", "lca")
    """
    python lca_analysis.py \
    --blast_type nt \
    --fpath ${tsv_name}/blast_nt.m9 \
    --filtered_blast_path $filtered_tsv \
    --outpath $lca_tsv \
    --ident_cutoff 0.5 \
    --align_len_cutoff 0.5
    """
}

// // sourmash_sketches.println()
// // sourmash_sketches.groupTuple(by: [0,3]).println()

// process sourmash_compare_sketches {
//     tag "${sketch_id}"

//     container 'czbiohub/nf-kmer-similarity'
//     publishDir "${params.outdir}/", mode: 'copy'
//     errorStrategy 'retry'
//   maxRetries 3

//     input:
//   set val(sketch_id), val(molecule), val(ksize), val(log2_sketch_size), file ("sketches/*.sig") \
//     from sourmash_sketches.groupTuple(by: [0, 3])

//     output:
//     file "similarities_${sketch_id}.csv"

//     script:
//     """
//     sourmash compare \
//         --ksize ${ksize[0]} \
//         --${molecule[0]} \
//         --csv similarities_${sketch_id}.csv \
//         --traverse-directory .
//     """

// }
