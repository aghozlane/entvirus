#!/usr/bin/env nextflow

/*workflow.onComplete {
    def subject = 'Annot virus'
    def recipient = 'amine.ghozlane@pasteur.fr'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}*/

workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}


workflow.onError = {
    println "Oops .. something went wrong"
}

params.help=false

def usage() {
    println("entvirus --in <reads_dir> --out <output_dir> --cpus <nb_cpus> --mode <clc,megahit,metacompass,ray,spades> -w <temp_work_dir> --annotated <yes,no> --focus <yes,no>")
}


if(params.help){
    usage()
    exit(1)
}


params.in="$baseDir/test/"
readChannel = Channel.fromFilePairs("${params.in}/*_R{1,2}.{fastq,fastq.dsrc2,fastq.gz}")
                    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.in}" }
                    //.subscribe { println it }

/*contigsChannel = Channel
                .fromPath("${params.in}/*.fasta")
                .map {file -> tuple(file.baseName.replaceAll(".fasta",""),file)}*/
params.mail = "amine.ghozlane@pasteur.fr"
params.databases = "/pasteur/scratch/amine/polston_databases"
params.cpus = 2
params.vp1 = "${params.databases}/vp1_seq_nov_17.fasta"
params.ncbi = "${params.databases}/ncbi_viruses.fna"
params.rvdb = "${params.databases}/rVDBv10.2.fasta"
params.uniprot = "${params.databases}/uniprot_taxonomy.fasta"
params.uniref = "${params.databases}/uniref_uniprot.fasta"
params.viral = "${params.databases}/viral_catalogue_poltson.fna"
params.out = "$baseDir/annotation/"
//params.nt = "/local/databases/fasta/nt"
params.nt = "${params.databases}/nt"
params.gitaxidnucl = "/local/databases/release/taxodb/gi_taxid_nucl.dmp"
params.numberBestannotation = 10
params.coverage = 0
params.mismatch = 1
params.contaminant = "/local/databases/index/bowtie/2.1.0/hg19.fa"
params.alienseq = "${params.databases}/alienTrimmerPF8contaminants.fasta"
params.minlength = 45
params.cleaned_reads = "${params.out}/cleaned_reads"
params.khmer_reads = "${params.out}/khmer_reads"
params.blastdir = "${params.out}/blast"
params.mode = "clc"
params.filter = 1
params.vp1info = "/pasteur/scratch/amine/polston_databases/vp1_info_nov_17.tsv"
params.taxadb = "/pasteur/scratch/amine/polston_databases/taxadb_nucl.sqlite"
params.annotated =  "no"
params.focus = "no"
params.readlength = 150
params.vp1coverage = 50
inDir = file(params.in)

myDir = file(params.out)
myDir.mkdirs()

cleanDir = file("${params.cleaned_reads}")
cleanDir.mkdirs()

khmerDir = file("${params.khmer_reads}")
khmerDir.mkdirs()

blastDir = file("${params.blastdir}")
blastDir.mkdirs()


process filtering {
    //publishDir "$myDir", mode: 'copy'
    cpus params.cpus
    memory "8G"

    input:
    set pair_id, file(reads) from readChannel

    output:
    set pair_id, file("unmapped/*.1"), file("unmapped/*.2") into unmappedChannel
    //set pair_id, file("unmapped/*.1"), file("unmapped/*.2") into redundantChannel

    shell:
    """
    #!/usr/bin/env bash
    case "!{reads[0]}" in
    *.dsrc2 )
        dsrc d -t!{params.cpus} !{reads[0]} file_R1.fastq
        dsrc d -t!{params.cpus} !{reads[1]} file_R2.fastq
    ;;
    *.gz )
        gunzip !{reads[0]} -c > file_R1.fastq
        gunzip !{reads[1]} -c > file_R2.fastq
    ;;
    *)
        ln -s !{reads[0]} file_R1.fastq
        ln -s !{reads[1]} file_R2.fastq
    ;;
    esac
    mkdir unmapped
    if [ !{params.focus} ==  "no" ]
    then
        bowtie2 -q -N !{params.mismatch} -1 file_R1.fastq -2 file_R2.fastq \
                -x !{params.contaminant} --un-conc unmapped/ -S /dev/null \
                -p !{params.cpus} --very-sensitive-local
    else
        bowtie2 -q -N !{params.mismatch} -1 file_R1.fastq -2 file_R2.fastq \
                -x !{params.viral} --al-conc unmapped/ -S /dev/null \
                -p !{params.cpus} --very-sensitive-local
    fi
    """
}


process trimming {
    //publishDir "$myDir", mode: 'copy'
    cpus params.cpus

    input:
    set pair_id, file(forward), file(reverse) from unmappedChannel


    output:
    set pair_id, file("*_R1.fastq"), file("*_R2.fastq") into trimChannel
    file("*_R1.fastq") into r1Channel
    file("*_R2.fastq") into r2Channel
    file("*.fastq.gz") into mappingChannel mode flatten

    script:
    """
    AlienTrimmer -if ${forward} -ir ${reverse} -of ${pair_id}_R1.fastq \
                 -or ${pair_id}_R2.fastq -os ${pair_id}_sgl.fastq \
                 -c ${params.alienseq} -l ${params.minlength}
    gzip -c ${pair_id}_R1.fastq >  ${pair_id}_R1.fastq.gz
    gzip -c ${pair_id}_R2.fastq > ${pair_id}_R2.fastq.gz
    """

}

mappingChannel.subscribe { it.copyTo(cleanDir) }

/*filterChannel = Channel.create()
nextChannel = Channel.create()
if( params.filter == 0 ) {
    trimChannel.into{ nextChannel }
    filterChannel.close()
    trimChannel.close()
}
else{
    trimChannel.into{ filterChannel }
    nextChannel.close()
    trimChannel.close()
}*/

process khmer {
    cpus params.cpus
    //memory "25G"
    memory "15G"

    input:
    //set pair_id, file(forward), file(reverse) from filterChannel
    set pair_id, file(forward), file(reverse) from trimChannel

    output:
    set pair_id, file("khmer/*_R1.fastq"), file("khmer/*_R2.fastq") into khmerChannel
    file("khmer/*.fastq.gz") into khmeroutChannel mode flatten

    shell:
    """
    #!/usr/bin/env bash
    mkdir khmer
    interleave-reads.py ${forward} ${reverse} --output interleaved.pe
    normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct \
        interleaved.pe --output output.pe.keep
    filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter \
        -T ${params.cpus}
    extract-paired-reads.py output.pe.filter --output-paired output.dn.pe \
        --output-single output.dn.se
    split-paired-reads.py output.dn.pe -1 khmer/${pair_id}_khmer_R1.fastq \
        -2 khmer/${pair_id}_khmer_R2.fastq
    gzip -c khmer/${pair_id}_khmer_R1.fastq >  khmer/${pair_id}_R1.fastq.gz
    gzip -c khmer/${pair_id}_khmer_R2.fastq > khmer/${pair_id}_R2.fastq.gz
    """
}

khmeroutChannel.subscribe{ it.copyTo(khmerDir) }

/*readtoassemblyChannel = Channel.create()
assemblyChannel = Channel.create()
readtoassemblyChannel.mix(nextChannel, khmerChannel)
                     .into{ assemblyChannel }*/
/*mappingChannel = Channel.create()
assemblyChannel = Channel.create()
khmerChannel.into { mappingChannel ; assemblyChannel }*/

/*Channel
  .fromFilePairs('*-{1,2}.txt')
    .flatMap { key, files -> [[key], files].combinations() }
      .collectFile()*/


process assembly {
    publishDir "$myDir", mode: 'copy'
    //memory "20G"
    memory "5G"
    //cpus params.cpus

    if(params.mode == "clc"){
        clusterOptions='--qos=normal -C clcbio -p common'
        //clusterOptions='--qos=clcgwb --x11 clcgenomicswb9'
        cpus params.cpus
    }
    else if(params.mode == "metacompass"){
        beforeScript ='source /local/gensoft2/adm/etc/profile.d/modules.sh;module use /pasteur/projets/policy01/Matrix/modules;export PATH=/pasteur/projets/policy01/Matrix/metagenomics/entvirus/bin/kmer-code-2013-trunk/Linux-amd64/bin/:$PATH'
        module = 'Python/3.6.0:samtools/1.3:snakemake/3.5.4:bowtie2/2.2.9'
        cpus params.cpus
    }
    else if(params.mode == "ray"){
        cpus 1
    }
    //else if(params.mode == "minia"){
    //    beforeScript ='source /local/gensoft2/adm/etc/profile.d/modules.sh;module use /pasteur/projets/policy01/Matrix/modules;export PYTHONPATH=/pasteur/projets/policy01/Matrix/metagenomics/python-lib/lib/python2.7/site-packages;export LD_LIBRARY_PATH=/pasteur/projets/policy01/Matrix/metagenomics/htslib/lib'
    //    module = 'Python/2.7.8:bwa/0.7.7'
    //}
    else{
        cpus params.cpus
    }

    input:
    //set pair_id, file(forward), file(reverse) from assemblyChannel
    set pair_id, file(forward), file(reverse) from khmerChannel
    //set pair_id, file(forward), file(reverse) from trimChannel

    output:
    set pair_id, file("assembly/*_{clc,megahit,metacompass,ray,spades}.fasta") into contigsChannel

    shell:
    """
    #!/usr/bin/env bash
    if [ !{params.mode} ==  "clc" ]
    then
        mkdir assembly
        clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q \
            -i !{forward} !{reverse} --cpus !{params.cpus}
        python !{baseDir}/bin/rename_fasta.py -i assembly/contigs.fasta \
                -s !{pair_id} -o assembly/!{pair_id}_clc.fasta
    elif [ !{params.mode} ==  "metacompass" ]
    then
        python3 !{baseDir}/bin/MetaCompass/go_metacompass.py -P !{forward},!{reverse}  \
            -t !{params.cpus} -o assembly/   -l !{params.readlength} --clobber \
            -r !{params.viral}
            python !{baseDir}/bin/rename_fasta.py -i assembly/metacompass.final.ctg.fa \
            -o assembly/!{pair_id}_metacompass.fasta -s !{pair_id}
    elif [ !{params.mode} ==  "ray" ]
    then
        Ray -k 31 -p !{forward} !{reverse} -o assembly/
        python !{baseDir}/bin/rename_fasta.py -i assembly/Contigs.fasta \
            -o assembly/!{pair_id}_ray.fasta -s !{pair_id}
    elif [ !{params.mode} ==  "megahit" ]
    then
        megahit -1 !{forward} -2 !{reverse} -o assembly/ -t !{params.cpus}
        python !{baseDir}/bin/rename_fasta.py -i assembly/final.contigs.fa \
            -o assembly/!{pair_id}_megahit.fasta -s !{pair_id}
    else
        mkdir assembly
        if [ !{params.focus} == "no" ]
        then
            spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpus} -o assembly/
        else
            spades.py -1 !{forward} -2 !{reverse} -t !{params.cpus} -o assembly/
        fi
        python !{baseDir}/bin/rename_fasta.py -i assembly/contigs.fasta \
                -s !{pair_id} -o assembly/!{pair_id}_spades.fasta
    fi
    """
}



// Perform all blast annotation
process blast {
    //publishDir "$myDir", mode: 'copy'
    //clusterOptions='--qos=normal -p common'
    cpus params.cpus
    //memory "25G"
    memory "5G"

    input:
    set contigsID, file(contigs) from contigsChannel

    output:
    //file("log.txt") into logChannel
    set contigsID, file(contigs), file("blast/*_vp1.tsv") into vp1blastChannel mode flatten
    //set contigsID, file(contigs), file("blast/*_vp1_vsearch.tsv") into vp1vsearchChannel mode flatten
    //file("blast_vp1/*.tsv") into vp1blasttocombineChannel
    set contigsID, file("blast/*_ncbi_all.tsv") into blastChannel mode flatten
    file("blast/*.tsv") into allblastChannel mode flatten

    shell:
    """
    mkdir blast
    blastn -query !{contigs}  -db !{params.vp1}  -num_threads !{params.cpus} \
           -out blast/!{contigsID}_vp1.tsv -max_target_seqs 1 \
           -outfmt '6 qseqid sseqid qlen length qstart qend sstart send pident qcovs evalue'\
           -task blastn -evalue 1E-3
    #vsearch --usearch_global !{contigs} --blast6out blast/!{contigsID}_vp1_vsearch.tsv --db !{params.vp1} \
    #       --id 0.5  --top_hits_only --userfields query+target+id+qrow  \
    #       --userout blast/!{contigsID}_vp1_vsearch.tsv
    blastn -query !{contigs}  -db !{params.viral}  -num_threads !{params.cpus} \
           -out blast/!{contigsID}_polston.tsv \
           -max_target_seqs !{params.numberBestannotation} -use_index true \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    blastn -query !{contigs}  -db !{params.ncbi}  -num_threads !{params.cpus} \
           -out blast/!{contigsID}_ncbi.tsv \
           -max_target_seqs !{params.numberBestannotation} -use_index true \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    blastn -query !{contigs}  -db !{params.rvdb}  -num_threads !{params.cpus} \
           -max_target_seqs !{params.numberBestannotation} -use_index true \
           -out blast/!{contigsID}_rvdb.tsv\
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    blastx -query !{contigs}  -db !{params.uniprot} -num_threads !{params.cpus}\
           -out blast/!{contigsID}_uniprot.tsv \
           -max_target_seqs !{params.numberBestannotation} \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    blastx -query !{contigs}  -db !{params.uniref}  -num_threads !{params.cpus}\
           -out blast/!{contigsID}_uniref.tsv \
           -max_target_seqs !{params.numberBestannotation} \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    blastn -query !{contigs}  -db !{params.nt}  -num_threads !{params.cpus} \
           -out blast/!{contigsID}_nt.tsv \
           -max_target_seqs !{params.numberBestannotation} \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    cat blast/!{contigsID}_ncbi.tsv blast/!{contigsID}_nt.tsv \
        blast/!{contigsID}_rvdb.tsv > blast/!{contigsID}_ncbi_all.tsv
    """
}

allblastChannel.subscribe { it.copyTo(blastDir) }


// Extract vp1 sequences
process vp1 {
    publishDir "$myDir", mode: 'copy'

    input:
    set contigsID, file(contigs), file(vp1blast) from vp1blastChannel

    output:
    file("vp1/*_vp1.fasta") into vp1Channel
    file("vp1_contigs/*_vp1_contigs.fasta") into vp1contigsChannel

    shell:
    """
    #!/usr/bin/env bash
    mkdir vp1 vp1_contigs
    nlines=\$(wc -l !{vp1blast} |cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
        # Extract VP1 sequence
        python !{baseDir}/bin/extract_sequence.py -q !{contigs} -b !{vp1blast} \
               -o vp1/!{contigsID}_vp1.fasta -a !{params.vp1info} \
               -c !{params.vp1coverage}
        # Extract VP1 contig
        python !{baseDir}/bin/grab_catalogue_sequence.py -i !{vp1blast} \
               -d !{contigs} -o vp1_contigs/!{contigsID}_vp1_contigs.fasta \
               -a !{params.vp1info} -v !{params.vp1coverage}
    else
        touch vp1_contigs/!{contigsID}_vp1_contigs.fasta vp1/!{contigsID}_vp1.fasta
    fi
    """
}

/*process vp1_vsearch {
    publishDir "$myDir", mode: 'copy'

    input:
    set contigsID, file(contigs), file(vp1vsearch) from vp1vsearchChannel

    output:
    file("vp1/*_vp1_vsearch.fasta") into vp1vChannel
    file("vp1_contigs/*_vp1_contigs_vsearch.fasta") into vp1contigsvChannel

    shell:
    """
    #!/usr/bin/env bash
    mkdir vp1 vp1_contigs
    nlines=\$(wc -l !{vp1blast} |cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
        # Extract VP1 sequence
        python !{baseDir}/bin/extract_sequence_vsearch.py -i !{vp1vsearch} \
               -o vp1/!{contigsID}_vp1_vsearch.fasta -a !{params.vp1info}
        # Extract VP1 contig
        python !{baseDir}/bin/grab_catalogue_sequence.py -i !{vp1vsearch} \
               -d !{contigs} -o vp1_contigs/!{contigsID}_vp1_contigs_vsearch.fasta
    else
        touch vp1_contigs/!{contigsID}_vp1_contigs_vsearch.fasta vp1/!{contigsID}_vp1_vsearch.fasta
    fi
    """
}*/

// Extract ncbi tree
process annotation {
    //clusterOptions='--qos=normal -p common'
    publishDir "$myDir", mode: 'copy'
    //memory "10G"

    beforeScript ='source /local/gensoft2/adm/etc/profile.d/modules.sh;module use /pasteur/projets/policy01/Matrix/modules;export PYTHONPATH=/pasteur/projets/policy01/Matrix/metagenomics/python-lib//lib/python3.6/site-packages/'
    module = 'Python/3.6.0'

    input:
    set contigsID, file(ncbi_annotation) from blastChannel

    output:
    file("annotation/*_annotation.tsv") into annotationChannel mode flatten
    file("annotation/*_taxonomy.tsv") into taxChannel mode flatten

    shell:
    """
    #!/usr/bin/env bash
    mkdir annotation
    nlines=\$(wc -l !{ncbi_annotation}|cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
        # Get the taxonomy
        python3 !{baseDir}/bin/get_taxonomy2.py -i !{ncbi_annotation} \
                -t !{params.gitaxidnucl} -d !{params.taxadb} \
                -o annotation/!{contigsID}_taxonomy.tsv
        python2 !{baseDir}/bin/ExtractNCBIDB.py -f !{ncbi_annotation} \
                -g annotation/!{contigsID}_taxonomy.tsv \
                -o annotation/!{ncbi_annotation.baseName}_annotation.tsv \
                -nb !{params.numberBestannotation} -fc !{params.coverage} -fi
    else
        touch annotation/!{ncbi_annotation.baseName}_annotation.tsv
        touch annotation/!{ncbi_annotation.baseName}_taxonomy.tsv
    fi
    """
}

combineChannel = annotationChannel.collectFile(name: 'combined.tsv')
//vp1contigs_vsearch = vp1contigsvChannel.collectFile(name: 'vp1contigs_vsearch.fasta')
vp1contigs = vp1contigsChannel.collectFile(name: 'vp1contigs.fasta')
//vp1annotation = vp1blasttocombineChannel.collectFile(name: 'vp1blast.tsv')

process buildIndex {
    cpus params.cpus
    //publishDir "$myDir", mode: 'copy'
    //clusterOptions='--qos=normal -p common'

    input:
    file vp1contigs

    output:
    file("vp1contigs.index*") into vp1contigs_index
    file(vp1contigs) into vp1contigs_fasta

    script:
    """
    bowtie2-build --threads ${params.cpus} ${vp1contigs} vp1contigs.index
    """
}

process abundance_vp1 {
    //queue = 'common'
    //clusterOptions='--qos=normal -p common'
    publishDir "$myDir", mode: 'copy'
    //errorStrategy 'finish'

    input:
    file(vp1contigs) from vp1contigs_index.first()
    file(forward) from r1Channel.toList()
    file(reverse) from r2Channel.toList()
    //file cleanDir


    output:
    file("abundance/count_matrix.tsv") into countChannel
    // change to shared
    shell:
    """
    mbma.py mapping --r1 !{forward} --r2 !{reverse} -o abundance \
           -db vp1contigs.index -e !{params.mail} -q !{params.queue} \
           -p !{params.partition} --bowtie2 \
           --shared -m PE -t !{params.cpus}
    mv abundance/comptage/count_matrix.txt abundance/count_matrix.tsv
    """
}


//process combine_vp1 {
    //input:
    //file(vp1) from vp1blasttocombineChannel.toList()

    //shell:
    //"""
    //python3 !{baseDir}/bin/extract_result.py -i params.in summary
    //"""
//}

process combine_annotation {
    publishDir "$myDir", mode: 'copy'

    input:
    file(ncbi_annotation) from combineChannel
    file(vp1contigs_fasta) from vp1contigs_fasta

    output:
    file("contigs_annotation.tsv") into finalChannel mode flatten
    file("vp1contigs_annotation.tsv") into vp1finalChannel
    file(vp1contigs_fasta) into vp1contigsfasta

    shell:
    """
    #!/bin/bash
    nlines=\$(wc -l !{ncbi_annotation}|cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
        # Combine the annotation
        python !{baseDir}/bin/combine_annotation.py -i !{ncbi_annotation} \
            -o contigs_annotation.tsv
        nannot=\$(wc -l contigs_annotation.tsv|cut -f 1 -d ' ')
        if [ \${nannot} -gt '0' ]
        then
            # Get All the vp1
            grep "^>"  !{vp1contigs_fasta} |cut -f 2 -d ">" > list
            nvp1=\$(wc -l list|cut -f 1 -d ' ')
            if [ \${nvp1} -gt '0' ]
            then
                # Extract vp1 annotation
                python !{baseDir}/bin/strange.py -i contigs_annotation.tsv \
                    -c list -o vp1contigs_annotation.tsv
            else
                touch vp1contigs_annotation.tsv
            fi
        else
            touch vp1contigs_annotation.tsv
        fi
    else
        touch vp1contigs_annotation.tsv contigs_annotation.tsv
    fi
    """
}


//contigs_annotation = finalChannel.collectFile(name: 'contigs_annotation.tsv')
//contigs_annotation.subscribe { it.copyTo(myDir) }
//vp1contigs_annotation = vp1finalChannel.collectFile(name: 'vp1contigs_annotation.tsv')
//vp1contigs_annotation.subscribe { it.copyTo(myDir) }

process summary {
    //publishDir "$myDir", mode: 'copy'

    input:
    file(vp1contigs_annot) from vp1finalChannel
    file myDir
    file inDir
    file countChannel

    output:
    file("result_summary.tsv") into resultChannel
    //file(vp1contigs_annot) into overChannel

    script:
    """
    python ${baseDir}/bin/extract_result.py -i $myDir -a ${params.vp1info}\
             -vp1 ${vp1contigs_annot} -o result_summary.tsv -r $inDir\
             -c ${countChannel} -l ${params.annotated} -v ${params.vp1coverage}
    """
}

resultChannel.subscribe { it.copyTo(myDir) }
println "Project : $workflow.projectDir"
//println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
