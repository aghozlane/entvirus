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
    println("entvirus --in <reads_dir> --out <output_dir> --cpus <nb_cpus> --mode <clc,megahit,metacompass,ray,spades> -w <temp_work_dir> --annotated <yes,no> --focus <yes,no> --abundance <yes,no>")
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
params.cpus = 2
params.cpus_phylogeny = 12
params.seq3d = "$baseDir/databases/EVABCD_3D_2018_04.fa"
params.seq5utr = "$baseDir/databases/EVABCD_5UTR_UG52_UC53_2018_04.fa"
params.vp1 = "$baseDir/databases/SEQ-EV-vp1-2018-03.fasta"
params.p1 = "$baseDir/databases/SEQ-EV-P1-2018-03.fa"
params.seqfull = "$baseDir/databases/EV_ABCD_full_2018_04.fa"
params.viral = "$baseDir/databases/EV_ABCD_full_2018_04.fa"
//params.viral = "$baseDir/databases/viral_catalogue_poltson.fna"
params.out = "$baseDir/annotation/"
params.nt = "/local/databases/fasta/nt"
//params.nt = "/pasteur/projets/policy01/Biomics/metagenomics/catalogue/nt"
params.gitaxidnucl = "/local/databases/release/taxodb/gi_taxid_nucl.dmp"
params.numberBestannotation = 10
params.coverage = 0
params.mismatch = 1
params.contaminant = "/local/databases/index/bowtie/2.1.0/hg19.fa"
params.alienseq = "$baseDir/databases/alienTrimmerPF8contaminants.fasta"
params.minlength = 45
params.cleaned_reads = "${params.out}/cleaned_reads"
params.khmer_reads = "${params.out}/khmer_reads"
params.blastdir = "${params.out}/blast"
params.mode = "clc"
params.filter = 1
params.info3d = "$baseDir/databases/EVABCD_3D_2018_04.tsv"
params.info5utr = "$baseDir/databases/EVABCD_5UTR_UG52_UC53_2018_04.tsv"
params.vp1info = "$baseDir/databases/SEQ-EV-vp1-2018-03.tsv"
params.p1info = "$baseDir/databases/SEQ-EV-P1-2018-03.tsv"
params.taxadb = "/local/databases/rel/taxadb/current/db/taxadb_full.sqlite"
params.annotated =  "no"
params.focus = "no"
params.readlength = 150
params.vp1coverage = 50
params.p1coverage = 50
params.coverage5utr = 50
params.coverage3d = 50
params.abundance = "no"
params.evalue = 1E-3
params.filter_matrix = 100
params.ent_serotype = "$baseDir/databases/enterovirus_species_type-2018-04.csv"
params.conserved_position = 0.8
params.root_seq_full = "$baseDir/databases/RV-A1B-B632-D00239.fa"
params.root_seq_p1 = "$baseDir/databases/RV-A1B-B632-D00239-P1.fa"
params.root_seq_vp1 = "$baseDir/databases/RV-A1B-B632-D00239-VP1.fa"
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
    file("raw/*_R1.fastq") into r1Channelresume
    file("raw/*_R2.fastq") into r2Channelresume
    //set pair_id, file("unmapped/*.1"), file("unmapped/*.2") into redundantChannel

    shell:
    """
    #!/usr/bin/env bash
    mkdir raw
    case "!{reads[0]}" in
    *.dsrc2 )
        dsrc d -t!{params.cpus} !{reads[0]} raw/${pair_id}_R1.fastq
        dsrc d -t!{params.cpus} !{reads[1]} raw/${pair_id}_R2.fastq
    ;;
    *.gz )
        gunzip -c !{reads[0]} > raw/!{pair_id}_R1.fastq
        gunzip -c !{reads[1]} > raw/!{pair_id}_R2.fastq
    ;;
    *)
        mv !{reads[0]} raw/!{pair_id}_R1.fastq
        mv !{reads[1]} raw/!{pair_id}_R2.fastq
    ;;
    esac
    mkdir unmapped
    if [ !{params.focus} ==  "no" ]
    then
        cp -L raw/!{pair_id}_R1.fastq unmapped/!{pair_id}.1
        cp -L raw/!{pair_id}_R2.fastq unmapped/!{pair_id}.2
        #bowtie2 -q -N !{params.mismatch} -1 raw/!{pair_id}_R1.fastq -2 raw/!{pair_id}_R2.fastq \
        #        -x !{params.contaminant} --un-conc unmapped/ -S /dev/null \
        #        -p !{params.cpus} --very-sensitive-local
    else
        bowtie2 -q -N !{params.mismatch} -1 raw/!{pair_id}_R1.fastq -2 raw/!{pair_id}_R2.fastq \
                -x !{params.viral} --al-conc unmapped/ -S /dev/null \
                -p !{params.cpus} --very-sensitive-local
    fi
    """
}


process trimming {
    //publishDir "$myDir", mode: 'copy'

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
                 -or ${pair_id}_R2.fastq \
                 -c ${params.alienseq} -l ${params.minlength}
    gzip -c ${pair_id}_R1.fastq > ${pair_id}_R1.fastq.gz
    gzip -c ${pair_id}_R2.fastq > ${pair_id}_R2.fastq.gz
    """

}

mappingChannel.subscribe { it.copyTo(cleanDir) }


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
    file("khmer/*_khmer_R1.fastq.gz") into procr1Channelresume
    file("khmer/*_khmer_R2.fastq.gz") into procr2Channelresume

    shell:
    """
    #!/usr/bin/env bash
    mkdir khmer
    interleave-reads.py !{forward} !{reverse} --output interleaved.pe
    normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savegraph graph.ct \
        interleaved.pe --output output.pe.keep
    filter-abund.py -V graph.ct output.pe.keep --output output.pe.filter \
        -T ${params.cpus}
    extract-paired-reads.py output.pe.filter --output-paired output.dn.pe \
        --output-single output.dn.se
    split-paired-reads.py output.dn.pe -1 khmer/!{pair_id}_khmer_R1.fastq \
        -2 khmer/!{pair_id}_khmer_R2.fastq
    pigz -c -p !{params.cpus} khmer/!{pair_id}_khmer_R1.fastq >  khmer/!{pair_id}_khmer_R1.fastq.gz
    pigz -c -p !{params.cpus} khmer/!{pair_id}_khmer_R2.fastq > khmer/!{pair_id}_khmer_R2.fastq.gz
    """
}

khmeroutChannel.subscribe{ it.copyTo(khmerDir) }


process assembly {
    publishDir "$myDir", mode: 'copy'
    //memory "20G"
    memory "5G"
    //cpus params.cpus

    if(params.mode == "clc"){
        //clusterOptions='--qos=normal -C clcbio -p common'
        clusterOptions='--qos=clcbio -p clcbio'
        cpus params.cpus
    }
    /*else if(params.mode == "metacompass"){
        beforeScript ='source /local/gensoft2/adm/etc/profile.d/modules.sh;module use /pasteur/projets/policy01/Matrix/modules;export PATH=/pasteur/projets/policy01/Matrix/metagenomics/entvirus/bin/kmer-code-2013-trunk/Linux-amd64/bin/:$PATH'
        module = 'Python/3.6.0:samtools/1.3:snakemake/3.5.4:bowtie2/2.2.9'
        cpus params.cpus
    }*/
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
    file("assembly/*_{clc,megahit,metacompass,ray,spades}.fasta") into contigsChannelresume

    shell:
    """
    #!/usr/bin/env bash
    if [ !{params.mode} ==  "clc" ]
    then
        mkdir assembly
        clc_assembler -o assembly/contigs.fasta -p fb ss 180 250 -q \
            -i !{forward} !{reverse} --cpus !{params.cpus}
        rename_fasta.py -i assembly/contigs.fasta \
                -s !{pair_id} -o assembly/!{pair_id}_clc.fasta
    elif [ !{params.mode} ==  "metacompass" ]
    then
        python3 !{baseDir}/bin/MetaCompass/go_metacompass.py -P !{forward},!{reverse}  \
            -t !{params.cpus} -o assembly/   -l !{params.readlength} --clobber \
            -r !{params.viral}
            rename_fasta.py -i assembly/metacompass.final.ctg.fa \
            -o assembly/!{pair_id}_metacompass.fasta -s !{pair_id}
    elif [ !{params.mode} ==  "ray" ]
    then
        Ray -k 31 -p !{forward} !{reverse} -o assembly/
        rename_fasta.py -i assembly/Contigs.fasta \
            -o assembly/!{pair_id}_ray.fasta -s !{pair_id}
    elif [ !{params.mode} ==  "megahit" ]
    then
        megahit -1 !{forward} -2 !{reverse} -o assembly/ -t !{params.cpus}
        rename_fasta.py -i assembly/final.contigs.fa \
            -o assembly/!{pair_id}_megahit.fasta -s !{pair_id}
    else
        mkdir assembly
        if [ !{params.focus} == "no" ]
        then
            spades.py --meta -1 !{forward} -2 !{reverse} -t !{params.cpus} -o assembly/
        else
            spades.py -1 !{forward} -2 !{reverse} -t !{params.cpus} -o assembly/
        fi
        rename_fasta.py -i assembly/scaffolds.fasta \
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
    set contigsID, file(contigs), file("*_vp1.tsv") into vp1blastChannel
    set contigsID, file(contigs), file("*_p1.tsv"), file("*_5utr.tsv"), file("*_3d.tsv") into otherblastChannel
    //file("blast/*_5utr.tsv") into blast5utrChannelresume
    //file("blast/*_3d.tsv") into blast3dChannelresume
    //file("blast/*_vp1.tsv") into vp1blastChannelresume
    //file("blast/*_p1.tsv") into p1blastChannelresume
    set contigsID, file("*_nt.tsv") into blastChannel
    file("*.tsv") into allblastChannel mode flatten
    file("*_p1.tsv") into p1bisblastChannel
    file("*_vp1.tsv") into vp1bisblastChannel
    file("*_5utr.tsv") into blast5utrbisChannel
    file("*_3d.tsv") into blast3dbisChannel

    shell:
    """
    #!/usr/bin/env bash
    blastn -query !{contigs}  -db !{params.vp1} -num_threads !{params.cpus} \
           -out !{contigsID}_vp1.tsv -max_target_seqs 1 \
           -outfmt '6 qseqid sseqid qlen length qstart qend sstart send pident qcovs evalue'\
           -task blastn -evalue !{params.evalue} -max_hsps 1
    blastn -query !{contigs}  -db !{params.seq5utr} -num_threads !{params.cpus} \
           -out !{contigsID}_5utr.tsv -max_target_seqs 1 \
           -outfmt '6 qseqid sseqid qlen length qstart qend sstart send pident qcovs evalue'\
           -task blastn -evalue !{params.evalue} -max_hsps 1
    blastn -query !{contigs}  -db !{params.seq3d} -num_threads !{params.cpus} \
           -out !{contigsID}_3d.tsv -max_target_seqs 1 \
           -outfmt '6 qseqid sseqid qlen length qstart qend sstart send pident qcovs evalue'\
           -task blastn -evalue !{params.evalue} -max_hsps 1
    blastn -query !{contigs} -db !{params.p1} -num_threads !{params.cpus} \
           -out !{contigsID}_p1.tsv -max_target_seqs 1 \
           -outfmt '6 qseqid sseqid qlen length qstart qend sstart send pident qcovs evalue'\
           -task blastn -evalue !{params.evalue} -max_hsps 1
    blastn -query !{contigs}  -db !{params.nt}  -num_threads !{params.cpus} \
           -out !{contigsID}_nt.tsv -evalue !{params.evalue}\
           -max_target_seqs !{params.numberBestannotation} \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    """
}

/*p1blastChannelresume = Channel.create()
vp1blastChannelresume = Channel.create()
blast5utrChannelresume = Channel.create()
blast3dChannelresume = Channel.create()*/
allblastChannel.subscribe { it.copyTo(blastDir) }
p1bisblastChannel.collectFile(name: 'p1blast.tsv')
                 .set{ p1blastChannelresume }
vp1bisblastChannel.collectFile(name: 'vp1blast.tsv')
                  .set{ vp1blastChannelresume }
blast5utrbisChannel.collectFile(name: '5utrblast.tsv')
                  .set{ blast5utrChannelresume }
blast3dbisChannel.collectFile(name: '3dblast.tsv')
                 .set { blast3dChannelresume }



// Extract vp1 sequences
process vp1 {
    publishDir "$myDir", mode: 'copy'

    input:
    set contigsID, file(contigs), file(vp1blast) from vp1blastChannel

    output:
    file("vp1/*_vp1.fasta") into vp1Channel
    set contigsID,file("vp1/*_vp1.fasta") into vp1fastaChannel
    file("vp1_contigs/*_vp1_contigs.fasta") into vp1contigsChannel
    file("vp1_contigs/*_vp1_contigs.fasta") into vp1contigsChannelresume
    file("vp1/*_vp1_list.txt") into vp1setChannel

    shell:
    """
    #!/usr/bin/env bash
    mkdir vp1 vp1_contigs
    nlines=\$(wc -l !{vp1blast} |cut -f 1 -d ' ')
    if [ "\${nlines}" -gt '0' ]
    then
        # Extract VP1 sequence
        extract_sequence.py -q !{contigs} -b !{vp1blast} \
               -o vp1/!{contigsID}_vp1.fasta -a !{params.vp1info} \
               -c !{params.vp1coverage} -t vp1
        # Extract VP1 contig
        grab_catalogue_sequence.py -i !{vp1blast} \
               -d !{contigs} -o vp1_contigs/!{contigsID}_vp1_contigs.fasta \
               -a !{params.vp1info} -v !{params.vp1coverage}
        nseq=\$(grep "^>" -c vp1/!{contigsID}_vp1.fasta)
        if [ "\${nseq}" -gt '0' ]
        then
            grep "^>" vp1/!{contigsID}_vp1.fasta |cut -f 1 -d " " | sed "s:>::g" |sed "s:_vp1::g" > vp1/!{contigsID}_vp1_list.txt
        else
            touch vp1/!{contigsID}_vp1_list.txt
        fi
    else
        touch vp1_contigs/!{contigsID}_vp1_contigs.fasta vp1/!{contigsID}_vp1.fasta vp1/!{contigsID}_vp1_list.txt
    fi
    """
}
vp1listChannel =  vp1setChannel.collectFile(name: 'vp1_list.txt')

vp1Channel.collectFile(name: 'vp1sequences.fasta')
          .into { vp1tosave ; vp1fasta}
vp1tosave.subscribe { it.copyTo(myDir) }

process blast_vp1 {
    publishDir "$myDir/blast/", mode: 'copy'
    memory "5G"
    cpus params.cpus

    input:
    set contigsID, file(vp1fasta) from vp1fastaChannel

    output:
    set contigsID, file("*_nt_vp1.tsv") into blastvp1ntChannel

    shell:
    """
    #!/usr/bin/env bash
    blastn -query !{vp1fasta}  -db !{params.nt}  -num_threads !{params.cpus} \
           -out !{contigsID}_nt_vp1.tsv -evalue !{params.evalue}\
           -max_target_seqs 1 -max_hsps 1 \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    """
}

process annotation_vp1 {
    label 'python3'
    //clusterOptions='--qos=normal -p common'
    publishDir "$myDir/annotation/", mode: 'copy'
    //memory "10G"

    //beforeScript ='source /local/gensoft2/adm/etc/profile.d/modules.sh;module use /pasteur/projets/policy01/Matrix/modules'
    //module = 'Python/3.6.0:taxadb/0.6.0'

    input:
    set contigsID, file(vp1_blast) from blastvp1ntChannel

    output:
    file("*_annotation_vp1.tsv") into vp1finalChannel

    shell:
    """
    #!/usr/bin/env bash
    nlines=\$(wc -l !{vp1_blast}|cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
        # Get the taxonomy
        get_taxonomy3.py -i !{vp1_blast} \
        -d !{params.taxadb} -o !{contigsID}_taxonomy_vp1.tsv
        ExtractNCBIDB2.py -f !{vp1_blast} \
        -g !{contigsID}_taxonomy_vp1.tsv -nb !{params.numberBestannotation} \
        -o !{vp1_blast.baseName}_annotation_vp1.tsv
    else
        touch !{vp1_blast.baseName}_annotation_vp1.tsv
    fi
    """
}
vp1finalChannelresume = vp1finalChannel.collectFile(name:"vp1combined.tsv")


process other_region {
    publishDir "$myDir", mode: 'copy'
    cache 'deep'

    input:
    set contigsID, file(contigs), file(p1blast), file(blast5utr), file(blast3d)  from otherblastChannel
    each file(vp1list) from vp1listChannel

    output:
    file("p1/*_p1.fasta") into p1Channel
    file("5utr/*_5utr.fasta") into fasta5utrChannel
    file("3d/*_3d.fasta") into fasta3dChannel
    set contigsID, file("p1/*_p1.fasta"), file("5utr/*_5utr.fasta"), file("3d/*_3d.fasta") into otherfastaChannel

    shell:
    """
    #!/usr/bin/env bash
    mkdir p1 5utr 3d
    nlines=\$(wc -l !{p1blast} |cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
       # Extract VP1 sequence
       extract_sequence.py -q !{contigs} -b !{p1blast} \
            -o p1/!{contigsID}_p1.fasta -a !{params.p1info} \
            -c !{params.p1coverage} -t p1 -l !{vp1list}
    else
        touch p1/!{contigsID}_p1.fasta
    fi
    nlines=\$(wc -l !{blast5utr} |cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
       # Extract 5UTR sequence
       extract_sequence.py -q !{contigs} -b !{blast5utr} \
            -o 5utr/!{contigsID}_5utr.fasta -a !{params.info5utr} \
            -c !{params.coverage5utr} -t 5utr -l !{vp1list}
    else
        touch 5utr/!{contigsID}_5utr.fasta
   fi
   nlines=\$(wc -l !{blast3d} |cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
       # Extract 3D sequence
       extract_sequence.py -q !{contigs} -b !{blast3d} \
            -o 3d/!{contigsID}_3d.fasta -a !{params.info3d} \
            -c !{params.coverage3d} -t 3d -l !{vp1list}
    else
        touch 3d/!{contigsID}_3d.fasta
    fi
   """
}

p1Channel.collectFile(name: 'p1sequences.fasta')
         .into{ p1mbma ; p1tosave ; p1fasta }
p1tosave.subscribe { it.copyTo(myDir) }

fasta5utrChannel.collectFile(name: '5utrsequences.fasta')
                .into { fna5utr ; tosave5utr }
tosave5utr.subscribe { it.copyTo(myDir) }

fasta3dChannel.collectFile(name:'3dsequences.fasta')
              .into { fna3d ; tosave3d }
tosave3d.subscribe { it.copyTo(myDir) }
process blast_other {
    publishDir "$myDir/blast/", mode: 'copy'
    memory "5G"
    cpus params.cpus
    cache 'deep'

    input:
    set contigsID, file(p1blast), file(blast5utr), file(blast3d) from otherfastaChannel

    output:
    //file("blast/*_nt_p1.tsv"), file("blast/*_nt_5utr.tsv"), file("blast/*_nt_3d.tsv") into otherblastChannelresume
    //set file("blast/*_nt_p1.tsv"), file("blast/*_nt_5utr.tsv"), file("blast/*_nt_3d.tsv") into otherblastChannelresume
    set contigsID, file("*_nt_p1.tsv"), file("*_nt_5utr.tsv"), file("*_nt_3d.tsv") into otherblastntChannel
    //file("blast/*.tsv") into otherblastncbiChannel

    shell:
    """
    #!/usr/bin/env bash
    blastn -query !{p1blast}  -db !{params.nt}  -num_threads !{params.cpus} \
           -out !{contigsID}_nt_p1.tsv -evalue !{params.evalue}\
           -max_target_seqs 1 -max_hsps 1 \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    blastn -query !{blast5utr}  -db !{params.nt}  -num_threads !{params.cpus} \
           -out !{contigsID}_nt_5utr.tsv -evalue !{params.evalue}\
           -max_target_seqs 1 -max_hsps 1 \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
    blastn -query !{blast3d}  -db !{params.nt}  -num_threads !{params.cpus} \
           -out !{contigsID}_nt_3d.tsv -evalue !{params.evalue}\
           -max_target_seqs 1 -max_hsps 1 \
           -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'
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
        extract_sequence_vsearch.py -i !{vp1vsearch} \
               -o vp1/!{contigsID}_vp1_vsearch.fasta -a !{params.vp1info}
        # Extract VP1 contig
        grab_catalogue_sequence.py -i !{vp1vsearch} \
               -d !{contigs} -o vp1_contigs/!{contigsID}_vp1_contigs_vsearch.fasta
    else
        touch vp1_contigs/!{contigsID}_vp1_contigs_vsearch.fasta vp1/!{contigsID}_vp1_vsearch.fasta
    fi
    """
}*/

// Extract ncbi tree
process annotation {
    label 'python3'
    //clusterOptions='--qos=normal -p common'
    publishDir "$myDir/annotation/", mode: 'copy'
    //memory "10G"

    //beforeScript ='source /local/gensoft2/adm/etc/profile.d/modules.sh;module use /pasteur/projets/policy01/Matrix/modules'
    //module = 'Python/3.6.0:taxadb/0.6.0'

    input:
    set contigsID, file(ncbi_annotation) from blastChannel

    output:
    file("*_annotation.tsv") into annotationChannel mode flatten
    file("*_taxonomy.tsv") into taxChannel mode flatten

    shell:
    """
    #!/usr/bin/env bash
    nlines=\$(wc -l !{ncbi_annotation}|cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
        # Get the taxonomy
        get_taxonomy3.py -i !{ncbi_annotation} \
        -d !{params.taxadb} -o !{contigsID}_taxonomy.tsv
        ExtractNCBIDB2.py -f !{ncbi_annotation} \
        -g !{contigsID}_taxonomy.tsv -nb !{params.numberBestannotation} \
        -o !{ncbi_annotation.baseName}_annotation.tsv
    else
        touch !{ncbi_annotation.baseName}_annotation.tsv
        touch !{ncbi_annotation.baseName}_taxonomy.tsv
    fi
    """
}

process annotation_other {
    label 'python3'
    publishDir "$myDir/annotation/", mode: 'copy'
    //beforeScript ='source /local/gensoft2/adm/etc/profile.d/modules.sh;module use /pasteur/projets/policy01/Matrix/modules'
    //module = 'Python/3.6.0:taxadb/0.6.0'

    input:
    set contigsID, file(blast_nt_p1), file(blast_nt_5utr), file(blast_nt_3d) from otherblastntChannel

    output:
    file("*_p1_annotation.tsv") into p1finalChannel
    file("*_5utr_annotation.tsv") into final5utrChannel
    file("*_3d_annotation.tsv") into final3dChannel

    shell:
    """
    #!/usr/bin/env bash
    nlines=\$(wc -l < !{blast_nt_p1})
    if [ \${nlines} -gt '0' ]
    then
        # Get the taxonomy
        get_taxonomy3.py -i !{blast_nt_p1} \
        -d !{params.taxadb} -o !{contigsID}_taxonomy_p1.tsv
        ExtractNCBIDB2.py -f !{blast_nt_p1} \
        -g !{contigsID}_taxonomy_p1.tsv -nb !{params.numberBestannotation} \
        -o !{blast_nt_p1.baseName}_annotation.tsv
    else
        touch !{blast_nt_p1.baseName}_annotation.tsv
    fi
    nlines=\$(wc -l < !{blast_nt_5utr})
    if [ \${nlines} -gt '0' ]
    then
        # Get the taxonomy
        get_taxonomy3.py -i !{blast_nt_5utr} \
        -d !{params.taxadb} -o !{contigsID}_taxonomy_5utr.tsv
        ExtractNCBIDB2.py -f !{blast_nt_5utr} \
        -g !{contigsID}_taxonomy_5utr.tsv -nb !{params.numberBestannotation} \
        -o !{blast_nt_5utr.baseName}_annotation.tsv
    else
        touch !{blast_nt_5utr.baseName}_annotation.tsv
    fi
    nlines=\$(wc -l < !{blast_nt_3d})
    if [ \${nlines} -gt '0' ]
    then
        # Get the taxonomy
        get_taxonomy3.py -i !{blast_nt_3d} \
        -d !{params.taxadb} -o !{contigsID}_taxonomy_3d.tsv
        ExtractNCBIDB2.py -f !{blast_nt_3d} \
        -g !{contigsID}_taxonomy_3d.tsv -nb !{params.numberBestannotation} \
        -o !{blast_nt_3d.baseName}_annotation.tsv
    else
        touch !{blast_nt_3d.baseName}_annotation.tsv
    fi
    """

}

p1finalChannelresume = p1finalChannel.collectFile(name: 'p1combined.tsv')
final5utrChannelresume = final5utrChannel.collectFile(name: 'combined5utr.tsv')
final3dChannelresume = final3dChannel.collectFile(name: 'combined3d.tsv')

combineChannel = annotationChannel.collectFile(name: 'combined.tsv')
//vp1contigs_vsearch = vp1contigsvChannel.collectFile(name: 'vp1contigs_vsearch.fasta')
vp1contigs = vp1contigsChannel.collectFile(name: 'vp1contigs.fasta')
//vp1contigs.into { vp1contigs_mbma; vp1contigs_s }
//vp1annotation = vp1blasttocombineChannel.collectFile(name: 'vp1blast.tsv')

process buildIndex {
    cpus params.cpus
    cache 'deep'
    //publishDir "$myDir", mode: 'copy'
    //clusterOptions='--qos=normal -p common'

    input:
    file p1mbma

    output:
    file("p1mbma.index*") into p1mbma_index
    //file(vp1contigs) into vp1contigs_fasta

    when:
    params.abundance == "yes"

    script:
    """
    bowtie2-build --threads ${params.cpus} ${p1mbma} p1mbma.index
    """
}

process abundance_vp1 {
    //queue = 'common'
    //clusterOptions='--qos=normal -p common'
    publishDir "$myDir", mode: 'copy'
    cache 'deep'
    //conda '/pasteur/homes/aghozlan/condaenvs/entvirus_env'
    //errorStrategy 'finish'

    input:
    file(p1mbma) from p1mbma_index.first()
    file(forward) from r1Channel.toList()
    file(reverse) from r2Channel.toList()
    //file cleanDir

    when:
    params.abundance == "yes"

    output:
    file("abundance/count_matrix.tsv") into countChannel
    file("abundance/count_matrix.tsv") into counttofilerChannel
    // change to shared
    shell:
    """
    #!/usr/bin/env bash
    mbma.py mapping --r1 !{forward} --r2 !{reverse} -o abundance \
           -db p1mbma.index -e !{params.mail} -q !{params.queue} \
           -p !{params.partition} --bowtie2 \
           --best -m PE -t !{params.cpus}
    mv abundance/comptage/count_matrix.txt abundance/count_matrix.tsv
    """
}

process filter_count {
    //queue = 'common'
    //clusterOptions='--qos=normal -p common'
    publishDir "$myDir/abundance/", mode: 'copy'
    cache 'deep'
    //errorStrategy 'finish'

    input:
    file(count) from counttofilerChannel
    

    when:
    params.abundance == "yes"

    output:
    file("abundance/count_matrix_shaman.tsv") into count_shamanChannel
    // change to shared
    shell:
    """
    #!/usr/bin/env bash
    filter_count_matrix.py !{count} !{params.filter_matrix} count_matrix_shaman.tsv
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
    cache 'deep'

    input:
    file(ncbi_annotation) from combineChannel
    file(vp1contigs_fasta) from vp1contigs

    output:
    file("contigs_annotation.tsv") into finalChannel
    file("vp1contigs_annotation.tsv") into vp1contigsfinalChannel
    file("vp1contigs_annotation.tsv") into vp1contigsfinalChannelabundance
    file(vp1contigs_fasta) into vp1contigsfasta

    shell:
    """
    #!/usr/bin/env bash
    nlines=\$(wc -l !{ncbi_annotation}|cut -f 1 -d ' ')
    if [ \${nlines} -gt '0' ]
    then
        # Combine the annotation
        combine_annotation.py -i !{ncbi_annotation} \
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
                strange.py -i contigs_annotation.tsv \
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
//file("vp1contigs_annotation.tsv") into vp1finalChannelabundance
//p1contigs_annotation.subscribe { it.copyTo(myDir) }
if ( params.abundance == "no") {
    process summary {
        publishDir "$myDir", mode: 'copy'
        cache 'deep'

        input:
        file(rawr1) from r1Channelresume.toList()
        file(rawr2) from r2Channelresume.toList()
        file(procr1) from procr1Channelresume.toList()
        file(procr2) from procr2Channelresume.toList()
        file(contigs) from contigsChannelresume.toList()
        file(vp1contigs) from vp1contigsChannelresume.toList()
        file(vp1blast) from vp1blastChannelresume //.toList()
        file(p1blast) from p1blastChannelresume //.toList()
        file(blast3d) from blast3dChannelresume //.toList()
        file(blast5utr) from blast5utrChannelresume //.toList()
        file(seqvp1_annot) from vp1finalChannelresume
        file(seqp1_annot) from p1finalChannelresume
        file(seq5utr_annot) from final5utrChannelresume
        file(seq3d_annot) from final3dChannelresume
        file(vp1contigs_annot) from vp1contigsfinalChannel

        output:
        file("result_summary.tsv") into resultChannel
        file("occurence/*.tsv") into annotabundChannel mode flatten
        //file(vp1contigs_annot) into overChannel

        shell:
        """
        mkdir occurence
        extract_result2.py -r1 !{rawr1} -r2 !{rawr2} -pr1 !{procr1} -pr2 !{procr2} -a !{params.vp1info}\
               -vp1 !{vp1contigs_annot} -o result_summary.tsv \
               -l !{params.annotated} -v !{params.vp1coverage}\
               -p !{params.p1info} -c !{contigs} -vp1c !{vp1contigs} \
               -3d !{params.info3d} -5utr !{params.info5utr} \
               -bvp1 !{vp1blast} -bp1 !{p1blast} \
               -b5utr !{blast5utr} -b3d !{blast3d} \
               -bntvp1 !{seqvp1_annot} -bntp1 !{seqp1_annot} \
               -bnt5utr !{seq5utr_annot} -bnt3d !{seq3d_annot} \
               -s !{params.ent_serotype} -oc occurence/count.tsv -oa occurence/annotation.tsv
        """
    }
}
else{
    process summary_abundance {
        publishDir "$myDir", mode: 'copy'
        cache 'deep'

        input:
        file(rawr1) from r1Channelresume.toList()
        file(rawr2) from r2Channelresume.toList()
        file(procr1) from procr1Channelresume.toList()
        file(procr2) from procr2Channelresume.toList()
        file(contigs) from contigsChannelresume.toList()
        file(vp1contigs) from vp1contigsChannelresume.toList()
        file(vp1blast) from vp1blastChannelresume
        file(p1blast) from p1blastChannelresume
        file(blast3d) from blast3dChannelresume
        file(blast5utr) from blast5utrChannelresume
        file(seqvp1_annot) from vp1finalChannelresume
        file(seqp1_annot) from p1finalChannelresume
        file(seq5utr_annot) from final5utrChannelresume
        file(seq3d_annot) from final3dChannelresume
        file(vp1contigs_annot) from vp1contigsfinalChannel
        file(countmatrix) from countChannel

        output:
        file("result_summary.tsv") into resultChannel
        file("abundance/annotation.tsv") into annotabundChannel
        file("occurence/*.tsv") into occurenceChannel mode flatten
        //file(vp1contigs_annot) into overChannel

        shell:
        """
        mkdir abundance occurence
        extract_result2.py -r1 !{rawr1} -r2 !{rawr2} -pr1 !{procr1} -pr2 !{procr2} -a !{params.vp1info}\
            -vp1 !{vp1contigs_annot} -o result_summary.tsv \
            -n !{countmatrix} -l !{params.annotated} -v !{params.vp1coverage}\
            -p !{params.p1info} -c !{contigs} -vp1c !{vp1contigs}\
            -3d !{params.info3d} -5utr !{params.info5utr} \
            -bvp1 !{vp1blast} -bp1 !{p1blast}\
            -b5utr !{blast5utr} -b3d !{blast3d} \
            -bntvp1 !{seqvp1_annot} -bntp1 !{seqp1_annot} \
            -bnt5utr !{seq5utr_annot} -bnt3d !{seq3d_annot} \
            -s !{params.ent_serotype} -os abundance/annotation.tsv -oc occurence/count.tsv -oa occurence/annotation.tsv
        """
    }
}

//resultChannel.subscribe { it.copyTo(myDir) }
//resultChannelabundance.subscribe { it.copyTo(myDir) }
//fastachan = Channel.from(["Contigs_with_VP1", ${vp1contigsfasta}], ["P1_sequences", ${p1fasta}], ["VP1_sequences", ${vp1fasta}])

process enterovirus_classification {
    publishDir "$myDir", mode: 'copy'
    cache 'deep'

    input:
    file(result_summary) from resultChannel
    file vp1contigsfasta
    file p1fasta
    file vp1fasta
    file fna5utr
    file fna3d
    //set key, fasta from fastachan

    output:
    //set "!{fasta.baseName}", file("*.fasta") into fastaserotype
    //file("enterovirus_classified/*.fasta") into fastaclassified mode flatten
    //file("enterovirus_classified/*.fasta") into classified mode flatten
    file("enterovirus_classified/*_{contigs,p1,vp1,5utr,3d}.fasta") optional true into fastaclassified mode flatten
    set file("enterovirus_classified/*_{contigs,p1,vp1,5utr,3d}_rooted.fasta"), file("enterovirus_classified/*_{contigs,p1,vp1,5utr,3d}_association.txt") optional true into fastaclassified_rooted mode flatten
    file("phylogeny/*.txt") optional true into itolChannel mode flatten

    shell:
    """
    mkdir enterovirus_classified/ phylogeny/
    enterovirus_classification.py -i !{result_summary} -f !{vp1contigsfasta} -o enterovirus_classified/ -t "Contigs_with_VP1" -itol phylogeny/ -c
    enterovirus_classification.py -i !{result_summary} -f !{p1fasta} -o enterovirus_classified/ -t "P1_sequences" -itol phylogeny/ -c
    enterovirus_classification.py -i !{result_summary} -f !{vp1fasta} -o enterovirus_classified/ -t "VP1_sequences" -itol phylogeny/ -c
    enterovirus_classification.py -i !{result_summary} -f !{fna5utr} -o enterovirus_classified/ -t "5UTR_sequences" -itol phylogeny/ -c
    enterovirus_classification.py -i !{result_summary} -f !{fna3d} -o enterovirus_classified/ -t "3D_sequences" -itol phylogeny/ -c
    # ROOTING trees with an other class
    enterovirus_classification.py -i !{result_summary} -f !{vp1contigsfasta} -o enterovirus_classified/ -t "Contigs_with_VP1" -r !{params.seqfull} -itol phylogeny/ -e !{params.ent_serotype} -c
    enterovirus_classification.py -i !{result_summary} -f !{p1fasta} -o enterovirus_classified/ -t "P1_sequences" -r !{params.p1} -itol phylogeny/ -e !{params.ent_serotype} -c
    enterovirus_classification.py -i !{result_summary} -f !{vp1fasta} -o  enterovirus_classified/ -t "VP1_sequences" -r !{params.vp1} -itol phylogeny/ -e !{params.ent_serotype} -c
    enterovirus_classification.py -i !{result_summary} -f !{fna5utr} -o enterovirus_classified/ -t "5UTR_sequences" -r !{params.seq5utr} -itol phylogeny/ -e !{params.ent_serotype} -c
    enterovirus_classification.py -i !{result_summary} -f !{fna3d} -o enterovirus_classified/ -t "3D_sequences" -r !{params.seq3d} -itol phylogeny/ -e !{params.ent_serotype} -c
    """
}

process multiple_alignment {
    publishDir "$myDir/msa/", mode: 'copy'
    cpus params.cpus
    cache 'deep'

    input:
    //set fastaID, file(fasta) from fastaserotype
    file(fasta) from fastaclassified

    output:
    //set fastaID, file("msa/*.ali") into msaserotype
    file("*.ali") optional true into msadata //mode flatten

    shell:
    """
    #!/usr/bin/env bash
    nseq=\$( grep -c "^>" !{fasta} )
    if [ "\${nseq}" -ge '200' ]
    then
        mafft --retree 1 --thread !{params.cpus} !{fasta} > !{fasta.baseName}.ali
    elif [ "\${nseq}" -gt '1' ]
    then
        mafft  --thread !{params.cpus}  --maxiterate 1000 --localpair !{fasta} > !{fasta.baseName}.ali
    fi
    """
}

process multiple_alignment_rooted {
    publishDir "$myDir/msa/", mode: 'copy', pattern: 'msa/*.ali'
    cpus params.cpus_phylogeny
    cache 'deep'

    input:
    //set fastaID, file(fasta) from fastaserotype
    set file(fasta), file(association) from fastaclassified_rooted

    output:
    //set fastaID, file("msa/*.ali") into msaserotype
    set file("*.ali"), file(association) optional true into msadata_rooted mode flatten

    shell:
    """
    #!/usr/bin/env bash
    nseq=\$( grep -c "^>" !{fasta} )
    if [ "\${nseq}" -ge '200' ]
    then
        mafft --retree 1 --thread !{params.cpus} !{fasta} > !{fasta.baseName}.ali
    elif [ "\${nseq}" -gt '3' ]
    then
        mafft  --thread !{params.cpus}  --maxiterate 1000 --localpair !{fasta} > !{fasta.baseName}.ali
    fi
    """
}

process filtering_alignment {
    publishDir "$myDir/msa/", mode: 'copy'
    cache 'deep'

    input:
    //set fastaID, file(msa) from msaserotype
    file(msa) from msadata

    output:
    //set fastaID, file("*_bmge.ali") into msafiltserotype
    file("*_bmge.ali") optional true into msadatafilt //mode flatten
    file("*_bmge_large.ali") optional true into msalargedatafilt

    shell:
    """
    #!/usr/bin/env bash
    nseq=\$( grep -c "^>" !{msa} )
    if [ "\${nseq}" -ge '200' ]
    then
        BMGE -i ${msa} -t DNA -m ID -h 1 -g !{params.conserved_position} -w 1 -b 1 -of !{msa.baseName}_bmge_large.ali
    elif [ "\${nseq}" -gt '3' ]
    then
        BMGE -i ${msa} -t DNA -m ID -h 1 -g !{params.conserved_position} -w 1 -b 1 -of !{msa.baseName}_bmge.ali
    fi
    """
}

process filtering_alignment_rooted {
    publishDir "$myDir/msa/", mode: 'copy', pattern: 'msa/*.ali'
    cache 'deep'

    input:
    //set fastaID, file(msa) from msaserotype
    set file(msa), file(association) from msadata_rooted

    output:
    //set fastaID, file("*_bmge.ali") into msafiltserotype
    set file("*_bmge.ali"), file(association) optional true into msadatafilt_rooted mode flatten
    set file("*_bmge_large.ali"), file(association) optional true into msalargedatafilt_rooted mode flatten

    shell:
    """
    #!/usr/bin/env bash
    mkdir -p msa
    nseq=\$( grep -c "^>" !{msa} )
    if [ "\${nseq}" -ge '200' ]
    then
        BMGE -i ${msa} -t DNA -m ID -h 1 -g !{params.conserved_position} -w 1 -b 1 -of !{msa.baseName}_bmge_large.ali
    elif [ "\${nseq}" -gt '3' ]
    then
        BMGE -i ${msa} -t DNA -m ID -h 1 -g !{params.conserved_position} -w 1 -b 1 -of !{msa.baseName}_bmge.ali
    fi
    """
}

process phylogeny {
    publishDir "$myDir/phylogeny/", mode: 'copy'
    cache 'deep'

    input:
    //set fastaID, file(msafilt) from msafiltserotype
    file(msafilt) from msadatafilt

    output:
    file("*.treefile") into phylogenyserotype //mode flatten

    shell:
    """
    #!/usr/bin/env bash
    iqtree -m GTR+I+G4 -s !{msafilt}
    """
}

process phylogeny_large {
    publishDir "$myDir/phylogeny/", mode: 'copy'
    cpus params.cpus_phylogeny
    cache 'deep'

    input:
    //set fastaID, file(msafilt) from msafiltserotype
    file(msafilt) from msalargedatafilt

    output:
    file("*.treefile") into phylogenyserotype_large //mode flatten

    shell:
    """
    #!/usr/bin/env bash
    iqtree -m GTR+I+G4 -nt !{params.cpus_phylogeny} -s !{msafilt}
    """
}

process phylogeny_rooted {
    publishDir "$myDir/phylogeny/", mode: 'copy'
    cache 'deep'

    input:
    //set fastaID, file(msafilt) from msafiltserotype
    set file(msafilt), file(association) from msadatafilt_rooted

    output:
    file("*_final.treefile") into phylogenyserotype_rooted //mode flatten

    shell:
    """
    #!/usr/bin/env bash
    iqtree -m GTR+I+G4 -s !{msafilt}
    name=\$(echo !{msafilt.baseName} | cut -f 1,2 -d "_")
    gotree reroot outgroup -r -i !{msafilt.baseName}.ali.treefile  -l \${name}_association.txt > !{msafilt.baseName}_final.treefile
    """
}

process phylogeny_rooted_large {
    publishDir "$myDir/phylogeny/", mode: 'copy'
    cpus params.cpus_phylogeny
    cache 'deep'

    input:
    //set fastaID, file(msafilt) from msafiltserotype
    set file(msafilt), file(association) from msalargedatafilt_rooted

    output:
    file("*.treefile") into phylogenyserotype_large_rooted //mode flatten

    shell:
    """
    #!/usr/bin/env bash
    iqtree -m GTR+I+G4 -nt !{params.cpus_phylogeny} -s !{msafilt}
    name=\$(echo !{msafilt.baseName} | cut -f 1,2 -d "_")
    gotree reroot outgroup -r -i !{msafilt.baseName}.ali.treefile  -l \${name}_association.txt > !{msafilt.baseName}_final.treefile
    """
}

println "Project : $workflow.projectDir"
//println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
