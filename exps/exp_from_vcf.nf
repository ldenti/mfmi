#!/usr/bin/env nextflow

params.fa = "/home/ld/data/mfmi/test.fa"
params.vcf = "/home/ld/data/mfmi/test.vcf.gz"

params.rlcsa_bin = "/home/ld/code/panpp/PANPP"
params.svdss_bin = "/home/ld/code/SVDSS_linux_x86-64"

process get_samples {
    input:
    val n

    output:
    path('samples.list')

    """
    bash $baseDir/scripts/get_samples.sh $params.vcf $n > samples.list
    """
}

process get_rlist {
    input:
    path(list)

    output:
    tuple val('r'), path('r.list')

    """
    head -n -1 $list > r.list
    """
}

process get_qlist {
    input:
    path(list)

    output:
    tuple val('q'), path('q.list')

    """
    tail -1 $list > q.list
    """
}

process get_fasta {
    input:
    tuple val(t), path(list)

    output:
    tuple val(t), path("${t}.fa")

    shell:
    '''
    while read idx
    do
        bcftools consensus --haplotype 1 --samples $idx --fasta-ref !{params.fa} !{params.vcf} | sed "s/>.*$/>$idx/g" 
    done < !{list} > !{t}.fa
    wc -l !{t}.fa
    '''
}

process mfmi_index {
    input:
    path(fa)
    
    output:
    path("index.mfmi")

    """
    $baseDir/../mfmi index -r -@1 $fa > index.mfmi
    """
}

process mfmi_query {
    input:
    tuple path(index), path(fa)

    output:
    path("mfmi.sfs")
    
    """
    $baseDir/../mfmi pingpong $index $fa > mfmi.sfs
    """
}

process rlcsa_index {
    input:
    path(fa)
    
    output:
    path("INDEX.rlcsa.array")

    """
    ${params.rlcsa_bin} index -r -@1 -i INDEX $fa
    """
}

process rlcsa_query {
    input:
    tuple path(index), path(fa)

    output:
    path("rlcsa.sfs")
    
    """
    ${params.rlcsa_bin} search INDEX $fa -a > rlcsa.sfs
    """
}

process svdss_index {
    input:
    path(fa)
    
    output:
    path("index.fmd")

    """
    ${params.svdss_bin} index --fasta $fa --index index.fmd
    """
}

process svdss_query {
    input:
    tuple path(index), path(fa)

    output:
    path("svdss.sfs")
    
    """
    ${params.svdss_bin} search --index $index --fastq $fa
    cat *.sfs > svdss.sfs
    """
}

/*
 * Define the workflow
 */
workflow {
    ( get_samples(1) | ( get_rlist & get_qlist ) | mix | get_fasta ).branch {
        ref: it[0] == 'r'
        query: it[0] == 'q'
    }.set { input }

    input.ref.view { "$it is ref" }
    input.query.view { "$it is query" }
    ref = input.ref.map{ it[1] }
    input.query.map { it[1] }.set{ query }
    ref.view { "ref: $it" }
    query.view { "query: $it" }

    mfmi_index(ref).concat(query) | collect | mfmi_query | view
    svdss_index(ref).concat(query) | collect | svdss_query | view
    rlcsa_index(ref).concat(query) | collect | rlcsa_query | view

}
