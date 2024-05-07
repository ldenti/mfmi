rule get_references:
    input:
        agc=AGC,
    output:
        rlist=pjoin(WD, "{n}", "reference.list"),
    params:
        n="{n}",
        seed=SEED,
    shell:
        """
        agc listset {input.agc} | (RANDOM={params.seed}; while read line ; do echo "$RANDOM $line" ; done ) | sort | head -n {params.n} | cut -f2 -d' ' > {output.rlist}
        """


rule get_fasta:
    input:
        agc=AGC,
        rlist=pjoin(WD, "{n}", "reference.list"),
    output:
        fa=pjoin(WD, "{n}", "first_entry.fa"),
        mfa=pjoin(WD, "{n}", "reference.fa"),
    shell:
        """
        while read line
        do
            agc getset {input.agc} $line
        done < {input.rlist} > {output.mfa}
        samtools faidx {output.mfa}

        agc getset {input.agc} $(head -1 {input.rlist}) > {output.fa}
        samtools faidx {output.fa}
        """


def get_fasta_names(wildcards):
    ck_output = checkpoints.split_fasta.get(**wildcards).output.mfa_dir
    (SMP,) = glob_wildcards(pjoin(ck_output, "{sample}.fa"))
    return expand(pjoin(ck_output, "{SAMPLE}.fa"), SAMPLE=SMP)


checkpoint split_fasta:
    input:
        agc=AGC,
        rlist=pjoin(WD, "{n}", "reference.list"),
    output:
        mfa_dir=directory(pjoin(WD, "{n}", "fastas")),
    shell:
        """
        mkdir -p {output.mfa_dir}
        while read line
        do
            agc getset {input.agc} $line > {output.mfa_dir}/$line.fa
            samtools faidx {output.mfa_dir}/$line.fa
        done < {input.rlist}
        """


rule simulate_reads:
    input:
        fa=pjoin(WD, "{n}", "first_entry.fa"),
    output:
        fq=pjoin(WD, "{n}", "perfect_reads.fq"),
    params:
        seed=SEED,
        n=1000000,
        l=100,
    shell:
        """
        wgsim -e 0 -N {params.n} -1 {params.l} -2 {params.l} -r 0 -R 0 -X 0 -S {params.seed} {input.fa} {output.fq} /dev/null
        """


rule get_query_ref:
    input:
        agc=AGC,
    output:
        fa=pjoin(WD, "{n}", "query_ref.fa"),
    params:
        seed=SEED,
    shell:
        """
        idx=$(agc listset {input.agc} | (RANDOM={params.seed}; while read line ; do echo "$RANDOM $line" ; done ) | sort | tail -1 | cut -f2 -d' ')
        agc getset {input.agc} $idx > {output.fa}
        samtools faidx {output.fa}
        """
