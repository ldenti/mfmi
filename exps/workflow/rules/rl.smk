# NOTE: index creates partial indexes based on input fastas without any check (so multiple concurrent runs on same fastas will break something for sure)
rule rl_index:
    input:
        get_fasta_names,
    output:
        index=pjoin(WD, "{n}", "rl-index.rlcsa.parameters"),
    params:
        index_prefix=pjoin(WD, "{n}", "rl-index"),
    log:
        log=pjoin(WD, "{n}", "logs", "rl_index.log"),
        time=pjoin(WD, "{n}", "times", "rl_index.time"),
    threads: THREADS
    shell:
        """
        /usr/bin/time -vo {log.time} ./ext/pp-rlcsa/rl index -@ {threads} -i {params.index_prefix} {input} 2> {log.log}
        """


rule rl_exact:
    input:
        index=pjoin(WD, "{n}", "rl-index.rlcsa.parameters"),
        fq=pjoin(WD, "{n}", "perfect_reads.fq"),
    output:
        txt=pjoin(WD, "{n}", "rl-exact.out"),
    params:
        index_prefix=pjoin(WD, "{n}", "rl-index"),
    log:
        time=pjoin(WD, "{n}", "times", "rl_exact.time"),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} ./ext/pp-rlcsa/rl exact {params.index_prefix} {input.fq} > {output.txt}
        """


rule rl_pp:
    input:
        index=pjoin(WD, "{n}", "rl-index.rlcsa.parameters"),
        fa=pjoin(WD, "{n}", "query_ref.fa"),
    output:
        txt=pjoin(WD, "{n}", "rl-specifics.txt"),
    params:
        index_prefix=pjoin(WD, "{n}", "rl-index"),
    log:
        time=pjoin(WD, "{n}", "times", "rl_pp.time"),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} ./ext/pp-rlcsa/rl pingpong {params.index_prefix} {input.fa} > {output.txt}
        """
