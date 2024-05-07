rule prl_index:
    input:
        mfa=pjoin(WD, "{n}", "reference.fa"),
    output:
        index=pjoin(WD, "{n}", "prl-index.rld"),
    log:
        log=pjoin(WD, "{n}", "logs", "prl_index.log"),
        time=pjoin(WD, "{n}", "times", "prl_index.time"),
    threads: THREADS
    shell:
        """
        /usr/bin/time -vo {log.time} ../mfmi index {input.mfa} > {output.index} 2> {log.log}
        """


rule prl_exact:
    input:
        index=pjoin(WD, "{n}", "prl-index.rld"),
        fq=pjoin(WD, "{n}", "perfect_reads.fq"),
    output:
        txt=pjoin(WD, "{n}", "prl-exact.out"),
    log:
        time=pjoin(WD, "{n}", "times", "prl_exact.time"),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} ../mfmi exact {input.index} {input.fq} > {output.txt}
        """


rule prl_pp:
    input:
        index=pjoin(WD, "{n}", "prl-index.rld"),
        fa=pjoin(WD, "{n}", "query_ref.fa"),
    output:
        txt=pjoin(WD, "{n}", "prl-specifics.txt"),
    log:
        time=pjoin(WD, "{n}", "times", "prl_pp.time"),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} ../mfmi pingpong {input.index} {input.fa} > {output.txt}
        """
