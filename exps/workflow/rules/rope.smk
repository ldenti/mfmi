rule rope_index:
    input:
        mfa=pjoin(WD, "{n}", "reference.fa"),
    output:
        index=pjoin(WD, "{n}", "rope-index.{MT}.rld"),
    params:
        t=lambda wildcards: "" if wildcards.MT == "st" else "-@",
    log:
        log=pjoin(WD, "{n}", "logs", "rope_index.{MT}.log"),
        time=pjoin(WD, "{n}", "times", "rope_index.{MT}.time"),
    threads: lambda wildcards: 1 if wildcards.MT == "st" else 4
    shell:
        """
        /usr/bin/time -vo {log.time} ./ext/pp-ropebwt2/rt2 index {params.t} {input.mfa} > {output.index}
        """


rule rope_compare_indexes:
    input:
        index=pjoin(WD, "{n}", "rope-index.st.rld"),
        index_mt=pjoin(WD, "{n}", "rope-index.mt.rld"),
    output:
        txt=pjoin(WD, "{n}", "rope-index.diff"),
    shell:
        """
        diff {input.index} {input.index_mt} > {output.txt}
        """


rule rope_exact:
    input:
        index=pjoin(WD, "{n}", "rope-index.st.rld"),
        fq=pjoin(WD, "{n}", "perfect_reads.fq"),
    output:
        txt=pjoin(WD, "{n}", "rope-exact.out"),
    log:
        time=pjoin(WD, "{n}", "times", "rope_exact.time"),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} ./ext/pp-ropebwt2/rt2 exact {input.index} {input.fq} > {output.txt}
        """


rule rope_pp:
    input:
        index=pjoin(WD, "{n}", "rope-index.st.rld"),
        fa=pjoin(WD, "{n}", "query_ref.fa"),
    output:
        txt=pjoin(WD, "{n}", "rope-specifics.txt"),
    log:
        time=pjoin(WD, "{n}", "times", "rope_pp.time"),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} ./ext/pp-ropebwt2/rt2 pingpong {input.index} {input.fa} > {output.txt}
        """
