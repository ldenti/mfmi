# MFMI

```
make
./mfmi index example/4.fa.gz > INDEX
./mfmi search INDEX example/perfect_reads.fq.gz
./mfmi pingpong INDEX example/reads.fq.gz > specific_strings.sfs
# diff specific_strings.sfs example/4.sfs
```
