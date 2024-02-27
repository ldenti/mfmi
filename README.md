# MFMI

```
make
./mfmi index -r example/4.fa.gz > INDEX
./mfmi search INDEX example/perfect_reads.fq.gz
./mfmi pingpong INDEX example/reads.fq.gz > specific_strings.sfs
```

- [X] code refactoring
- [X] rank on interval
- [X] iterate from begin/end of vector
- [X] ~tests~ example (+ simple comparison with [panpp](github.com/ldenti/panpp))
- [X] rlcsa merging
- [X] "iterative" rlcsa construction
- [X] more fastas
- [X] dump/load rlcsa
- [X] fmd-index

### Index construction - ideas
We can keep the 4 bit vectors and:
- [ ] parallelize index construction
- [ ] implement "single-bit" rope and rle (7+1 for first block..)
- [ ] improve array sorting
- [ ] improve rope merging
or we can try to check if we can use:
- [ ] single 3bit rope (`$ACGTN``)+ update with marked positions
