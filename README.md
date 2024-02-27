# MFMI

```
make
./mfmi index -r example/4.fa.gz > INDEX
./mfmi search INDEX example/reads.fq.gz
```

- [X] code refactoring
- [X] rank on interval
- [X] iterate from begin/end of vector
- [ ] tests
- [ ] bit rope and rle (7+1)
- [X] rlcsa merging
- [X] "iterative" rlcsa construction
- [X] more fastas
- [X] dump/load rlcsa
- [ ] parallelize index construction
- [ ] improve array sorting
- [ ] improve rope merging
- [ ] fmd-index
