# PingPong using rlcsa index

```
git submodule update --init --recursive rlcsa
cd rlcsa
make -j2
cd ..
make -j2

./rl index -i INDEX ../../../example/4.fa.gz
./rl exact INDEX ../../../example/perfect_reads.fq.gz
./rl pingpong INDEX ../../../example/reads.fq.gz > specific_strings.sfs
```
