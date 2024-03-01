# PingPong using ropebwt2 index

```
git submodule update --init --recursive rlcsa
cd rlcsa
make -j2
cd ..
make -j2

./rl index -r -i index ../../../example/4.fa.gz
./rl pingpong index ../../../example/reads.fq.gz > specific_strings.sfs
```