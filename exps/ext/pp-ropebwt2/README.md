# PingPong using ropebwt2 index

```
git submodule update --init --recursive ropebwt2
cd ropebwt2
make -j2
cd ..
make -j2

./rt2 index ../../../example/4.fa.gz > index.fmd
./rt2 pingpong index.fmd ../../../example/reads.fq.gz > specific_strings.sfs
```
