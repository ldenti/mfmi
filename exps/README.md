# Experiments

```
git submodule update --init --recursive

cd ext/pp-rlcsa
cd rlcsa
make -j2
cd ..
make -j2

cd ../pp-ropebwt2/ropebwt2
make -j2
cd ..
make -j2

cd ../..

# Edit config/config.yaml
snakemake -c16 -p [-n]
```