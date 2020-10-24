# Stochastic Molecular Dynamics!

[![Build Status](https://travis-ci.org/richardtjornhammar/stochMD.svg?branch=master)](https://travis-ci.org/richardtjornhammar/stochMD)

This is a molecular dynamics package supporting grand canonical simulations of
systems with bonded and non-bonded interactions. The velocity re-scaling for the
temperature coupling is stochastic.

# Build 

To build the program either issue the command:

```
gcc src/*.c -lm -o stochMD
```
or:

```
cmake .
make
```

# Run a grand canonical argon simulation :
To run the program issue commands like so (MD and MC):
```
./stochMD -i example/argon.inp -c example/single.xyz -p example/single.top -o argon.nfo MC
```

# License

GPL-3
