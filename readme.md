
# Overview
rnabridge-denova implements an efficient algorithm to bridge paired-end RNA-seq reads, i.e.,
to find the entire fragments for the given paired-end reads.
See [rnabridge-test](https://github.com/Shao-Group/rnabridge-test) for the evaluation of this tool.

# Installation

## Install Bifrost
rnabridge-denova uses additional library Bifrost for de Bruijn graph construction, the instrcution for downloading and installing Bifrost is here(https://github.com/pmelsted/bifrost).

If you add add the option `-DCMAKE_INSTALL_PREFIX=/pathtobifrost` in to the `cmake` command in building Bifrost, the install path for Bifrost will not be the default path. Remember to set the environment variables C_INCLUDE_PATH, CPLUS_INCLUDE_PATH, LD_LIBRARY_PATH, LIBRARY_PATH and PATH correctly.


```
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/pathtobifrost/include/
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/pathtobifrostinclude/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pathtobifrost/lib/
export LIBRARY_PATH=$LIBRARY_PATH:/pathtobifrost/lib/
export PATH=$PATH:/pathtobifrost/lib/:/pathtobifrost/bin/
```

## Compile rnabridge-align

Use the following to compile rnabridge-align:
```
cd src
make
```

# Usage

The usage of `rnabridge-denova` is:

```
./rnabridge-denova pathtoread1.fq pathtoread2.fq pathtooutputbridge.fa
```

pathtoread1.fq and pathtoread2.fq are your input paired-end RNA-seq reads. pathtooutputbridge.fa is the path for the output file, which contains all the brdiges in fasta format.
