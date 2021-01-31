
# Overview
rnabridge-denova implements an efficient algorithm to bridge paired-end RNA-seq reads, i.e.,
to find the entire fragments for the given paired-end reads.
See [rnabridge-test](https://github.com/Shao-Group/rnabridge-test) for the evaluation of this tool.

# Installation

## Install Bifrost
rnabridge-denova uses additional library Bifrost for de Bruijn graph construction, the instrcution for downloading and installing Bifrost is here(https://github.com/pmelsted/bifrost).

If the install path for Bifrost is not the default path. Remember to set the environment variables C_INCLUDE_PATH, CPLUS_INCLUDE_PATH, LD_LIBRARY_PATH, LIBRARY_PATH and PATH correctly.

Assuming the header files (.h) are located at the path /usr/local/include/, the following command set the environment variables C_INCLUDE_PATH and CPLUS_INCLUDE_PATH correctly for the time of the session:

```
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include/
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/include/
```
executing the binary of Bifrost fails because libbifrost.so or libbifrost.a is not found

Assuming that libbifrost.(so|dylib|a) is located at the path /usr/local/lib/, the following command set the environment variables LD_LIBRARY_PATH, LIBRARY_PATH and PATH correctly for the time of the session:


```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib/
export PATH=$PATH:/usr/local/lib/
```

Assuming that the binary executable file Bifrost is located at the path /usr/local/bin/, the following command set the environment variable PATH correctly for the time of the session:
```
export PATH=$PATH:/usr/local/bin/
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
Bifrost build -t 4 -k 31 -c -i -d -s /pathtoread1.fq -s /pathtoread1.fq -o ./dbg
./src/graph_align update -t 4 -g ./dbg.gfa -f ./dbg.bfg_colors -q ./graph -x ./query
./src/find_bridge ./graph ./query ./bridge
./src/bridge2fa ./bridge ./bridge.fa
```


