# PCEGS
Identification of protein complexes throuth integrative gene ontology parent-child semantic similarity and dynamic protein-protein interactions

## usage
### step 0: logging

install spglod <https://github.com/gabime/spdlog>
```bash
$ git clone https://github.com/gabime/spdlog.git
$ cd spdlog && mkdir build && cd build
$ cmake .. && cmake --build .
```

### step 1: complie
create build directory and compile

```bash
$ mkdir build && cd build && cmake ..
$ make
```
### step 2: run pcegs
pcegs data_name result alpha bete

example
```bash
$ cd bin
$ ./pcegs data_name result 0.5 0.5
```

## illustration
### data
Gene expression, gene theory and essential proteins are placed in the folder dataset/Yeast/DAG.

PPIs are placed in the dataset/Yeast/PPI
### input
By default, PPI data in Yeast is used. Place it under the dataset/Yeast/PPI folder. Parameter ignore file suffix.
### output
The prediction results are stored in the results directory by default. Each line represents a protein complex, and the proteins are separated by tabs.
<!-- ### evaluation -->