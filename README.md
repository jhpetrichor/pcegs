# PCEGS
Identification of protein complexes throuth integrative gene ontology parent-child semantic similarity and dynamic protein-protein interactions

## usage
### step 0: loggor

install spglod<https://github.com/gabime/spdlog>
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