Hybrid Edge Partitioner
=============================================

Implementation of HEP, as published in SIGMOD 2021.

Please cite the paper as follows:

> Ruben Mayer and Hans-Arno Jacobsen. 2021. Hybrid Edge Partitioner: Partitioning Large Power-Law Graphs under Memory Constraints. In Proceedings of the 2021 International Conference on Management of Data (SIGMOD ’21), June 20–25, 2021, Virtual Event, China. ACM, New York, NY, USA, 14 pages. https://doi.org/10.1145/3448016.3457300

Compilation and Usage
---------------------

We tested our program on Ubuntu 18.04, and it requires the following
packages: `cmake`, `glog`, `gflags`, `boost`:
```
sudo apt-get install libgoogle-glog-dev libgflags-dev libboost-all-dev
```

Compilation:
```
git clone <this repository>
cd <repository_name>
mkdir release && cd release
cmake ..
make -j8
```

Usage:
```
$ ./main --help
main: -filename <path to the input graph> [-method hep] [-hdf <threshold / \tau>] [-p <number of partitions>] 
```

**Example.** Partition the Orkut graph into 8 parts using HEP with \tau = 10.0:
```
$ ./main -p 8 -method hep -hdf 10.0 -filename /path/to/com-orkut.ungraph.txt
```

Acknowledgements
---------------------
Parts of the implementation are based on the NE reference implementation provided by Qin Liu: https://github.com/ansrlab/edgepart

This refers to the following classes. Some of them were adapted to fit the different data formats in HEP/NE++ compared to NE.
conversions.cpp/hpp
dense_bitset.hpp
edgepart.hpp
graph.cpp/hpp 
min_heap.hpp
util.cpp/hpp

We also integrated (for comparison to NE++) the original implementation of NE:
ne_graph.cpp/hpp
ne_min_heap.hpp
ne_partitioner.cpp/hpp
