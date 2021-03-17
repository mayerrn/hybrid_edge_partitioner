Hybrid Edge Partitioner
=============================================


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

