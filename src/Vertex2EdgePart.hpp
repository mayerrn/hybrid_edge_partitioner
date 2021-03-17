#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <parallel/algorithm>

#include "util.hpp"
#include "dense_bitset.hpp"
#include "edgepart.hpp"
#include "partitioner.hpp"

class Vertex2EdgePart : public Partitioner {
private:
	const size_t BUFFER_SIZE = 64 * 1024 / sizeof(edge_t);
	std::string basefilename;

    vid_t num_vertices;
    size_t num_edges;
    int p;

    int fin;
    off_t filesize;

    std::vector<int> vertex2partition; // maps vertex id to partition id
	std::vector<dense_bitset> is_mirrors;
    std::vector<size_t> counter;

    // Removes \n from the end of line
    void FIXLINE(char *s)
    {
        int len = (int)strlen(s) - 1;
        if (s[len] == '\n')
            s[len] = 0;
    }

    void initDataStructures()
    {
    	vertex2partition.resize(num_vertices + 1, -1);
    	counter.resize(p, 0);
    	is_mirrors.resize(p, dense_bitset(num_vertices + 1));
    }

public:
	Vertex2EdgePart(std::string basefilename);
	virtual ~Vertex2EdgePart();

	void split();

	void readVertexPartitioning();

	int findEdgePartition(vid_t u, vid_t v);
};
