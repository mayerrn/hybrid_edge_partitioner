/*
 * Vertex2EdgePart.cpp
 *
 *  Created on: May 14, 2020
 *      Author: root
 */

#include "Vertex2EdgePart.hpp"

Vertex2EdgePart::Vertex2EdgePart(std::string basefilename) {
	this->basefilename = basefilename;
	p = FLAGS_p;
}

Vertex2EdgePart::~Vertex2EdgePart() {
	// TODO Auto-generated destructor stub
}

void Vertex2EdgePart::readVertexPartitioning(){
	// open the partitioning file
	// assumption: basefilename + ".part"
	std::string partfilename = basefilename + ".part." + std::to_string(p); //METIS file name format

	 std::string line;
	 std::ifstream partfile(partfilename);
	 uint32_t current_vertex = 0;
	 int result;
	 while(std::getline(partfile, line)){
		 result = std::stoi(line);
		 current_vertex++; // as in METIS format, first vertex id is 1, not 0!
		 vertex2partition[current_vertex] = result;
	 }
}

int Vertex2EdgePart::findEdgePartition(vid_t u, vid_t v){
	int target_u = vertex2partition[u];
	int target_v = vertex2partition[v];
	if (u == v){
		return target_u;
	}
	else{
		// flip a coin: edge goes to either u's or v's partition
		if (rand() % 2 == 0){
			return target_u;
		}
		else{
			return target_v;
		}
	}
}


void Vertex2EdgePart::split(){
    // read the metis adjacency list
    FILE *inf = fopen(basefilename.c_str(), "r");
    if (inf == NULL) {
            LOG(FATAL) << "Could not load:" << basefilename
                       << " error: " << strerror(errno) << std::endl;
        }
        LOG(INFO) << "Reading in adjacency list format!" << std::endl;

        int maxlen = 1000000000;
        char *s = (char *)malloc(maxlen);

        size_t bytesread = 0;

        char delims[] = " \t";
        size_t linenum = 0;

        while (fgets(s, maxlen, inf) != NULL) {
        	linenum++;

            FIXLINE(s);
            bytesread += strlen(s);

            if (s[0] == '#')
                continue; // Comment
            if (s[0] == '%')
                continue; // Comment


            if (linenum == 1){
            	char *t = strtok(s, delims);
            	if (t == NULL)
            		LOG(FATAL) << "First line must contain num verts and num edges" << std::endl;; // empty line

            	num_vertices = atoi(t);
                t = strtok(NULL, delims);
                if (t != NULL){
                	num_edges = atol(t);
                }
                else {
                	LOG(FATAL) << "First line must contain num verts and num edges" << std::endl;
                }
                LOG(INFO) << "Vertices: " << num_vertices << ", Edges: " << num_edges << std::endl;

                initDataStructures();
                readVertexPartitioning();

                continue; //done with first line
            }

            vid_t from = linenum - 1; // because first line contained the num of verts and edges
            char *t = strtok(s, delims);
            if (t == NULL)
            	continue;

           do {
        	   vid_t to = atoi(t);
        	   if (from != to && from < to) { //ignore one direction, because METIS format contains both directions of an undirected edge
        		   int target_p = findEdgePartition(from, to);
        		   // "assign" edge to target_p
        		   counter[target_p]++;
        		   is_mirrors[target_p].set_bit_unsync(from);
        		   is_mirrors[target_p].set_bit_unsync(to);
        	   }
           }
           while((t = strtok(NULL, delims)) != NULL);


        }
    free(s);
    fclose(inf);

	size_t max_occupied = *std::max_element(counter.begin(), counter.end());
	LOG(INFO) << "balance: " << (double)max_occupied / ((double)num_edges / p);
	size_t total_mirrors = 0;
	rep (i, p)
		total_mirrors += is_mirrors[i].popcount();
	LOG(INFO) << "total mirrors: " << total_mirrors;
	LOG(INFO) << "replication factor: " << (double)total_mirrors / num_vertices;
}
