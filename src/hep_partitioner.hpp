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

#include <set>

#include "util.hpp"
#include "min_heap.hpp"
#include "dense_bitset.hpp"
#include "edgepart.hpp"
#include "partitioner.hpp"
#include "graph.hpp"

/* Hybrid Edge Partitioner (HEP) */
class HepPartitioner : public Partitioner
{
  private:
    const double BALANCE_RATIO = 1.0;
    double lambda;
    bool extended_metrics;

    bool two_ps; // whether to use 2ps restreaming for the second phase of HEP
    bool stream_random; // whether to use random streaming for the second phase of HEP

    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges, assigned_edges, num_h2h_edges;
    int p, bucket; 
    double average_degree;
    size_t capacity;
    size_t capacity_in_memory; // capacity per partition of in memory partitioning
    size_t invalidated_edges_count; // number of edges removed in clean up phase overall
    size_t min_size = 0; // currently smallest partition
    size_t max_size = 0; // currently largest partition

    bool write_out_partitions = false; // whether the partitions should be written to the out-file or not
    bool write_low_degree_edgelist = false; // whether edges incident to a low-degree vertex should be written out to a file. useful if this sub-graph should be analyzed separately.

    std::vector<edge_t> edges;
    mem_graph_t mem_graph; // graph for in-memory processing
    double high_degree_factor;
    MinHeap<vid_t, vid_t> min_heap;
    std::vector<size_t> occupied;
    std::vector<dense_bitset> is_boundarys; 
    dense_bitset is_in_a_core;
    dense_bitset is_high_degree;
    dense_bitset has_high_degree_neighbor;
    std::vector<size_t> count; // degrees of vertices//(num_vertices, 0);


    vid_t search_index_free_vertex = 0;


    edgepart_writer<vid_t, uint16_t> writer;
    

    void in_memory_clean_up_neighbors(vid_t vid, dense_bitset & is_core, dense_bitset & is_boundary);


    void assign_edge(int bucket, vid_t from, vid_t to)
    {        
    	if (write_out_partitions){
    		writer.write_edge_assignment(from, to, bucket);
    	}
        assigned_edges++;
        occupied[bucket]++;

        is_boundarys[bucket].set_bit_unsync(from);
        is_boundarys[bucket].set_bit_unsync(to);

    }

    void in_memory_add_boundary(vid_t vid){

    	bool bucket_full = false;
    	auto &is_core = is_in_a_core, &is_boundary = is_boundarys[bucket];

		if (is_boundary.get(vid)){
			return; //vid is already in boundary!
		}

		is_boundary.set_bit_unsync(vid);

		if (is_high_degree.get(vid)){ // high degree vertices are treated as if they were in the core
			is_in_a_core.set_bit_unsync(vid);
			return; // high degree vertices are ignored in the normal expansion. we do not at all look at their neighbors. we do also not add them to min heap.
		}

		bool vid_is_in_core = is_core.get(vid);

		if (!vid_is_in_core) {
			min_heap.insert(mem_graph[vid].size(), vid); //is not in the core yet: Potential next candidate --> insert in MinHeap
		}
		auto &neighbors = mem_graph[vid].adj;
		vid_t count = 0;
		for(; count < mem_graph[vid].size_out(); count++) //for the adj_out neighbors
		{
			if (occupied[bucket] >= capacity){
				bucket_full = true;
			} // full, stop adding vertices to the boundary of this bucket

			vid_t &u = neighbors[count];

			if (is_high_degree.get(u)){ // high degree vertices are always considered to be in c
				if (!bucket_full){
					assign_edge(bucket, vid , u ); //assign edge --> vid is the left vertex
					if (!vid_is_in_core) {
						min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
					}
				}
				else{ //bucket is full; assign to next bucket
					assign_edge(bucket + 1, vid , u );
				}

			}
			else{
				if(is_core.get(u)){ //If the neighbor of vid is in core
					if (!bucket_full){
						assign_edge(bucket, vid , u ); //assign edge --> vid is the left vertex
						if (!vid_is_in_core) {
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
						}
					}
					else {
						//bucket is full; assign to next bucket
						assign_edge(bucket + 1, vid , u );
					}
				}else if (is_boundary.get(u)) {
					if (!bucket_full){
						assign_edge(bucket, vid , u );
						min_heap.decrease_key(u, 1, mem_graph[u].size());
						if (!vid_is_in_core) {
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
						}
					}
					else {
						//bucket is full; assign to next bucket
						assign_edge(bucket + 1, vid , u );
					}
				}
			}

		}
		for(; count < mem_graph[vid].size(); count++) //for the adj_in neighbors
		{

			if (occupied[bucket] >= capacity){
				bucket_full = true;
			} // full, stop adding vertices to the boundary

			vid_t &u = neighbors[count];

			if (is_high_degree.get(u)){ // high degree vertices are always considered to be in c
				if (!bucket_full){
					assign_edge(bucket, u , vid ); //assign edge --> vid is the right vertex
					if (!vid_is_in_core) {
						min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
					}
				}
				else {
					//bucket is full; assign to next bucket
					assign_edge(bucket + 1, u , vid );
				}
			}
			else{
				if(is_core.get(u)){
					if (!bucket_full){
						assign_edge(bucket, u, vid  ); //vid is on the right side
						if (!vid_is_in_core) {
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
						}
					}
					else {
						//bucket is full; assign to next bucket
						assign_edge(bucket + 1, u , vid );
					}
				}else if (is_boundary.get(u)) {
					if (!bucket_full){
						assign_edge(bucket, u, vid  );
						min_heap.decrease_key(u, 1, mem_graph[u].size());
						if (!vid_is_in_core) {
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
						}
					}
					else {
						//bucket is full; assign to next bucket
						assign_edge(bucket + 1, u , vid );
					}
				}
			}
		}
    }


    void in_memory_occupy_vertex(vid_t vid, vid_t d){

    	CHECK(!is_in_a_core.get(vid)) << "add " << vid << " to core again";

    	is_in_a_core.set_bit_unsync(vid);
    	if (d == 0){ //also applies to high degree vertices, as their in-memory degree is 0
    	    return;
    	}
    	in_memory_add_boundary(vid);
    	for(vid_t i = 0; i < mem_graph[vid].size(); i++) { //Set all neighbors of vid to boundary
    		in_memory_add_boundary(mem_graph[vid].adj[i]);
    	}
    }

    bool in_memory_get_free_vertex(vid_t &vid){

       vid = search_index_free_vertex;

      	/*
       	* find a vertex to start expansion with
       	*/
       	while ((mem_graph[vid].size() == 0 || is_in_a_core.get(vid)) && vid < num_vertices) {
       	   	vid++;
       	}

       	search_index_free_vertex = vid;
       	if (vid == num_vertices){
       		return false;
       	} // searched to the end, did not find free vertex
       	else{
           	return true;
       	}
    }


    void load_in_memory(std::string basefilename, std::ifstream &fin);
    void partition_in_memory();
    void in_memory_assign_remaining();

    double compute_partition_score(vid_t u, vid_t v, int bucket_id); // returns HDRF score for edge (u,v) on partition <bucket_id>
    int best_scored_partition(vid_t u, vid_t v); // returns bucket id where score is best for edge (u,v)

    size_t count_mirrors();
    void compute_stats();

    void random_streaming();
    void hdrf_streaming();


  public:
    HepPartitioner(std::string basefilename);
    void split();
    //needed for testing
    std::vector<dense_bitset> get_boundary(){
        return is_boundarys;
    }

	vid_t getNumVertices() const {
		return num_vertices;
	}

	void setNumVertices(vid_t numVertices) {
		num_vertices = numVertices;
	}
};
