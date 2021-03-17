#pragma once
#include <set>

#include <vector>
#include <parallel/algorithm>

#include <iostream>
#include <fstream>
#include <unordered_map>
#include "util.hpp"
#include "dense_bitset.hpp"

/*
 * Entry in the index array. Contains a reference to the column array and the lengths (in and out) of the entry.
 * Differentiating in and out length allows for preserving direction in a directed graph.
 */
class mem_adjlist_t{
	public:
	vid_t *adj; // link into the column array
    vid_t len_out;
    vid_t len_in;

  public:

    mem_adjlist_t() :   adj(NULL), len_out(0), len_in(0) {}
    mem_adjlist_t( vid_t *adj) :  adj(adj), len_out(0), len_in(0) {}

    vid_t *begin() { return adj; }
    vid_t *end() { return adj + len_out + len_in; }

    size_t size() const {return len_out + len_in;}
    size_t size_out() const {return len_out;}
    size_t size_in() const {return len_in;}

    void push_back_out(vid_t data){
    	adj[len_out++] = data;
    };

    void push_back_in(vid_t data, vid_t offset){
    	adj[len_in++ + offset] = data; // incoming neighbors are written after the outgoing neighbors.
     };


    void erase_out(vid_t to_erase_pos){
    	adj[to_erase_pos] = adj[size_out()-1];
    	adj[size_out()-1] = adj[size()-1];
    	pop_back_out();
    }

    void erase_in(vid_t to_erase_pos){
        adj[to_erase_pos] = adj[size()-1];
        pop_back_in();
    }

    void pop_back_out(){len_out--;}
    void pop_back_in(){len_in--;}
};

/*
 * in-memory graph
 */
class mem_graph_t{

public:
  vid_t num_vertices;
  size_t nedges;
  vid_t *neighbors;
  std::vector<mem_adjlist_t> vdata;
  vid_t high_degree_threshold; // starting from which degree is a node considered a high-degree node (if exceeded)
  double high_degree_factor; // average degree * hdf = hdt
  std::fstream h2h_file; // file that keeps edges between two high-degree vertices on external memory
  std::fstream low_degree_file; // file that keeps edges incident to a low-degree vertex on external memory
  size_t num_h2h_edges;


public:
  mem_graph_t() : num_vertices(0), nedges(0), high_degree_factor(0), high_degree_threshold(0), num_h2h_edges(0) {  }

  ~mem_graph_t()
  {
        if (neighbors){
            free(neighbors);
        }
  }

  void resize(vid_t _num_vertices)
  {
      num_vertices = _num_vertices;
      vdata.resize(num_vertices);
  }

  size_t num_edges() const { return nedges; }

  size_t stream_build(std::ifstream &fin, size_t num_edges, dense_bitset &is_high_degree, dense_bitset &has_high_degree_neighbor, std::vector<size_t> &count, bool write_low_degree_edgelist);
  mem_adjlist_t &operator[](size_t idx) { return vdata[idx]; };

  const mem_adjlist_t &operator[](size_t idx) const {
		return vdata[idx];
	}


};

