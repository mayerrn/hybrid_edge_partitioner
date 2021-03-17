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

#include "util.hpp"
#include "ne_min_heap.hpp"
#include "dense_bitset.hpp"
#include "edgepart.hpp"
#include "partitioner.hpp"
#include "ne_graph.hpp"

/* Neighbor Expansion (NE) */
class NePartitioner : public Partitioner
{
  private:
    const double BALANCE_RATIO = 1.00;

    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges, assigned_edges;
    int p, bucket;
    double average_degree;
    size_t capacity;

    std::vector<edge_t> edges;
    ne_graph_t adj_out, adj_in;
    NeMinHeap<vid_t, vid_t> min_heap;
    std::vector<size_t> occupied;
    std::vector<vid_t> degrees;
    std::vector<int8_t> master;
    std::vector<dense_bitset> is_cores, is_boundarys;

    size_t min_size = 0; // currently smallest partition
    size_t max_size = 0; // currently largest partition
    double lambda;

    bool write_out_partitions = false; // whether the partitions should be written to the out-file or not
    bool hybrid_partitioning = false; // whether to apply HEP-style hybrid partitioning using NE
    bool hybrid_random_stream = false; // in case of hybrid_partitioning, whether to do random streaming instead of HDRF

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    edgepart_writer<vid_t, uint16_t> writer;

    int check_edge(const edge_t *e)
    {
        rep (i, bucket) {
            auto &is_boundary = is_boundarys[i];
            if (is_boundary.get(e->first) && is_boundary.get(e->second) &&
                occupied[i] < capacity) {
                return i;
            }
        }

        rep (i, bucket) {
            auto &is_core = is_cores[i], &is_boundary = is_boundarys[i];
            if ((is_core.get(e->first) || is_core.get(e->second)) &&
                occupied[i] < capacity) {
                if (is_core.get(e->first) && degrees[e->second] > average_degree)
                    continue;
                if (is_core.get(e->second) && degrees[e->first] > average_degree)
                    continue;
                is_boundary.set_bit(e->first);
                is_boundary.set_bit(e->second);
                return i;
            }
        }

        return p;
    }

    void assign_edge(int bucket, vid_t from, vid_t to)
    {
    	if (write_out_partitions){
    		writer.write_edge_assignment(from, to, bucket);
    	}
        assigned_edges++;
        occupied[bucket]++;
        degrees[from]--;
        degrees[to]--;
    }

    void add_boundary(vid_t vid)
    {
        auto &is_core = is_cores[bucket], &is_boundary = is_boundarys[bucket];

        if (is_boundary.get(vid))
            return;
        is_boundary.set_bit_unsync(vid);

        if (!is_core.get(vid)) {
            min_heap.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
        }

        rep (direction, 2) {
            ne_adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
            for (size_t i = 0; i < neighbors.size();) {
                if (edges[neighbors[i].v].valid()) {
                    vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                    if (is_core.get(u)) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid);
                        min_heap.decrease_key(vid);
                        edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else if (is_boundary.get(u) &&
                               occupied[bucket] < capacity) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid);
                        min_heap.decrease_key(vid);
                        min_heap.decrease_key(u);
                        edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else
                        i++;
                } else {
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                }
            }
        }
    }

    void occupy_vertex(vid_t vid, vid_t d)
    {
        CHECK(!is_cores[bucket].get(vid)) << "add " << vid << " to core again";
        is_cores[bucket].set_bit_unsync(vid);

        if (d == 0)
            return;

        add_boundary(vid);

        for (auto &i : adj_out[vid])
            if (edges[i.v].valid())
                add_boundary(edges[i.v].second);
        adj_out[vid].clear();

        for (auto &i : adj_in[vid])
            if (edges[i.v].valid())
                add_boundary(edges[i.v].first);
        adj_in[vid].clear();
    }

    bool get_free_vertex(vid_t &vid)
    {
        vid = dis(gen);
        vid_t count = 0;
        while (count < num_vertices &&
               (adj_out[vid].size() + adj_in[vid].size() == 0 ||
                adj_out[vid].size() + adj_in[vid].size() >
                    2 * average_degree ||
                is_cores[bucket].get(vid))) {
            vid = (vid + ++count) % num_vertices;
        }
        if (count == num_vertices)
            return false;
        return true;
    }

    void assign_remaining();
    double compute_partition_score(vid_t u, vid_t v, int bucket_id); // returns HDRF score for edge (u,v) on partition <bucket_id>
    int best_scored_partition(vid_t u, vid_t v); // returns bucket id where score is best for edge (u,v)

    double compute_partition_score_tiebreaking_balance(vid_t u, vid_t v, int bucket_id); // returns HDRF score, but using balancing only as a  tie breaker (no balancing score)
    int best_scored_partition_tiebreaking_balance(vid_t u, vid_t v); // returns bucket id where score is best for edge (u,v), using the tiebreaking adaptation of HDRF scoring

//    void assign_master(); // do not assign master, as HEP does not either
    size_t count_mirrors();

    size_t random_streaming(); // returns number of streamed edges from input h2h file
    size_t hdrf_streaming(); // returns number of streamed edges from input h2h file

  public:
    NePartitioner(std::string basefilename);
    void split();
};
