#include "conversions.hpp"
#include "util.hpp"
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include "hep_partitioner.hpp"

HepPartitioner::HepPartitioner(std::string basefilename) // @suppress("Class members should be properly initialized")
    : basefilename(basefilename), writer(basefilename)
{

    Timer convert_timer;
    convert_timer.start();

    Converter *converter = new Converter(basefilename);
    convert(basefilename, converter); //converts the original edgelist in a suitable format for this application
    delete converter;

    convert_timer.stop();
    LOG(INFO) << "convert time: " << convert_timer.get_time(); 
    total_time.start();
    LOG(INFO) << "initializing partitioner";
    std::ifstream fin(binedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate); //are the edges in the format from convert(basefilename, converter)
    auto filesize = fin.tellg();
    LOG(INFO) << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);
    //Before the actual edges are coming in the binedgelist-file, there is an entry of the number of vertices and an entry of number of edges
    fin.read((char *)&num_vertices, sizeof(num_vertices)); 
    fin.read((char *)&num_edges, sizeof(num_edges));
    LOG(INFO) << "num_vertices: " << num_vertices
              << ", num_edges: " << num_edges;
    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t), filesize);


    count.resize(num_vertices, 0);

    p = FLAGS_p; //number of partitions. Value comes from console
    lambda = FLAGS_lambda; //for weighing in balancing score in streaming
    extended_metrics = FLAGS_extended_metrics; // displaying extended metrics
    write_out_partitions = FLAGS_write_results; // writing out partitions to file
    write_low_degree_edgelist = FLAGS_write_low_degree_edgelist; // writing low degree edgelist to file
    stream_random = FLAGS_random_streaming; // random streaming
    average_degree = (double)num_edges * 2 / num_vertices;
    assigned_edges = 0; //will take track of how many edges are assigned to a bucket so far
    capacity = (double)num_edges * BALANCE_RATIO / p + 1; //will be used to as stopping criterion later
    occupied.assign(p, 0);  //Will count how many edges are in one partition

    is_boundarys.assign(p, dense_bitset(num_vertices)); //shows if a vertex is in S of a bucket
    is_in_a_core = dense_bitset(num_vertices); //Shows if a vertex is in ANY C
    is_high_degree = dense_bitset(num_vertices); // whether a vertex has a high degree and is handled differently
    has_high_degree_neighbor = dense_bitset(num_vertices); // whether the vertex has a high degree neighbor (important in assign_remaining function)
//    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading and constructing...";

    high_degree_factor = FLAGS_hdf;

    load_in_memory(basefilename, fin);
    capacity_in_memory = ((double)num_edges - num_h2h_edges) * BALANCE_RATIO / p + 1;

    read_timer.stop();
    LOG(INFO) << "time used for graph input and construction: " << read_timer.get_time();

};


void HepPartitioner::compute_stats(){ // these are the extended stats, including degree distributions etc.
	// average degree of vertices in C and in S\C
	size_t total_degree_C = 0, total_degree_S = 0;
	vid_t vertex_count_C = 0, vertex_count_S = 0;
	vid_t max_degree = 0;
	for (vid_t i = 0; i < num_vertices; i++){
		if (count[i] > max_degree){
			max_degree = count[i];
		}
		if (is_in_a_core.get(i)){
			total_degree_C += count[i];
			vertex_count_C++;
		}
		else{ // for those not in core, check whether they are in the boundary of any of the first k-1 partitions (the last partition is not built based on expansion)
			for (int j = 0; j < p-1; j++){
				if (is_boundarys[j].get(i)){
					total_degree_S += count[i];
					vertex_count_S++;
					break;
				}
			}
		}
	}
	double avg_deg_C = (double) total_degree_C / vertex_count_C;
	double avg_deg_S = (double) total_degree_S / vertex_count_S;
	double invalidation_fraction = (double) invalidated_edges_count / (num_edges * 2);
	LOG(INFO) << "normalized avg degree C " << avg_deg_C / average_degree << std::endl;
	LOG(INFO) << "normalized avg degree S " << avg_deg_S / average_degree << std::endl;
	LOG(INFO) << "fraction of edges invalidated in clean up phase " << invalidation_fraction << std::endl;
//	LOG(INFO) << "total degree " << total_degree_C + total_degree_S << " total vertex count " << vertex_count_C + vertex_count_S << std::endl;

	// computing the distribution of degree to replication factor
	std::vector<size_t> replication_factor_per_vertex_degree(max_degree + 1, 0);
	std::vector<vid_t> num_vertices_per_vertex_degree(max_degree + 1, 0);
	for (vid_t i = 0; i < num_vertices; i++){
		vid_t rep_factor = 0;
		for (dense_bitset &is_boundary : is_boundarys){
			if (is_boundary.get(i)){
				rep_factor++;
			}
		}
		vid_t degree = count[i];
		replication_factor_per_vertex_degree[degree] += rep_factor;
		num_vertices_per_vertex_degree[degree]++;
	}
	std::cout << "Format: degree <whitespace> replication_factor" << std::endl;
	size_t bucket_min_degree = 1;
	size_t bucket_max_degree = 10;
	vid_t vertices_in_bucket = 0;
	size_t replication_in_bucket = 0;
	for (vid_t k = 1; k <= max_degree; k++){

		vertices_in_bucket += num_vertices_per_vertex_degree[k];
		replication_in_bucket += replication_factor_per_vertex_degree[k];


		if (k == bucket_max_degree){ // make new bucket
			// finalize it print results
			double result = 0.0;
			if (vertices_in_bucket != 0){
				result = (double) replication_in_bucket / vertices_in_bucket;
			}
			std::cout << bucket_min_degree << " " << bucket_max_degree << " " << result << std::endl;
			std::cout << bucket_min_degree << " " << bucket_max_degree << " " << (double) vertices_in_bucket / num_vertices << std::endl;
			bucket_min_degree = bucket_max_degree + 1;
			bucket_max_degree = bucket_max_degree * 10;
			vertices_in_bucket = 0;
			replication_in_bucket = 0;
		}
	}

}

void HepPartitioner::load_in_memory(std::string basefilename, std::ifstream &fin){
	mem_graph.high_degree_factor = high_degree_factor;
	mem_graph.h2h_file.open(h2hedgelist_name(basefilename), std::ios_base::binary | std::ios_base::out ); // *.h2h_edgelist file
	if (write_low_degree_edgelist){
		mem_graph.low_degree_file.open(lowedgelist_name(basefilename), std::ios_base::binary | std::ios_base::out ); // *.low_edgelist file;
	}
	mem_graph.resize(num_vertices);
	num_h2h_edges = mem_graph.stream_build(fin, num_edges, is_high_degree, has_high_degree_neighbor, count, write_low_degree_edgelist);
	mem_graph.h2h_file.close(); //flushed
	if (write_low_degree_edgelist){
		mem_graph.low_degree_file.close(); //flushed
	}
}


void HepPartitioner::in_memory_assign_remaining(){

	LOG(INFO) << "Assigned edges before assign_remaining: " << assigned_edges << std::endl;

	repv (vid, num_vertices){

		if (!is_in_a_core.get(vid)){

			auto &neighbors = mem_graph[vid].adj;
			vid_t i = 0;
			for(; i < mem_graph[vid].size_out(); i++)
			{
				int target = best_scored_partition(vid, neighbors[i]);
				assign_edge(target, vid, neighbors[i]);
			}

			// in case the vertex has high degree neighbors, the edges from
			// those have not been assigned yet, as the hd vertices were ignored
			// in the expansion. Hence, we have to assign those explicitly.
			if (has_high_degree_neighbor.get(vid)){
				for(; i < mem_graph[vid].size(); i++) //for the adj_in neighbors
				{
					if (is_high_degree.get(neighbors[i])){
						int target = best_scored_partition(neighbors[i], vid);
						assign_edge(target, neighbors[i], vid);
					}
				}

			}
		}
	}

	LOG(INFO) << "Assigned edges before streaming: " << assigned_edges << std::endl;
	LOG(INFO) << "Assigning edges between high-degree vertices" << std::endl;

	if (stream_random){
		random_streaming();
	}
	else {
		hdrf_streaming();
	}

}



void HepPartitioner::random_streaming(){

	LOG(INFO) << "Streaming randomly." << std::endl;
	// assign the edges between two high degree vertices

	mem_graph.h2h_file.open(h2hedgelist_name(basefilename), std::ios_base::binary | std::ios_base::in );
	mem_graph.h2h_file.seekg(0, std::ios::beg);

	std::vector<edge_t> tmp_edges; // temporary buffer to read edges from file
	size_t chunk_size;
	size_t left_h2h_edges = mem_graph.num_h2h_edges;

	if (left_h2h_edges >= 100000){
		chunk_size = 100000; // batch read of so many edges
	}
	else {
		chunk_size = left_h2h_edges;
	}
	tmp_edges.resize(chunk_size);

	//	LOG(INFO) << "Chunk size is " << chunk_size << endl;


	while (left_h2h_edges > 0){ // edges to be read
		mem_graph.h2h_file.read((char *)&tmp_edges[0], sizeof(edge_t) * chunk_size);
		for (size_t i = 0; i < chunk_size; i++){

			bucket = std::rand() % p; // random bucket

			assign_edge(bucket, tmp_edges[i].first, tmp_edges[i].second);
			is_boundarys[bucket].set_bit_unsync(tmp_edges[i].first);
			is_boundarys[bucket].set_bit_unsync(tmp_edges[i].second);

		}

		left_h2h_edges -= chunk_size;
		if (left_h2h_edges < chunk_size){ // adapt chunk size for last batch read
		 chunk_size = left_h2h_edges;
		}
	}
}


void HepPartitioner::hdrf_streaming(){

	LOG(INFO) << "Streaming using HDRF algorithm." << std::endl;
	// assign the edges between two high degree vertices

	mem_graph.h2h_file.open(h2hedgelist_name(basefilename), std::ios_base::binary | std::ios_base::in );
	mem_graph.h2h_file.seekg(0, std::ios::beg);

	std::vector<edge_t> tmp_edges; // temporary buffer to read edges from file
	size_t chunk_size;
	size_t left_h2h_edges = mem_graph.num_h2h_edges;

	if (left_h2h_edges >= 100000){
		chunk_size = 100000; // batch read of so many edges
	}
	else {
		chunk_size = left_h2h_edges;
	}
	tmp_edges.resize(chunk_size);


		/*
		 * init min_size
		 */
	min_size = UINT64_MAX;
	for (int i = 0; i < p; i++){
		if (occupied[i] < min_size){
			min_size = occupied[i];
		}
	}

	while (left_h2h_edges > 0){ // edges to be read
		mem_graph.h2h_file.read((char *)&tmp_edges[0], sizeof(edge_t) * chunk_size);
		for (size_t i = 0; i < chunk_size; i++){

			bucket = best_scored_partition(tmp_edges[i].first, tmp_edges[i].second); // according to HDRF scoring

			assign_edge(bucket, tmp_edges[i].first, tmp_edges[i].second);
			is_boundarys[bucket].set_bit_unsync(tmp_edges[i].first);
			is_boundarys[bucket].set_bit_unsync(tmp_edges[i].second);

			if (occupied[bucket] > max_size){
				max_size = occupied[bucket];
			}
			if (occupied[bucket] == min_size){
				int min_sized_bucket_count = 0;
				for (int i = 0; i < p; i++){
					if (occupied[i] == min_size){
						min_sized_bucket_count++;
					}
				}
				if (min_sized_bucket_count == 1){
					min_size++;
				}
			}
		}

		left_h2h_edges -= chunk_size;
		if (left_h2h_edges < chunk_size){ // adapt chunk size for last batch read
		 chunk_size = left_h2h_edges;
		}
	}
}

size_t HepPartitioner::count_mirrors()
{
    size_t result = 0;
    rep (i, p)
    	LOG(INFO) << "Partition " << i << " : " << is_boundarys[i].popcount() << " vertices " << std::endl;

    rep (i, p)
        result += is_boundarys[i].popcount();

    return result;
}



void HepPartitioner::in_memory_clean_up_neighbors(vid_t vid, dense_bitset & is_core, dense_bitset & is_boundary){
	mem_adjlist_t &neighbors = mem_graph[vid];

	const size_t num_neigh_out = neighbors.size_out();
	const size_t num_neigh_size = neighbors.size();
	size_t i = 0;
	size_t j = 0;

	for(; j< num_neigh_out; j++){
	    //naive vid_t u = neighbors.get_neighbor(neighbors_file, i);
		vid_t u = neighbors.adj[i];
	    if (is_core.get(u)){ // neighbor u is in core, so edge is removed
	    	invalidated_edges_count++;
	    	neighbors.erase_out(i);
	    }
	    else if (is_boundary.get(u)){ // neighbor u is in boundary, so edge is removed
	        invalidated_edges_count++;
	    	neighbors.erase_out(i);
	    }
	    else if (is_high_degree.get(u)){
	    	invalidated_edges_count++;
	    	neighbors.erase_out(i);
	    }
	    else{
	        i++;
	    }
	 }

	 for(; j< num_neigh_size; j++){
		 vid_t u = neighbors.adj[i];
	     if (is_core.get(u)){ // neighbor u is in core, so edge is removed
	    	 invalidated_edges_count++;
	         neighbors.erase_in(i);
	     }
	     else if (is_boundary.get(u)){ // neighbor u is in boundary, so edge is removed
	    	 invalidated_edges_count++;
	         neighbors.erase_in(i);
	     }
		 else if (is_high_degree.get(u)){
			 invalidated_edges_count++;
			 neighbors.erase_in(i);
		 }
         else{
        	 i++;
	     }
	}
}


void HepPartitioner::partition_in_memory(){

	bool expansion_finished = false;

    for (bucket = 0; bucket < p - 1; bucket++) {
        LOG(INFO) << bucket << ", ";

        //DLOG(INFO) << "sample size: " << adj_out.num_edges();
        while (occupied[bucket] < capacity_in_memory) {
            vid_t d, vid = 0;
            if (!min_heap.get_min(d, vid)) {
                if (!in_memory_get_free_vertex(vid)) {
                    LOG(INFO) << "partition " << bucket
                               << " stop: no free vertices";
                    expansion_finished = true;
                    break;
                }

                d= mem_graph[vid].size(); // a high degree vertex will not be chosen by get free vertex. also will not be in min heap.


            } else {

                min_heap.remove(vid);

            }

            in_memory_occupy_vertex(vid, d);

        }

        /*
         * clean up the adjacency lists of vertices from the S set from the last bucket
         */
        auto &is_core = is_in_a_core, &is_boundary = is_boundarys[bucket];

        std::vector<std::pair<vid_t, vid_t>> heap = min_heap.getHeap();
        vid_t size = min_heap.getSize(); // do not iterate over heap more than "size" entries
        /*
         * vid is the vertex in S
         * u is the neighbor of vid in the currently examined edge
         */
        vid_t count = 0;
        for (std::vector<std::pair<vid_t, vid_t>>::iterator it = heap.begin(); it != heap.end(); ++it){
        	if (count >= size){
        		break; //stop here
        	}

        	count++;
        	vid_t vid = it->second;

        	in_memory_clean_up_neighbors(vid, is_core, is_boundary);

        }

        min_heap.clear();

        if (expansion_finished){
        	break;
        }

    }

    in_memory_assign_remaining();
    LOG(INFO) << "Finished partitioning" << std::endl;
    writer.fout.close();; // flushing out the unwritten edge assignments

}


double HepPartitioner::compute_partition_score(vid_t u, vid_t v, int bucket_id) {
	if (occupied[bucket_id] >= capacity){
//		cout << "partition " << bucket_id << " is full with " << occupied[bucket_id] << endl;
		return -1.0; // partition is full, do not choose it
	}
	size_t degree_u = count[u];
	size_t degree_v = count[v];
	size_t sum = degree_u + degree_v;
	double gu = 0.0, gv = 0.0;
	if (is_boundarys[bucket_id].get(u)){
		gu = degree_u;
		gu/=sum;
		gu = 1 + (1-gu);
	}
	if (is_boundarys[bucket_id].get(v)){
		 gv = degree_v;
		 gv /= sum;
		 gv = 1 + (1-gv);
	}

	double bal = (max_size - occupied[bucket_id]) / (1.0 + max_size - min_size);

	double score = gu + gv + lambda * bal;
	return score;
}

int HepPartitioner::best_scored_partition(vid_t u, vid_t v) {
	double best_score = -1.0;
	int best_partition = 0;
	for (int i = 0; i < p; i++){
		double score = compute_partition_score(u, v, i);
//		cout << "score for partition " << i << " is " << score << endl;
		if (score > best_score){
			best_score = score;
			best_partition = i;
		}
	}
	return best_partition;
}

void HepPartitioner::split()
{
    // int abort_counter = 0;
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer compute_timer;

    min_heap.reserve(num_vertices);

    LOG(INFO) << "partitioning...";
    compute_timer.start();


    partition_in_memory();

    compute_timer.stop();

    LOG(INFO) << "expected edges in each partition: " << num_edges / p;
    rep (i, p)
        LOG(INFO) << "edges in partition " << i << ": " << occupied[i];
    size_t max_occupied = *std::max_element(occupied.begin(), occupied.end());
    LOG(INFO) << "balance: " << (double)max_occupied / ((double)num_edges / p);
    size_t total_mirrors = count_mirrors();
    LOG(INFO) << "total mirrors: " << total_mirrors;
    LOG(INFO) << "replication factor: " << (double)total_mirrors / num_vertices;
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();

    CHECK_EQ(assigned_edges, num_edges);

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();

    /*
     * compute some stats about the partitioned graph (for further analysis)
     * if extended_metrics flag is set
     */
    if (extended_metrics){
    	compute_stats();
    }
}
