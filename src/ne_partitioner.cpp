#include "ne_partitioner.hpp"
#include "conversions.hpp"

NePartitioner::NePartitioner(std::string basefilename)
    : basefilename(basefilename), rd(), gen(rd()), writer(basefilename)
{
    Timer convert_timer;
    convert_timer.start();
    Converter *converter = new Converter(basefilename);
    convert(basefilename, converter);
    delete converter;
    convert_timer.stop();
    LOG(INFO) << "convert time: " << convert_timer.get_time();

    total_time.start();
    LOG(INFO) << "initializing partitioner";

    hybrid_partitioning = FLAGS_hybrid_NE;
    if (hybrid_partitioning){
    	LOG(INFO) << "Hybrid Partitioning is performed." << std::endl;
    }

    std::string filename;
    if (hybrid_partitioning){
    	filename = lowedgelist_name(basefilename);
    }
    else {
    	filename = binedgelist_name(basefilename);
    }
    std::ifstream fin(filename,
                      std::ios::binary | std::ios::ate);
    LOG(INFO) << "File name is " << filename << std::endl;
    auto filesize = fin.tellg();
    LOG(INFO) << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);

    fin.read((char *)&num_vertices, sizeof(num_vertices));
    fin.read((char *)&num_edges, sizeof(num_edges));

    LOG(INFO) << "num_vertices: " << num_vertices
              << ", num_edges: " << num_edges;
    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t), filesize);

    p = FLAGS_p;
    write_out_partitions = FLAGS_write_results; // writing out partitions to file

    lambda = FLAGS_lambda;
    average_degree = (double)num_edges * 2 / num_vertices;
    assigned_edges = 0;
    capacity = (double)num_edges * BALANCE_RATIO / p + 1;
    occupied.assign(p, 0);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    is_cores.assign(p, dense_bitset(num_vertices));
    is_boundarys.assign(p, dense_bitset(num_vertices));
//    master.assign(num_vertices, -1); // we do not use master list, as HEP does not either
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading...";
    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";
    adj_out.build(edges);
    adj_in.build_reverse(edges);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);

    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    read_timer.stop();
    LOG(INFO) << "time used for graph input and construction: " << read_timer.get_time();
};

void NePartitioner::assign_remaining()
{
    auto &is_boundary = is_boundarys[p - 1], &is_core = is_cores[p - 1];
    repv (u, num_vertices)
        for (auto &i : adj_out[u])
            if (edges[i.v].valid()) {
                assign_edge(p - 1, u, edges[i.v].second);
                is_boundary.set_bit_unsync(u);
                is_boundary.set_bit_unsync(edges[i.v].second);
            }

    repv (i, num_vertices) {
        if (is_boundary.get(i)) {
            is_core.set_bit_unsync(i);
            rep (j, p - 1)
                if (is_cores[j].get(i)) {
                    is_core.set_unsync(i, false);
                    break;
                }
        }
    }
}


// we do not assign master vertex, as HEP does not either
/*void NePartitioner::assign_master()
{
    std::vector<vid_t> count_master(p, 0);
    std::vector<vid_t> quota(p, num_vertices);
    long long sum = p * num_vertices;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::vector<dense_bitset::iterator> pos(p);
    rep (b, p)
        pos[b] = is_boundarys[b].begin();
    vid_t count = 0;
    while (count < num_vertices) {
        long long r = distribution(gen) * sum;
        int k;
        for (k = 0; k < p; k++) {
            if (r < quota[k])
                break;
            r -= quota[k];
        }
        while (pos[k] != is_boundarys[k].end() && master[*pos[k]] != -1)
            pos[k]++;
        if (pos[k] != is_boundarys[k].end()) {
            count++;
            master[*pos[k]] = k;
            writer.save_vertex(*pos[k], k);
            count_master[k]++;
            quota[k]--;
            sum--;
        }
    }
    int max_masters =
        *std::max_element(count_master.begin(), count_master.end());
    LOG(INFO) << "master balance: "
              << (double)max_masters / ((double)num_vertices / p);
}
*/

size_t NePartitioner::count_mirrors()
{
    size_t result = 0;
    rep (i, p)
        result += is_boundarys[i].popcount();
    return result;
}

void NePartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer compute_timer;

    min_heap.reserve(num_vertices);

    LOG(INFO) << "partitioning...";
    compute_timer.start();
    for (bucket = 0; bucket < p - 1; bucket++) {
        std::cerr << bucket << ", ";
        DLOG(INFO) << "sample size: " << adj_out.num_edges();
        while (occupied[bucket] < capacity) {
            vid_t d, vid;
            if (!min_heap.get_min(d, vid)) {
                if (!get_free_vertex(vid)) {
                    DLOG(INFO) << "partition " << bucket
                               << " stop: no free vertices";
                    break;
                }
                d = adj_out[vid].size() + adj_in[vid].size();
            } else {
                min_heap.remove(vid);
                /* CHECK_EQ(d, adj_out[vid].size() + adj_in[vid].size()); */
            }

            occupy_vertex(vid, d);
        }
        min_heap.clear();
        rep (direction, 2)
            repv (vid, num_vertices) {
                ne_adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size();) {
                    if (edges[neighbors[i].v].valid()) {
                        i++;
                    } else {
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    }
                }
            }
    }
    bucket = p - 1;
    std::cerr << bucket << std::endl;
    assign_remaining();
//    assign_master(); // HEP also does not compute master vertices, so for comparison we skip this step in NE

    size_t num_h2h_edges = 0;
    if (hybrid_partitioning){
    	// perform HDRF or random streaming on the high-degree edge file
    	if (FLAGS_random_streaming){
    		num_h2h_edges = random_streaming();
    	} else {
    		num_h2h_edges = hdrf_streaming();
    	}
    }


    LOG(INFO) << "Assigned " << num_h2h_edges << " edges by streaming." << std::endl;
    num_edges = num_edges + num_h2h_edges;
    compute_timer.stop();
    LOG(INFO) << "expected edges in each partition: " << num_edges / p;
    rep (i, p)
        DLOG(INFO) << "edges in partition " << i << ": " << occupied[i];
    size_t max_occupied = *std::max_element(occupied.begin(), occupied.end());
    LOG(INFO) << "balance: " << (double)max_occupied / ((double)num_edges / p);
    size_t total_mirrors = count_mirrors();
    LOG(INFO) << "total mirrors: " << total_mirrors;
    LOG(INFO) << "replication factor: " << (double)total_mirrors / num_vertices;
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();

    CHECK_EQ(assigned_edges, num_edges);

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
}

size_t NePartitioner::hdrf_streaming(){

	LOG(INFO) << "Streaming using HDRF algorithm." << std::endl;
	// assign the edges between two high degree vertices
    std::ifstream h2h_file(h2hedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate);
    auto filesize = h2h_file.tellg();
    LOG(INFO) << "file size: " << filesize;
    h2h_file.seekg(0, std::ios::beg);

	size_t total_h2h_edges = filesize / sizeof(size_t);

	// need to adapt capacity. old capacity did not take into account the streaming edges.
    capacity = (double)(num_edges + total_h2h_edges) * BALANCE_RATIO / p + 1;

	std::vector<edge_t> tmp_edges; // temporary buffer to read edges from file
	size_t chunk_size;
	size_t left_h2h_edges = total_h2h_edges;

	if (left_h2h_edges >= 100000){
		chunk_size = 100000; // batch read of so many edges
	}
	else {
		chunk_size = left_h2h_edges;
	}
	tmp_edges.resize(chunk_size);

	//	LOG(INFO) << "Chunk size is " << chunk_size << endl;

		/*
		 * init min_size
		 */
	min_size = UINT64_MAX;
	for (int i = 0; i < p; i++){
		if (occupied[i] < min_size){
			min_size = occupied[i];
		}
	}

//	LOG(INFO) << "Smallest partition has size " << min_size << std::endl;

	while (left_h2h_edges > 0){ // edges to be read
		h2h_file.read((char *)&tmp_edges[0], sizeof(edge_t) * chunk_size);
		for (size_t i = 0; i < chunk_size; i++){

			bucket = best_scored_partition(tmp_edges[i].first, tmp_edges[i].second); // according to HDRF scoring

	//			LOG(INFO) << "Assigning to p " << bucket << " : " << tmp_edges[i].first << ", "  << tmp_edges[i].second << endl;

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


	return total_h2h_edges;
}

double NePartitioner::compute_partition_score(vid_t u, vid_t v, int bucket_id) {
	if (occupied[bucket_id] >= capacity){
//		LOG(INFO) << "partition " << bucket_id << " is full with " << occupied[bucket_id] << std::endl;
		return -1.0; // partition is full, do not choose it
	}
	size_t degree_u = degrees[u];
	size_t degree_v = degrees[v];
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

//	cout << "gu " << gu << " gv " << gv << endl;
//	cout << "max_size " << max_size << " occupied " << occupied[bucket_id] << " min_size " << min_size << " bal " << bal <<endl;

	double score = gu + gv + lambda * bal;
	return score;
}


double NePartitioner::compute_partition_score_tiebreaking_balance(vid_t u, vid_t v, int bucket_id) {
	if (occupied[bucket_id] >= capacity){
//		cout << "partition " << bucket_id << " is full with " << occupied[bucket_id] << endl;
		return -1.0; // partition is full, do not choose it
	}
	size_t degree_u = degrees[u];
	size_t degree_v = degrees[v];
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

	double score = gu + gv;
	return score;
}

int NePartitioner::best_scored_partition_tiebreaking_balance(vid_t u, vid_t v) {
	double best_score = -1.0;
	int best_partition = 0;
	for (int i = 0; i < p; i++){
		double score = compute_partition_score_tiebreaking_balance(u, v, i);
//		cout << "score for partition " << i << " is " << score << endl;
		if (score == best_score){
			// use the tiebreaker rule: the smaller partition wins
			if (occupied[i] < occupied[best_partition]){
				best_score = score;
				best_partition = i;
			}
		}
		else if (score > best_score){
			best_score = score;
			best_partition = i;
		}
	}
	return best_partition;
}

int NePartitioner::best_scored_partition(vid_t u, vid_t v) {
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

size_t NePartitioner::random_streaming(){

	LOG(INFO) << "Streaming randomly." << std::endl;
	// assign the edges between two high degree vertices
    std::ifstream h2h_file(h2hedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate);
    auto filesize = h2h_file.tellg();
    LOG(INFO) << "file size: " << filesize;
    h2h_file.seekg(0, std::ios::beg);

	size_t total_h2h_edges = filesize / sizeof(size_t);


	std::vector<edge_t> tmp_edges; // temporary buffer to read edges from file
	size_t chunk_size;
	size_t left_h2h_edges = total_h2h_edges;

	if (left_h2h_edges >= 100000){
		chunk_size = 100000; // batch read of so many edges
	}
	else {
		chunk_size = left_h2h_edges;
	}
	tmp_edges.resize(chunk_size);

	//	LOG(INFO) << "Chunk size is " << chunk_size << endl;


	while (left_h2h_edges > 0){ // edges to be read
		h2h_file.read((char *)&tmp_edges[0], sizeof(edge_t) * chunk_size);
		for (size_t i = 0; i < chunk_size; i++){

			bucket = std::rand() % p; // random bucket

	//			LOG(INFO) << "Assigning to p " << bucket << " : " << tmp_edges[i].first << ", "  << tmp_edges[i].second << endl;

			assign_edge(bucket, tmp_edges[i].first, tmp_edges[i].second);
			is_boundarys[bucket].set_bit_unsync(tmp_edges[i].first);
			is_boundarys[bucket].set_bit_unsync(tmp_edges[i].second);

		}

		left_h2h_edges -= chunk_size;
		if (left_h2h_edges < chunk_size){ // adapt chunk size for last batch read
		 chunk_size = left_h2h_edges;
		}
	}

	return total_h2h_edges;
}



