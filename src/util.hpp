#pragma once

#include <utility>
#include <chrono>
#include <stdint.h>
#include <sys/stat.h>

#include <gflags/gflags.h>
#include <glog/logging.h>

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define repv(i, n) for (vid_t i = 0; i < n; ++i)
#define rept(i, n) for(int i=1; i<=(int)(n); ++i)

DECLARE_int32(p);
DECLARE_double(hdf);
DECLARE_string(filename);
DECLARE_string(filetype);
DECLARE_double(lambda);
DECLARE_bool(write_low_degree_edgelist);
DECLARE_bool(write_results);
DECLARE_bool(extended_metrics);
DECLARE_bool(random_streaming);
DECLARE_bool(hybrid_NE);

typedef uint32_t vid_t;
const vid_t INVALID_VID = -1;

struct edge_t {
    vid_t first, second;
    edge_t() : first(0), second(0) {}
    edge_t(vid_t first, vid_t second) : first(first), second(second) {}
    const bool valid() { return first != INVALID_VID; }
    void remove() { first = INVALID_VID; }
};

void preada(int f, char *buf, size_t nbytes, size_t off);
void reada(int f, char *buf, size_t nbytes);
void writea(int f, char *buf, size_t nbytes);

inline std::string h2hedgelist_name(const std::string &basefilename)
{
	std::stringstream ss;
	ss << basefilename << ".h2h_edgelist";
	return ss.str();
}

inline std::string lowedgelist_name(const std::string &basefilename)
{
	std::stringstream ss;
	ss << basefilename << ".low_edgelist";
	return ss.str();
}

inline std::string binedgelist_name(const std::string &basefilename)
{
    std::stringstream ss;
    ss << basefilename << ".binedgelist";
    return ss.str();
}
inline std::string degree_name(const std::string &basefilename)
{
    std::stringstream ss;
    ss << basefilename << ".degree";
    return ss.str();
}

inline std::string partitioned_name(const std::string &basefilename)
{
    std::stringstream ss;
    ss << basefilename << ".edgepart." << FLAGS_p;
    return ss.str();
}

inline bool is_exists(const std::string &name)
{
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

class Timer
{
  private:
    std::chrono::system_clock::time_point t1, t2;
    double total;

  public:
    Timer() : total(0) {}
    void reset() { total = 0; }
    void start() { t1 = std::chrono::system_clock::now(); }
    void stop()
    {
        t2 = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = t2 - t1;
        total += diff.count();
    }
    double get_time() { return total; }
};
