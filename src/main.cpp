#include <string>

#include "hep_partitioner.hpp"
#include "ne_partitioner.hpp"
#include "util.hpp"
#include "Vertex2EdgePart.hpp"

DECLARE_bool(help);
DECLARE_bool(helpshort);

DEFINE_int32(p, 10, "number of partitions");
DEFINE_string(filename, "", "the file name of the input graph");
DEFINE_string(filetype, "edgelist",
              "the type of input file (HEP supports only 'edgelist')");
DEFINE_string(method, "hep", "partition method: hep or ne. v2e is a special case that just transforms a given vertex partitioning to an edge partitioning.");
DEFINE_bool(write_results, false, "Should the result be written to the output file or not. Take care, writing is in ASCII format is is really really slow in the current implementation.");
DEFINE_bool(write_low_degree_edgelist, false, "Should the list of edges incident to a low-degree vertex be written out to a file?");
DEFINE_double(hdf, 100, "High-degree factor: hdf * average_degree = high-degree threshold (hdth). Called \tau in the paper. Vertices with than hdth neighbors are treated specially in fast NE");
DEFINE_double(lambda, 1.1, "Lambda value to weigh in balancing score in streaming partitioning via HDRF");
DEFINE_bool(extended_metrics, false, "Display extended metrics in the result");
DEFINE_bool(random_streaming, false, "Use random streaming instead of HDRF in the second phase of HEP.");
DEFINE_bool(hybrid_NE, false, "Perform hybrid partitioning in HEP-style, but use NE instead of NE++ for the first phase.");


int main(int argc, char *argv[])
{
    std::string usage = "-filename <path to the input graph> "
                        "[-filetype <edgelist|adjlist>] "
                        "[-p <number of partitions>] ";
    google::SetUsageMessage(usage);
    google::ParseCommandLineNonHelpFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = 1; // output log to stderr
    if (FLAGS_help) {
        FLAGS_help = false;
        FLAGS_helpshort = true;
    }
    
    google::HandleCommandLineHelpFlags();
    Timer timer;
    timer.start();
    Partitioner *partitioner = NULL;
    if (FLAGS_method == "hep")
        partitioner = new HepPartitioner(FLAGS_filename);
    else if (FLAGS_method == "ne")
        	partitioner = new NePartitioner(FLAGS_filename);
    else if (FLAGS_method == "v2e")
    	partitioner = new Vertex2EdgePart(FLAGS_filename);
    else
        LOG(ERROR) << "unkown method: " << FLAGS_method;
    LOG(INFO) << "partition method: " << FLAGS_method;
    partitioner->split();

    timer.stop();
    LOG(INFO) << "total time: " << timer.get_time();
}
