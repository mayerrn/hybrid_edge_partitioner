#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "util.hpp"

template <typename vertex_type, typename proc_type>
struct edgepart_writer {
    std::ofstream fout;

    edgepart_writer(const std::string &basefilename)
           : fout(partitioned_name(basefilename))
       {

       }

    void write_edge_assignment(vertex_type from, vertex_type to, proc_type proc)
    {
        fout << from << " " << to << " " << proc << std::endl;
    }
};
