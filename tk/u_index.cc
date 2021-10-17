#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <set>
#include <string>
#include <fstream>
#include <iostream>

#include "macros.h"
#include "u_io.h"
#include "u_label.h"
#include "u_spc.h"

bool ReadBool(const char tmp) {
    ASSERT('y' == tmp || 'n' == tmp);
    return 'y' == tmp;
}

int main(int argc, char** argv) {
    // initialize the options
    std::string gfilename;
    std::string qfilename;
    std::string afilename;
    std::string lfilename;
    std::string bi_en;
    std::string osname;
    spc::USPCIndex::OrderScheme os = spc::USPCIndex::OrderScheme::kInvalid;
    {
        int option = -1;
        while (-1 != (option = getopt(argc, argv, "g:q:a:l:o:b:"))) {
            switch (option) {
                case 'g':
                    gfilename = optarg; break;
                case 'q':
                    qfilename = optarg; break;
                case 'a':
                    afilename = optarg; break;
                case 'l':
                    lfilename = optarg; break;
                case 'b':
                    bi_en = optarg; break;
                case 'o':
                    osname = optarg;
                    if ("degree" == osname) {
                        os = spc::USPCIndex::OrderScheme::kDegree;
                    } else if ("bidegree" == osname){
                        os = spc::USPCIndex::OrderScheme::kBDegree;
                    } else{
                        os = spc::USPCIndex::OrderScheme::kInvalid;
                    }
                    break;
            }
        }
    }

    // read the graph
    uint32_t n, m;
    spc::Graph graph;
    spc::Graph rev_graph;
    GraphRead(gfilename, graph, rev_graph, n, m);
    spc::USPCIndex spc;
    spc.set_os(os);

    // if there is a query file, run naive BFS instead of building index
    if (!qfilename.empty()) {
        std::string file_index = gfilename.substr(6,1);
        std::string query_files[5] = {"_high.txt","_midhigh.txt","_midlow.txt","_low.txt","_bottom.txt",};
        for (auto qf : query_files) {
            std::string qfilename_cstr = qfilename + file_index + qf;
            FILE* file = fopen(qfilename_cstr.c_str(), "r");

            uint32_t num_queries = 0;
            fscanf(file, "%" SCNu32, &num_queries);
            std::vector<uint32_t> queries;
            for (uint32_t q = 0; q < num_queries; ++q) {
                uint32_t vertex;
                fscanf(file, "%" SCNu32, &vertex);
                queries.push_back(vertex);
            }
            fclose(file);

            std::ofstream outfile;
            std::string answer_file;
            answer_file = afilename + file_index + "_bfs" + qf;

            outfile.open(answer_file.c_str());
            for (const auto query : queries) {
                uint32_t result_d;
                uint32_t result_c;
                const auto beg = std::chrono::steady_clock::now();
                // BFS
                uint32_t sd = UINT32_MAX;
                uint32_t spc = 0;
                std::vector<uint32_t> D(n, UINT32_MAX);
                std::vector<uint32_t> C(n, 0);
                C[query] = 1;
                std::queue<uint32_t> Q({query});

                while (!Q.empty()) {
                    uint32_t v = Q.front(); Q.pop();

                    if (v == query && D[v] < sd) {
                        result_d = D[v];
                        result_c = C[v];
                        break;
                    }

                    for (const uint32_t w : graph[v]) {
                        uint32_t DD = D[v];
                        if (v == query) DD = 0;
                        if (D[w] < DD + 1) continue;
                        else if (D[w] > DD + 1) {
                            D[w] = DD + 1;
                            C[w] = C[v];
                            Q.push(w);
                        }
                        else if (D[w] == DD + 1) {
                            C[w] += C[v];
                        }  
                    }
                }
                
                if (D[query] == UINT32_MAX) {
                    result_d = 0;
                    result_c = 0;
                }
                // BFS end

                const auto end = std::chrono::steady_clock::now();
                const auto dif = end - beg;
                outfile << query << "\t" << result_d << "\t" << result_c << "\t" << std::chrono::duration<double, std::micro>(dif).count() << std::endl;
            }
            outfile.close();
        }
        return 0;
    } // naive BFS ends

    const auto beg = std::chrono::steady_clock::now();
    // build index
    if (bi_en == "y")
        spc.BuildIndex(graph, rev_graph);
    else
        spc.BuildIndex_Ori(graph, rev_graph);

    const auto end = std::chrono::steady_clock::now();
    const auto dif = end - beg;

    spc.MergeIndexAll();
    uint64_t label_cnt;
    if (bi_en == "y")
        label_cnt = spc.IndexWrite(lfilename);
    else
        label_cnt = spc.IndexWrite_Ori(lfilename);
    
    // info file write
    std::ofstream outfile;
    std::string info = "info/";
    info.append(gfilename.substr(6));
    outfile.open(info.c_str());

    if (osname == "degree" || osname == "bidegree") {
        outfile << "Order: Degree\n";
    } else {
        outfile << "Order: Other\n";
    }
    outfile << "Index construction costs: " << std::chrono::duration<double, std::micro>(dif).count() << std::endl;
    outfile << "Index number: " << label_cnt << std::endl;

    outfile.close();
}
