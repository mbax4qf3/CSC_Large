#include <unistd.h>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <typeinfo>

#include "macros.h"
#include "u_io.h"
#include "u_label.h"
#include "u_spc.h"

bool ReadBool(const char tmp) {
    ASSERT('y' == tmp || 'n' == tmp);
    return 'y' == tmp;
}

int main(int argc, char** argv) {
    std::string lfilename;  // label file
    std::string ifilename;  // inserted edge file
    std::string dfilename;  // deleted edge file
    std::string qfilename;  // query file
    std::string afilename;  // answer file
    std::string minimality; // minimality enable
    std::string qtype;      // query type

    int option = -1;
    while (-1 != (option = getopt(argc, argv, "l:i:d:m:q:a:t:"))) {
        switch (option) {
            case 'l':
                lfilename = optarg; break;
            case 'i':
                ifilename = optarg; break;
            case 'd':
                dfilename = optarg; break;
            case 'm':
                minimality = optarg; break;
            case 'q':
                qfilename = optarg; break;
            case 'a':
                afilename = optarg; break;
            case 't':
                qtype = optarg; break; // o-ori; b-bi; i-inc; d-dec; l-lack(removed edge graph)
        }
    }

    // graph files are placed in graph/
    std::string file_index = lfilename.substr(6,1);
    
    // read index
    spc::USPCQuery uspc;

    if (qtype == "b") {
        uspc.IndexRead(lfilename);
    }
    else if (qtype == "o") {
        uspc.IndexRead_Ori(lfilename);
    }
    else if (qtype == "i") {
        uspc.IndexRead(lfilename);
    }
    else if (qtype == "d") {
        uspc.IndexRead(lfilename);
    } 
    else if (qtype == "l") {
        uspc.IndexRead(lfilename);
    }

    // update index if inserting edge
    if (qtype == "i") {
        // query first
        // read queries
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

            // answer queries and write results
            std::ofstream outfile;
            std::string answer_file;
            answer_file = afilename + file_index + "_bi" + qf;

            outfile.open(answer_file.c_str());
            for (const auto query : queries) {
                uint32_t result_d;
                uint32_t result_c;
                const auto beg = std::chrono::steady_clock::now();

                std::tie(result_d, result_c) = uspc.Count_Bigraph(query);

                const auto end = std::chrono::steady_clock::now();
                const auto dif = end - beg;
                outfile << query << "\t" << result_d << "\t" << result_c << "\t" << std::chrono::duration<double, std::micro>(dif).count() << std::endl;
            }
            outfile.close();
        }

        // read inserted edges
        std::vector<std::pair<uint32_t, uint32_t>> inc_edge;
        uint32_t num_inc = 0;

        FILE* file1 = fopen(ifilename.c_str(), "r");
        fscanf(file1, "%" SCNu32, &num_inc);
        
        for (uint32_t q = 0; q < num_inc; ++q) {
            uint32_t v1,v2;
            fscanf(file1, "%" SCNu32 " %" SCNu32, &v1, &v2);
            inc_edge.push_back(std::make_pair(v1,v2));
        }
        fclose(file1);

        // incremental update
        std::vector<std::tuple<uint32_t,uint32_t,uint32_t>> cnt_list;

        int mini_fg = 0;
        if (minimality == "m") mini_fg = 1;

        const auto beg_i = std::chrono::steady_clock::now();
        for (auto inc_e : inc_edge) {
            cnt_list.push_back(uspc.IncEdge(inc_e.first,inc_e.second,mini_fg));
        }
        const auto end_i = std::chrono::steady_clock::now();
        const auto dif_i = end_i - beg_i;
        
        // write info: avg time; label add, label renew, label removed(minimality)
        std::string info;
        if (minimality == "m") {
            info = "info/" + file_index + "m_incEdgeCnt.txt";
        } else if (qtype == "i"){
            info = "info/" + file_index + "_incEdgeCnt.txt";
        }

        std::ofstream outfile;
        outfile.open(info.c_str());
        outfile << "Avg update time: " << std::chrono::duration<double, std::micro>(dif_i).count() / num_inc << std::endl;
        for (auto lc : cnt_list) {
            outfile << std::get<0>(lc) << " " << std::get<1>(lc) << " " << std::get<2>(lc) << std::endl;
        }
        outfile.close();
    } // incremental update ends

    // update index if deleting edge
    if (qtype == "d") {
        // decremental update
        std::vector<std::tuple<uint32_t,uint32_t,uint32_t,uint32_t>> cnt_list;

        std::vector<std::pair<uint32_t, uint32_t>> dec_edge;
        uint32_t num_dec = 0;

        FILE* file1 = fopen(dfilename.c_str(), "r");
        fscanf(file1, "%" SCNu32, &num_dec);
        
        for (uint32_t q = 0; q < num_dec; ++q) {
            uint32_t v1,v2;
            fscanf(file1, "%" SCNu32 " %" SCNu32, &v1, &v2);
            dec_edge.push_back(std::make_pair(v1,v2));
        }
        fclose(file1);

        std::ofstream timefile;
        timefile.open("info/" + file_index +"_dec_time.txt");
        int no = 1;
        const auto beg_d_all = std::chrono::steady_clock::now();

        for (auto dec_e : dec_edge) {
            const auto beg_d = std::chrono::steady_clock::now();

            cnt_list.push_back(uspc.DecEdge(dec_e.first,dec_e.second));

            const auto end_d = std::chrono::steady_clock::now();
            const auto dif_d = end_d - beg_d;
            // write each update time
            timefile << std::chrono::duration<double, std::micro>(dif_d).count() << std::endl;
        }

        const auto end_d_all = std::chrono::steady_clock::now();
        const auto dif_d_all = end_d_all - beg_d_all;

        // write average update time
        timefile << "Avg: " << std::chrono::duration<double, std::micro>(dif_d_all).count() / num_dec << std::endl;
        timefile.close();

        std::string info;
        if (qtype == "d") {
            info = "info/" + file_index + "_decEdgeCnt.txt";
        }

        std::ofstream outfile;
        outfile.open(info.c_str());

        // for each deleted edge, write the size of affected hubs set A and set B, number of labels removed and number of labels added
        for (auto lc : cnt_list) {
            outfile << std::get<0>(lc) << " " << std::get<1>(lc) << ", " << std::get<2>(lc) << " " << std::get<3>(lc) << std::endl;
        }
        outfile.close();

        uspc.IndexWrite("label/" + file_index + "_dec.txt");
        
    } // decremental update ends
    
    // read queries
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

        // answer queries and write results
        std::ofstream outfile;
        std::string answer_file;
        if (qtype == "b") 
            answer_file = afilename + file_index + "_bi" + qf;
        else if (qtype == "o")
            answer_file = afilename + file_index +"_ori" + qf; 
        else if (qtype == "i")
            answer_file = afilename + file_index +"_bi_inc" + qf;
        else if (qtype == "d") {
            answer_file = afilename + file_index +"_bi_dec" + qf; 
        }
        else if (qtype == "l") {
            answer_file = afilename + file_index +"_bi_lack" + qf; 
        }
        if (minimality == "m") 
            answer_file = afilename + file_index +"m_bi_inc" + qf; 

        outfile.open(answer_file.c_str());
        for (const auto query : queries) {
            uint32_t result_d;
            uint32_t result_c;
            const auto beg = std::chrono::steady_clock::now();

            if (qtype == "i" || qtype == "b" || qtype == "d" || qtype == "l") 
                std::tie(result_d, result_c) = uspc.Count_Bigraph(query);
            else
                std::tie(result_d, result_c) = uspc.Count_Ori_Cycle(query);

            const auto end = std::chrono::steady_clock::now();
            const auto dif = end - beg;
            outfile << query << "\t" << result_d << "\t" << result_c << "\t" << std::chrono::duration<double, std::micro>(dif).count() << std::endl;
        }
        outfile.close();
    }
}