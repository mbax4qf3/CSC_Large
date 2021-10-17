#include "u_io.h"

#include <algorithm>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <set>

#include "macros.h"

void GraphRead(const std::string& filename,
                             spc::Graph& graph,
                             spc::Graph& rev_graph,
                             uint32_t& n, uint32_t& m) {
    ASSERT(graph.empty());
    ASSERT(rev_graph.empty());
    FILE* gfile = fopen(filename.c_str(), "r");
    fscanf(gfile, "%" SCNu32 " %" SCNu32, &n, &m);
    spc::NormalV(n);
    /* construct Graph */
    graph.resize(n);
    rev_graph.resize(n);
    std::set<std::pair<uint32_t, uint32_t>> edges;
    for (uint32_t e = 0; e < m; ++e) {
        uint32_t v1, v2;
        fscanf(gfile, "%" SCNu32 " %" SCNu32, &v1, &v2);
        if (!(v1 < n && v2 < n))
            printf("%u %u\n",v1,v2);
        ASSERT(v1 < n && v2 < n);
        graph[v1].push_back(v2);
        rev_graph[v2].push_back(v1);
        /* check */
        ASSERT(v1 != v2); // no self loop
        if (0 != edges.count({v1, v2})) continue; // no duplicated edge
        edges.insert({v1, v2});
    }
    for (uint32_t v = 0; v < n; ++v) {
        std::sort(graph[v].begin(), graph[v].end());
        std::sort(rev_graph[v].begin(), rev_graph[v].end());
    }
    fclose(gfile);
}
