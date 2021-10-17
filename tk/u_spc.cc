#include "u_spc.h"
#include <stdlib.h>
#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <queue>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <tuple>
#include <stack>
#include <set>
#include <ctime>

#include "macros.h"

namespace spc {

std::tuple<uint32_t, uint32_t> USPCIndex::BFS_Count(const Graph& const_graph, uint32_t node) 
{
    uint32_t sd = UINT32_MAX;
    uint32_t spc = 0;
    std::vector<uint32_t> D(n_, UINT32_MAX);
    std::vector<uint32_t> C(n_, 0);
    std::queue<uint32_t> Q;

    for (auto vq_nbr:G_[node]) {
        Q.push(vq_nbr);
        D[vq_nbr] = 1;
        C[vq_nbr] = 1;
    }

    while (!Q.empty()) {
        uint32_t v = Q.front(); Q.pop();

        if (v == node) {
            return std:: make_tuple(D[v],C[v]);
        }
        for (const uint32_t w : const_graph[v]) {
            if (D[w] > D[v] + 1) {
                D[w] = D[v] + 1;
                C[w] = C[v];
                Q.push(w);
            }
            else if (D[w] == D[v] + 1) {
                C[w] += C[v];
            }
        }
    }
    return std::make_tuple(0,0);
}

// construct an index for the graph "const_graph", Original graph construction
void USPCIndex::BuildIndex_Ori(const Graph& const_graph, const Graph& rev_graph) 
{
    ASSERT(dL_.empty() && cL_.empty() && G_.empty() && rev_dL_.empty() && rev_cL_.empty() && rev_G_.empty());
    
    G_ = const_graph; rev_G_ = rev_graph;     // graph and reverse graph
    n_ = G_.size();                           // # of nodes
    
    // initialization
    deg_G_.resize(n_);                        //degree of each node
    dL_.resize(n_); cL_.resize(n_);           // cannonical and non-can in-label
    rev_dL_.resize(n_); rev_cL_.resize(n_);   // cannonical and non-can out-label
    order_.resize(n_); rank_.resize(n_);      // order: order[rank]=node; rank: rank[node]=rank
    // calculate the degree for each pair of node; 
    // let deg_G_[i] == deg_G_[i + n_/2], so that the pair shares the same(consecutive) rank 
    for (uint32_t i = 0; i < n_; ++i) {
        deg_G_[i] = G_[i].size() * rev_G_[i].size();
    }
    (this->*of_[os_])(deg_G_);
    OrderRank();

    // some auxiliary structures
    std::vector<uint32_t> dLu(n_, UINT32_MAX);    // current temp label
    std::vector<uint32_t> D(n_, UINT32_MAX);      // distance during BFS
    std::vector<uint32_t> C(n_, 0);               // counting during BFS

    // hub pushing
    for (size_t i = 0; i < n_; ++i) {
        const uint32_t u = order_[i];

        //****************************** FORWARD ******************************

        // for fast distance computation 
        for (const auto e : rev_dL_[u]) dLu[LEExtractV(e)] = LEExtractD(e);

        // bfs
        std::vector<uint32_t> reset({u});
        std::queue<uint32_t> Q({u});
        D[u] = 0; C[u] = 1;
        while (!Q.empty()) {
            const uint32_t v = Q.front(); Q.pop();
            const uint32_t dSoFar = Distance(dLu, dL_[v]);
            if (D[v] > dSoFar) continue;

            // add a corresponding entry 
            NormalD(D[v]); NormalC(C[v]);
            (D[v] < dSoFar? dL_[v] : cL_[v]).push_back(LEMerge(u, D[v], C[v]));

            // bfs for each neighbor, jump if required
            for (const uint32_t w : G_[v]) {
                if (rank_[w] <= rank_[u]) continue;
                if (UINT32_MAX == D[w]) {
                    D[w] = D[v] + 1;
                    C[w] = C[v];
                    Q.push(w);
                    reset.push_back(w);
                } else if (D[w] == D[v] + 1) {
                    if (likely(kUBC - C[v] >= C[w])) C[w] += C[v];
                    else C[w] = kUBC;
                }
            }
        }

        // clear 
        for (const uint32_t v : reset) {
            D[v] = UINT32_MAX; C[v] = 0;
        }

        for (const auto e : rev_dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;

        //****************************** REVERSE ******************************

        // for fast distance computation
        for (const auto e : dL_[u]) dLu[LEExtractV(e)] = LEExtractD(e);

        // bfs
        std::vector<uint32_t> rev_reset({u});
        std::queue<uint32_t> rev_Q({u});
        D[u] = 0; C[u] = 1;
        while (!rev_Q.empty()) {
            const uint32_t v = rev_Q.front(); rev_Q.pop();
            const uint32_t dSoFar = Distance(dLu, rev_dL_[v]);
            if (D[v] > dSoFar) continue;

            // add a corresponding entry
            NormalD(D[v]); NormalC(C[v]);
            (D[v] < dSoFar? rev_dL_[v] : rev_cL_[v]).push_back(LEMerge(u, D[v], C[v]));

            // bfs for each neighbor, jump if required
            for (const uint32_t w : rev_G_[v]) {
                if (rank_[w] <= rank_[u]) continue;
                if (UINT32_MAX == D[w]) {
                    D[w] = D[v] + 1;
                    C[w] = C[v];
                    rev_Q.push(w);
                    rev_reset.push_back(w);
                } else if (D[w] == D[v] + 1) {
                    if (likely(kUBC - C[v] >= C[w])) C[w] += C[v];
                    else C[w] = kUBC;
                }
            }
        }

        // clear
        for (const uint32_t v : rev_reset) {
            D[v] = UINT32_MAX; C[v] = 0;
        }

        for (const auto e : dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;
    }
}

// construct a bi-index for the graph "const_graph"; Bipartite Conversion construction
void USPCIndex::BuildIndex(const Graph& const_graph, const Graph& rev_graph) 
{
    ASSERT(dL_.empty() && cL_.empty() && G_.empty() && rev_dL_.empty() && rev_cL_.empty() && rev_G_.empty());
    
    G_ = const_graph; rev_G_ = rev_graph;     // graph and reverse graph
    n_ = G_.size();                           // # of nodes
    
    // initialization 
    // inv_in_.resize(n_);
    // inv_out_.resize(n_);
    deg_G_.resize(n_/2);                      //degree of each node
    dL_.resize(n_); cL_.resize(n_);           // cannonical and non-can in-label
    rev_dL_.resize(n_); rev_cL_.resize(n_);   // cannonical and non-can out-label
    order_.resize(n_/2); rank_.resize(n_);    // order: order[rank]=node; rank: rank[node]=rank
    // calculate the degree for each pair of node 
    // let deg_G_[i] == deg_G_[i + n_/2], so that the pair shares the same(consecutive) rank 
    for (uint32_t i = 0; i < n_/2; ++i) {
        deg_G_[i] = G_[i].size() * rev_G_[i + n_/2].size();
    }
    (this->*of_[os_])(deg_G_);
    OrderRank();

    // some auxiliary structures 
    std::vector<uint32_t> dLu(n_, UINT32_MAX);    // current temp label
    std::vector<uint32_t> D(n_, UINT32_MAX);      // distance during BFS
    std::vector<uint32_t> C(n_, 0);               // counting during BFS
    // hub pushing 
    for (size_t i = 0; i < n_; ++i) {
        const uint32_t u = order_[i];
        if (u < n_/2) {
            dL_[u].push_back(LEMerge(u, 0, 1));
            rev_dL_[u].push_back(LEMerge(u, 0, 1));
            continue;
        }

        //****************************** FORWARD ******************************

        // for fast distance computation 
        for (const auto e : rev_dL_[u-n_/2]) dLu[LEExtractV(e)] = LEExtractD(e);

        // bfs 
        std::vector<uint32_t> reset({u});
        std::queue<uint32_t> Q({u});
        D[u] = 0; C[u] = 1;
        while (!Q.empty()) {
            const uint32_t v = Q.front(); Q.pop();
            const uint32_t dSoFar = Distance(dLu, dL_[v-n_/2]); // big_hub to small_current
            if (D[v] > dSoFar) continue;

            // add a corresponding entry 
            NormalD(D[v]); NormalC(C[v]);
            // (D[v] < dSoFar? dL_[v] : cL_[v]).push_back(LEMerge(u, D[v], C[v]));
            // if (u != v) {inv_in_[u].push_back(v);}
            
            // jump from new v to old v directly 
            if (v >= n_/2) {
                (D[v] < dSoFar ? dL_[v-n_/2] : cL_[v-n_/2]).push_back(LEMerge(u, D[v]+1, C[v]));
                // if (u != v-n_/2) {inv_in_[u].push_back(v-n_/2);}
                D[v-n_/2] = D[v]+1;
                C[v-n_/2] = C[v];
                reset.push_back(v-n_/2);
            } 

            // bfs for each neighbor, jump if required
            uint32_t cur_v = v;
            if (v >= n_/2) cur_v = v - n_/2;
            for (const uint32_t w : G_[cur_v]) {
                if (rank_[w] <= rank_[u]) continue;
                if (UINT32_MAX == D[w]) {
                    D[w] = D[cur_v] + 1;
                    C[w] = C[cur_v];
                    Q.push(w);
                    reset.push_back(w);
                } else if (D[w] == D[cur_v] + 1) {
                    if (likely(kUBC - C[cur_v] >= C[w])) C[w] += C[cur_v];
                    else C[w] = kUBC;
                }
            }
        }

        // clear
        for (const uint32_t v : reset) {
            D[v] = UINT32_MAX; C[v] = 0;
        }

        for (const auto e : rev_dL_[u-n_/2]) dLu[LEExtractV(e)] = UINT32_MAX;

        //****************************** REVERSE ******************************

        // for fast distance computation 
        for (const auto e : dL_[u-n_/2]) dLu[LEExtractV(e)] = LEExtractD(e);

        // bfs 
        std::vector<uint32_t> rev_reset({u});
        std::queue<uint32_t> rev_Q({u});
        D[u] = 0; C[u] = 1;
        while (!rev_Q.empty()) {
            const uint32_t v = rev_Q.front(); rev_Q.pop();
            if (u != v) {
	            const uint32_t dSoFar = Distance(dLu, rev_dL_[v]) - 1;
	            if (D[v] > dSoFar) continue;

	            // add a corresponding entry 
	            NormalD(D[v]); NormalC(C[v]);
	            (D[v] < dSoFar? rev_dL_[v] : rev_cL_[v]).push_back(LEMerge(u, D[v], C[v]));
	            // if (u != v) {inv_out_[u].push_back(v);}
	            // jump from old v to new v directly 
	            if (v < n_/2 && v+n_/2 != u) {
	                // (D[v] < dSoFar? rev_dL_[v+n_/2] : rev_cL_[v+n_/2]).push_back(LEMerge(u, D[v]+1, C[v]));
	                // if (u != v+n_/2) {inv_out_[u].push_back(v+n_/2);}
	                D[v+n_/2] = D[v]+1;
	                C[v+n_/2] = C[v];
	                rev_reset.push_back(v+n_/2);
	            }
	        }

            // bfs for each neighbor, jump if required
            uint32_t cur_v = v;
            if (v < n_/2 && v+n_/2 != u) cur_v = v + n_/2;
            for (const uint32_t w : rev_G_[cur_v]) {
                if (rank_[w] <= rank_[u]) continue;
                if (UINT32_MAX == D[w]) {
                    D[w] = D[cur_v] + 1;
                    C[w] = C[cur_v];
                    rev_Q.push(w);
                    rev_reset.push_back(w);
                } else if (D[w] == D[cur_v] + 1) {
                    if (likely(kUBC - C[cur_v] >= C[w])) C[w] += C[cur_v];
                    else C[w] = kUBC;
                }
            }
        }

        // clear 
        for (const uint32_t v : rev_reset) {
            D[v] = UINT32_MAX; C[v] = 0;
        }

        for (const auto e : dL_[u-n_/2]) dLu[LEExtractV(e)] = UINT32_MAX;
    }
}

// Find the position of the particular label entry
uint32_t USPCQuery::Find_Can_Noncan_DC(uint32_t hub, uint32_t v, int in_out) 
{
    if (in_out == 0) {
        int32_t left = 0; int32_t right = cL_[v].size() - 1; int32_t mid;
        while (left <= right) {
            mid = (left + right) / 2;
            if (rank_[hub] == rank_[LEExtractV(cL_[v][mid])]) return mid;
            else if (rank_[hub] < rank_[LEExtractV(cL_[v][mid])]) right = mid - 1;
            else left = mid + 1;
        }
    } else {
        int32_t left = 0; int32_t right = rev_cL_[v].size() - 1; int32_t mid;
        while (left <= right) {
            mid = (left + right) / 2;
            if (rank_[hub] == rank_[LEExtractV(rev_cL_[v][mid])]) return mid;
            else if (rank_[hub] < rank_[LEExtractV(rev_cL_[v][mid])]) right = mid - 1;
            else left = mid + 1;
        }
    }

    return -1;
}

// Count SPC(v1,v2)
std::pair<uint32_t, uint64_t> USPCQuery::Count(uint32_t v1, uint32_t v2) const 
{
    // ASSERT(v1 != v2);
    if (v1 == v2) {
        return std::make_pair(0,1);
    }
    // count the # of shortest paths
    uint32_t sp_d = UINT32_MAX;
    uint64_t sp_c = 0;
    size_t p1 = 0, p2 = 0;

    while (p1 < rev_cL_[v1].size() && p2 < cL_[v2].size()) {
        const uint32_t w1 = LEExtractV(rev_cL_[v1][p1]);
        const uint32_t w2 = LEExtractV(cL_[v2][p2]);
        
        if (rank_[w1] < rank_[w2]) ++p1;
        else if (rank_[w1] > rank_[w2]) ++p2;
        else {
            const uint32_t d = LEExtractD(rev_cL_[v1][p1]) +
                    LEExtractD(cL_[v2][p2]);
            if (d < sp_d) {
                sp_d = d;
                sp_c = static_cast<uint64_t>(LEExtractC(rev_cL_[v1][p1])) *
                        LEExtractC(cL_[v2][p2]);
            } else if (d == sp_d) {
                uint64_t c = static_cast<uint64_t>(LEExtractC(rev_cL_[v1][p1])) *
                        LEExtractC(cL_[v2][p2]);
                sp_c += c;
            }
            ++p1; ++p2;        
        }
    }

    std::pair<uint32_t,uint64_t> dc = std::make_pair(sp_d,sp_c);
    return dc;
}

// Count SCCnt(node) (BI)
std::tuple<uint32_t, uint32_t> USPCQuery::Count_Bigraph(uint32_t node) const 
{
    uint32_t v1 = node;
    uint32_t v2 = node+n_/2;

    std::vector<uint32_t> hub;
    std::vector<uint32_t> dis;
    std::vector<uint32_t> cnt;
    if (G_[v1].size() == 0 || rev_G_[v2].size() == 0) {
        return std::make_tuple(0, 0);
    }
    // count the # of shortest paths
    uint32_t sp_d = UINT32_MAX;
    uint64_t sp_c = 0;
    size_t p1 = 0, p2 = 0;
    while (p1 < rev_cL_[v1].size() && p2 < cL_[v2].size()) {
        const uint32_t w1 = LEExtractV(rev_cL_[v1][p1]);
        const uint32_t w2 = LEExtractV(cL_[v2][p2]);

        if (rank_[w1] < rank_[w2]) ++p1;
        else if (rank_[w1] > rank_[w2]) ++p2;
        else {
            const uint32_t d = LEExtractD(rev_cL_[v1][p1]) +
                    LEExtractD(cL_[v2][p2]);
            if (d < sp_d) {
                sp_d = d;
                sp_c = static_cast<uint64_t>(LEExtractC(rev_cL_[v1][p1])) *
                        LEExtractC(cL_[v2][p2]);
                hub.clear(); dis.clear(); cnt.clear();
                hub.push_back(w1);
                dis.push_back(LEExtractD(rev_cL_[v1][p1])); dis.push_back(LEExtractD(cL_[v2][p2]));
                cnt.push_back(LEExtractC(rev_cL_[v1][p1])); cnt.push_back(LEExtractC(cL_[v2][p2]));
            } else if (d == sp_d) {
                uint64_t c = static_cast<uint64_t>(LEExtractC(rev_cL_[v1][p1])) *
                        LEExtractC(cL_[v2][p2]);
                sp_c += c;
                hub.push_back(w1);
                dis.push_back(LEExtractD(rev_cL_[v1][p1])); dis.push_back(LEExtractD(cL_[v2][p2]));
                cnt.push_back(LEExtractC(rev_cL_[v1][p1])); cnt.push_back(LEExtractC(cL_[v2][p2]));
            }
            ++p1; ++p2;
        }
    }

    if (sp_d == UINT32_MAX) sp_d = 0;
    return std::make_tuple((sp_d+1)/2, sp_c);
}

// Count SCCnt(node) (ORI)
std::tuple<uint32_t, uint32_t> USPCQuery::Count_Ori_Cycle(uint32_t node) const 
{
    if (G_[node].size() == 0 || rev_G_[node].size() == 0) {
        return std::make_tuple(0, 0);
    }

    // count the # of shortest paths
    uint32_t sp_d = UINT32_MAX;
    uint64_t sp_c = 0;
    if (G_[node].size() <= rev_G_[node].size()) {
        
        for (size_t i = 0; i < G_[node].size(); ++i) {
            size_t p1 = 0, p2 = 0;
            uint32_t vo = G_[node][i];

            while (p1 < rev_cL_[vo].size() && p2 < cL_[node].size()) {
                const uint32_t w1 = LEExtractV(rev_cL_[vo][p1]);
                const uint32_t w2 = LEExtractV(cL_[node][p2]);

                if (rank_[w1] < rank_[w2]) ++p1;
                else if (rank_[w1] > rank_[w2]) ++p2;
                else {
                    const uint32_t d = LEExtractD(rev_cL_[vo][p1]) +
                            LEExtractD(cL_[node][p2]);
                    if (d < sp_d) {
                        sp_d = d;
                        sp_c = static_cast<uint64_t>(LEExtractC(rev_cL_[vo][p1])) *
                                LEExtractC(cL_[node][p2]);
                    } else if (d == sp_d) {
                        uint64_t c = static_cast<uint64_t>(LEExtractC(rev_cL_[vo][p1])) *
                                LEExtractC(cL_[node][p2]);
                        sp_c += c;
                    }
                    ++p1; ++p2;
                }
            } //while

        } // for each neighbor
        if (sp_d < UINT32_MAX) {
            sp_d += 1;
        }
    } else {
            for (size_t i = 0; i < rev_G_[node].size(); ++i) {
            size_t p1 = 0, p2 = 0;
            uint32_t vi = rev_G_[node][i];

            while (p1 < rev_cL_[node].size() && p2 < cL_[vi].size()) {
                const uint32_t w1 = LEExtractV(rev_cL_[node][p1]);
                const uint32_t w2 = LEExtractV(cL_[vi][p2]);

                if (rank_[w1] < rank_[w2]) ++p1;
                else if (rank_[w1] > rank_[w2]) ++p2;
                else {
                    const uint32_t d = LEExtractD(rev_cL_[node][p1]) +
                            LEExtractD(cL_[vi][p2]);
                    if (d < sp_d) {
                        sp_d = d;
                        sp_c = static_cast<uint64_t>(LEExtractC(rev_cL_[node][p1])) *
                                LEExtractC(cL_[vi][p2]);
                    } else if (d == sp_d) {
                        uint64_t c = static_cast<uint64_t>(LEExtractC(rev_cL_[node][p1])) *
                                LEExtractC(cL_[vi][p2]);
                        sp_c += c;
                    }
                    ++p1; ++p2;
                }
            } //while

        } // for each neighbor
        if (sp_d < UINT32_MAX) {
            sp_d += 1;
        }
    }
    
    if (sp_d == UINT32_MAX) sp_d = 0;
    return std::make_tuple(sp_d, sp_c);
}

// Count SPC(v1,v2) and return true if use self as hub
std::tuple<uint32_t, uint64_t, bool> USPCQuery::Count_self(uint32_t v1, uint32_t v2) const 
{
    bool self_ = false;

    if (v1 == v2) {
        return std::make_tuple(0,1,true);
    }
    // count the # of shortest paths
    uint32_t sp_d = UINT32_MAX;
    uint64_t sp_c = 0;
    size_t p1 = 0, p2 = 0;

    while (p1 < rev_cL_[v1].size() && p2 < cL_[v2].size()) {
        const uint32_t w1 = LEExtractV(rev_cL_[v1][p1]);
        const uint32_t w2 = LEExtractV(cL_[v2][p2]);
        
        if (rank_[w1] < rank_[w2]) ++p1;
        else if (rank_[w1] > rank_[w2]) ++p2;
        else {
            const uint32_t d = LEExtractD(rev_cL_[v1][p1]) +
                    LEExtractD(cL_[v2][p2]);
            if (d < sp_d) {
                if (self_ == false && (w1 == v1 || w1 == v2)) self_ = true;
                if (self_ == true) self_ = false;
                sp_d = d;
                sp_c = static_cast<uint64_t>(LEExtractC(rev_cL_[v1][p1])) *
                        LEExtractC(cL_[v2][p2]);
            } else if (d == sp_d) {
                uint64_t c = static_cast<uint64_t>(LEExtractC(rev_cL_[v1][p1])) *
                        LEExtractC(cL_[v2][p2]);
                sp_c += c;
            }
            ++p1; ++p2;        
        }
    }

    std::tuple<uint32_t,uint64_t,bool> dcs = std::make_tuple(sp_d,sp_c,self_);
    return dcs;
}

// Count SPC(v1,v2) and return hubs
std::tuple<uint32_t, uint64_t, std::set<uint32_t>,bool> USPCQuery::Count_hub(uint32_t v1, uint32_t v2) const 
{
    uint32_t sp_d = UINT32_MAX;
    uint64_t sp_c = 0;

    std::set<uint32_t> hubs;
    bool self_ = false;

    if (v1 == v2) {
        sp_d = 0; sp_c = 1;
        return std::make_tuple(sp_d,sp_c,hubs,true);
    }

    size_t p1 = 0, p2 = 0;

    while (p1 < rev_cL_[v1].size() && p2 < cL_[v2].size()) {
        const uint32_t w1 = LEExtractV(rev_cL_[v1][p1]);
        const uint32_t w2 = LEExtractV(cL_[v2][p2]);
        
        if (rank_[w1] < rank_[w2]) ++p1;
        else if (rank_[w1] > rank_[w2]) ++p2;
        else {
            const uint32_t d = LEExtractD(rev_cL_[v1][p1]) +
                    LEExtractD(cL_[v2][p2]);
            if (d < sp_d) {
                hubs.clear();
                hubs.insert(w1);
                if ((!self_) && (w1 == v1 || w1 == v2)) {self_ = true;}
                else if (self_) {self_ = false;}
                sp_d = d;
                sp_c = static_cast<uint64_t>(LEExtractC(rev_cL_[v1][p1])) *
                        LEExtractC(cL_[v2][p2]);
            } else if (d == sp_d) {
                hubs.insert(w1);
                if ((!self_) && (w1 == v1 || w1 == v2)) {self_ = true;}

                uint64_t c = static_cast<uint64_t>(LEExtractC(rev_cL_[v1][p1])) *
                        LEExtractC(cL_[v2][p2]);
                sp_c += c;
            }
            ++p1; ++p2;        
        }
    }

    std::tuple<uint32_t,uint64_t,std::set<uint32_t>,bool> dcs = std::make_tuple(sp_d,sp_c,hubs,self_);
    return dcs;
}

// Count SPC(v1,v2)
std::uint32_t USPCQuery::Count_dis(uint32_t v1, uint32_t v2) const 
{
    if (v1 == v2) {
        return 0;
    }
    // count the # of shortest paths
    uint32_t sp_d = UINT32_MAX;
    size_t p1 = 0, p2 = 0;

    while (p1 < rev_cL_[v1].size() && p2 < cL_[v2].size()) {
        const uint32_t w1 = LEExtractV(rev_cL_[v1][p1]);
        const uint32_t w2 = LEExtractV(cL_[v2][p2]);
        
        if (rank_[w1] < rank_[w2]) ++p1;
        else if (rank_[w1] > rank_[w2]) ++p2;
        else {
            const uint32_t d = LEExtractD(rev_cL_[v1][p1]) + LEExtractD(cL_[v2][p2]);
            if (d < sp_d) {
                sp_d = d;
            }
            ++p1; ++p2;        
        }
    }

    return sp_d;
}

// Merge index before store in disk
void USPCIndex::MergeIndexAll()
{
    // merge non-canonical and canonical labels
    for (uint32_t i = 0; i < n_; ++i) {
        // out label
        std::vector<LabelEntry> rev_mL;
        rev_mL.reserve(rev_dL_[i].size() + rev_cL_[i].size());
        size_t rdi = 0, rci = 0;
        uint32_t label_cnt = 0;
        while (rdi < rev_dL_[i].size() && rci < rev_cL_[i].size()) {
            if (rank_[LEExtractV(rev_dL_[i][rdi])] < rank_[LEExtractV(rev_cL_[i][rci])]) {
                rev_mL.push_back(rev_dL_[i][rdi++]);
                label_cnt++;
            } else {
                rev_mL.push_back(rev_cL_[i][rci++]);
                label_cnt++;
            }
        }
        while (rdi < rev_dL_[i].size()) {
            rev_mL.push_back(rev_dL_[i][rdi++]);
            label_cnt++;
        }
        while (rci < rev_cL_[i].size()) {
            rev_mL.push_back(rev_cL_[i][rci++]);
            label_cnt++;
        }
        rev_mL.resize(label_cnt);
        rev_cL_[i] = rev_mL;

        // in label
        uint32_t in = i;
        std::vector<LabelEntry> mL;
        mL.reserve(dL_[in].size() + cL_[in].size());
        size_t di = 0, ci = 0;
        label_cnt = 0;
        while (di < dL_[in].size() && ci < cL_[in].size()) {
            if (rank_[LEExtractV(dL_[in][di])] < rank_[LEExtractV(cL_[in][ci])]) {
                mL.push_back(dL_[in][di++]);
                label_cnt++;
            } else {
                mL.push_back(cL_[in][ci++]);
                label_cnt++;
            }
        }
        while (di < dL_[in].size()) {
            mL.push_back(dL_[in][di++]);
            label_cnt++;
        }
        while (ci < cL_[in].size()) {
            mL.push_back(cL_[in][ci++]);
            label_cnt++;
        }
        mL.resize(label_cnt);
        cL_[in] = mL; // in -> i
    }
    decltype(dL_)().swap(dL_);
    decltype(rev_dL_)().swap(rev_dL_);
}

// Read index from disk (ORI)
void USPCQuery::IndexRead_Ori(const std::string& filename) 
{
    ASSERT(dL_.empty() && cL_.empty() && G_.empty() && rev_dL_.empty() && rev_cL_.empty() && rev_G_.empty());
    FILE* file = fopen(filename.c_str(), "rb");
    ASSERT(fread(&n_, sizeof(n_), 1, file) == 1);
    G_.resize(n_);
    rev_G_.resize(n_);

    // read graph
    for (uint32_t u = 0; u < n_; ++u) {
        uint32_t s = 0;
        ASSERT(fread(&s, sizeof(s), 1, file) == 1);
        G_[u].resize(s);
        ASSERT(fread(G_[u].data(), sizeof(G_[u].back()), s, file) == s);
    }
    for (uint32_t u = 0; u < n_; ++u) {
        uint32_t s = 0;
        ASSERT(fread(&s, sizeof(s), 1, file) == 1);
        rev_G_[u].resize(s);
        ASSERT(fread(rev_G_[u].data(), sizeof(rev_G_[u].back()), s, file) == s);
    }

    rank_.resize(n_);
    ASSERT(fread(rank_.data(), sizeof(rank_.back()), n_, file) == n_);

    // initialization
    cL_.resize(n_);
    rev_cL_.resize(n_);
    // read labels
    uint64_t num_labels = 0;
    uint32_t c_s, rc_s;
    // read canonical and non-canonical labels
    // n_ -> n_ / 2
    for (uint32_t i = 0; i < n_; ++i) {

        ASSERT(fread(&rc_s, sizeof(rc_s), 1, file) == 1);
        rev_cL_[i].resize(rc_s);
        ASSERT(fread(rev_cL_[i].data(), sizeof(rev_cL_[i].back()), rc_s, file) == rc_s);

        num_labels += rc_s;

        ASSERT(fread(&c_s, sizeof(c_s), 1, file) == 1);
        cL_[i].resize(c_s);
        ASSERT(fread(cL_[i].data(), sizeof(cL_[i].back()), c_s, file) == c_s);

        num_labels += c_s;
    }

    fclose(file);
}

// Read index from disk (BI)
void USPCQuery::IndexRead(const std::string& filename) 
{
    ASSERT(dL_.empty() && cL_.empty() && G_.empty() && rev_dL_.empty() && rev_cL_.empty() && rev_G_.empty());
    FILE* file = fopen(filename.c_str(), "rb");
    ASSERT(fread(&n_, sizeof(n_), 1, file) == 1);
    G_.resize(n_);
    rev_G_.resize(n_);

    // read graph
    for (uint32_t u = 0; u < n_; ++u) {
        uint32_t s = 0;
        ASSERT(fread(&s, sizeof(s), 1, file) == 1);
        G_[u].resize(s);
        ASSERT(fread(G_[u].data(), sizeof(G_[u].back()), s, file) == s);
    }
    for (uint32_t u = 0; u < n_; ++u) {
        uint32_t s = 0;
        ASSERT(fread(&s, sizeof(s), 1, file) == 1);
        rev_G_[u].resize(s);
        ASSERT(fread(rev_G_[u].data(), sizeof(rev_G_[u].back()), s, file) == s);
    }

    rank_.resize(n_);
    ASSERT(fread(rank_.data(), sizeof(rank_.back()), n_, file) == n_);

    /* 
    inv_in_.resize(n_);
    inv_out_.resize(n_);
    for (uint32_t u = 0; u < n_; ++u) {
        uint32_t s = 0;
        ASSERT(fread(&s, sizeof(s), 1, file) == 1);
        inv_in_[u].resize(s);
        ASSERT(fread(inv_in_[u].data(), sizeof(inv_in_[u].back()), s, file) == s);
    }
    for (uint32_t u = 0; u < n_; ++u) {
        uint32_t s = 0;
        ASSERT(fread(&s, sizeof(s), 1, file) == 1);
        inv_out_[u].resize(s);
        ASSERT(fread(inv_out_[u].data(), sizeof(inv_out_[u].back()), s, file) == s);
    }
	*/

    // initialization
    cL_.resize(n_);
    rev_cL_.resize(n_);
    // read labels
    uint64_t num_labels = 0;
    uint32_t c_s, rc_s;
    // read canonical and non-canonical labels
    // n_ -> n_ / 2
    for (uint32_t i = 0; i < n_ / 2; ++i) {

        ASSERT(fread(&rc_s, sizeof(rc_s), 1, file) == 1);
        rev_cL_[i].resize(rc_s);
        ASSERT(fread(rev_cL_[i].data(), sizeof(rev_cL_[i].back()), rc_s, file) == rc_s);
        
        uint64_t dis1 = LEMerge(0,1,0).v_d_c;

        rev_cL_[i+n_/2].resize(rc_s);
        rev_cL_[i+n_/2] = rev_cL_[i];
        rev_cL_[i+n_/2].pop_back();
        if (rev_cL_[i+n_/2].size() >= 1) {
            if (LEExtractV(rev_cL_[i+n_/2].back()) == (i+n_/2)) {
                rev_cL_[i+n_/2].pop_back(); 
            }
        }
        for (uint32_t j = 0; j < rev_cL_[i+n_/2].size();++j) {
            uint32_t v    = LEExtractV(rev_cL_[i+n_/2][j]);
            uint32_t d    = LEExtractD(rev_cL_[i+n_/2][j]);
            uint32_t c    = LEExtractC(rev_cL_[i+n_/2][j]);
            rev_cL_[i+n_/2][j] = LEMerge(v,d+1,c);
        }
        rev_cL_[i+n_/2].push_back(LEMerge(i+n_/2,0,1));

        num_labels += rc_s;

        ASSERT(fread(&c_s, sizeof(c_s), 1, file) == 1);
        cL_[i].resize(c_s);
        ASSERT(fread(cL_[i].data(), sizeof(cL_[i].back()), c_s, file) == c_s);

        cL_[i+n_/2].resize(c_s);
        cL_[i+n_/2] = cL_[i];
        cL_[i+n_/2].pop_back();
        
        for (uint32_t j = 0; j < cL_[i+n_/2].size();++j) {
            uint32_t v    = LEExtractV(cL_[i+n_/2][j]);
            uint32_t d    = LEExtractD(cL_[i+n_/2][j]);
            uint32_t c    = LEExtractC(cL_[i+n_/2][j]);
            cL_[i+n_/2][j] = LEMerge(v,d-1,c);
        }

        num_labels += c_s;
    }

    fclose(file);
}

// Write index to disk (ORI)
uint64_t USPCIndex::IndexWrite_Ori(const std::string& filename) const 
{
    ASSERT(0 != n_);
    FILE* file = fopen(filename.c_str(), "wb");

    // write original graph
    fwrite(&n_, sizeof(n_), 1, file);
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = G_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(G_[u].data(), sizeof(G_[u].back()), s, file);
    }

    // write reversed graph
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = rev_G_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(rev_G_[u].data(), sizeof(rev_G_[u].back()), s, file);
    }

    // write rank info
    fwrite(rank_.data(), sizeof(rank_.back()), n_, file);

    uint64_t num_labels = 0;
    // read canonical and non-canonical labels
    for (uint32_t i = 0; i < n_; ++i) {

        const uint32_t rev_c_s = rev_cL_[i].size();
        fwrite(&rev_c_s, sizeof(rev_c_s), 1, file);
        fwrite(rev_cL_[i].data(), sizeof(rev_cL_[i].back()), rev_c_s, file);
        
        num_labels += rev_c_s;

        const uint32_t c_s = cL_[i].size();
        fwrite(&c_s, sizeof(c_s), 1, file);
        fwrite(cL_[i].data(), sizeof(cL_[i].back()), c_s, file);

        num_labels += c_s;

    }

    fclose(file);
    return num_labels;
}

// Write index to disk (BI)
uint64_t USPCIndex::IndexWrite(const std::string& filename) const 
{
    ASSERT(0 != n_);
    FILE* file = fopen(filename.c_str(), "wb");
    // write original graph
    fwrite(&n_, sizeof(n_), 1, file);
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = G_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(G_[u].data(), sizeof(G_[u].back()), s, file);
    }
    // write reversed graph
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = rev_G_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(rev_G_[u].data(), sizeof(rev_G_[u].back()), s, file);
    }

    // write rank info
    fwrite(rank_.data(), sizeof(rank_.back()), n_, file);

    /*
    // write inverted index hubs
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = inv_in_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(inv_in_[u].data(), sizeof(inv_in_[u].back()), s, file);
    }
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = inv_out_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(inv_out_[u].data(), sizeof(inv_out_[u].back()), s, file);
    }
    */

    uint64_t num_labels = 0;
    // read canonical and non-canonical labels
    // n_ -> n_ / 2
    for (uint32_t i = 0; i < n_/2; ++i) {

        const uint32_t rev_c_s = rev_cL_[i].size();
        fwrite(&rev_c_s, sizeof(rev_c_s), 1, file);
        fwrite(rev_cL_[i].data(), sizeof(rev_cL_[i].back()), rev_c_s, file);
        
        num_labels += rev_c_s;

        const uint32_t c_s = cL_[i].size();
        fwrite(&c_s, sizeof(c_s), 1, file);
        fwrite(cL_[i].data(), sizeof(cL_[i].back()), c_s, file);

        num_labels += c_s;

    }

    fclose(file);
    return num_labels;
}

uint64_t USPCQuery::IndexWrite(const std::string& filename) const 
{
    ASSERT(0 != n_);
    FILE* file = fopen(filename.c_str(), "wb");
    // write original graph
    fwrite(&n_, sizeof(n_), 1, file);
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = G_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(G_[u].data(), sizeof(G_[u].back()), s, file);
    }
    // write reversed graph
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = rev_G_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(rev_G_[u].data(), sizeof(rev_G_[u].back()), s, file);
    }

    // write rank info
    fwrite(rank_.data(), sizeof(rank_.back()), n_, file);

    /*
    // write inverted index hubs
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = inv_in_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(inv_in_[u].data(), sizeof(inv_in_[u].back()), s, file);
    }
    for (uint32_t u = 0; u < n_; ++u) {
        const uint32_t s = inv_out_[u].size();
        fwrite(&s, sizeof(s), 1, file);
        fwrite(inv_out_[u].data(), sizeof(inv_out_[u].back()), s, file);
    }
    */ 
    uint64_t num_labels = 0;
    // read canonical and non-canonical labels
    // n_ -> n_ / 2
    for (uint32_t i = 0; i < n_/2; ++i) {

        const uint32_t rev_c_s = rev_cL_[i].size();
        fwrite(&rev_c_s, sizeof(rev_c_s), 1, file);
        fwrite(rev_cL_[i].data(), sizeof(rev_cL_[i].back()), rev_c_s, file);
        
        num_labels += rev_c_s;

        const uint32_t c_s = cL_[i].size();
        fwrite(&c_s, sizeof(c_s), 1, file);
        fwrite(cL_[i].data(), sizeof(cL_[i].back()), c_s, file);

        num_labels += c_s;

    }

    fclose(file);
    return num_labels;
}

// Get distance (used in index construction)
uint32_t USPCIndex::Distance(const std::vector<uint32_t>& dLu,
                                                         const std::vector<LabelEntry>& dLv) const 
{
    uint32_t d = UINT32_MAX;
    for (const auto e : dLv) {
        const uint32_t v = LEExtractV(e);
        if (UINT32_MAX == dLu[v]) continue;
        const uint32_t dd = dLu[v] + LEExtractD(e);
        if (dd < d) d = dd;
    }
    return d;
}

// Get distance (used in fast distance query)
uint32_t USPCQuery::Distance(const std::vector<uint32_t>& dLu,
                                                         const std::vector<LabelEntry>& dLv) const 
{
    uint32_t d = UINT32_MAX;
    for (const auto e : dLv) {
        const uint32_t v = LEExtractV(e);
        if (UINT32_MAX == dLu[v]) continue;
        const uint32_t dd = dLu[v] + LEExtractD(e);
        if (dd < d) d = dd;
    }
    return d;
}

// Degree order
void USPCIndex::DegreeOrder(const std::vector<uint32_t>& graph_deg) 
{
    for (uint32_t i = 0; i < n_; ++i) {
        order_[i] = i;
    }
    
    std::stable_sort(order_.begin(), order_.end(),
                        [&graph_deg](const uint32_t v1, const uint32_t v2) {
                            return graph_deg[v1] > graph_deg[v2];
                        });
}

// Bi-Degree order
void USPCIndex::BiDegreeOrder(const std::vector<uint32_t>& graph_deg) 
{
    for (uint32_t i = 0; i < n_/2; ++i) {
        order_[i] = i;
    }
    
    std::stable_sort(order_.begin(), order_.end(),
                        [&graph_deg](const uint32_t v1, const uint32_t v2) {
                            return graph_deg[v1] > graph_deg[v2];
                        });

    std::vector<uint32_t> tmp_order(n_);
    for (uint32_t i = 0; i < n_/2; ++i) {
        tmp_order[i*2+1] = order_[i];
        tmp_order[i*2] = order_[i] + n_/2;
    }
    order_.resize(n_);
    order_ = tmp_order;
}

// print labels
void USPCQuery::printLabel(const uint32_t v, int flag) 
{
    std::string out_label;
    std::string in_label;
    std::string in_label_v;
    std::string out_label_vc;
    if (flag == 0) {
        out_label = std::to_string(v) + "_ori_out.txt";
        in_label = std::to_string(v) + "_ori_in.txt";
    } else if (flag == 1) {
        out_label = std::to_string(v) + "_bi_out.txt";
        in_label = std::to_string(v + n_/2) + "_bi_in.txt";
    } else if (flag == 2) {
        out_label = std::to_string(v) + "_incBefore_out.txt";
        in_label = std::to_string(v + n_/2) + "_incBefore_in.txt";
    } else if (flag == 3) {
        out_label = std::to_string(v) + "_dec_out.txt";
        in_label = std::to_string(v + n_/2) + "_dec_in.txt";
        in_label_v = std::to_string(v) + "_dec_in.txt";
        out_label_vc = std::to_string(v + n_/2) + "_dec_out.txt";
    }
    else {
        out_label = std::to_string(v) + "_inc_out.txt";
        in_label = std::to_string(v + n_/2) + "_inc_in.txt";
    }

    FILE* file1 = fopen(out_label.c_str(), "w");
    FILE* file2 = fopen(in_label.c_str(), "w");
    FILE* file3;
    FILE* file4;
    if (flag == 3) {
        file3 = fopen(in_label_v.c_str(), "w");
        file4 = fopen(out_label_vc.c_str(), "w");
        for (int i = 0; i < cL_[v].size(); ++i) {
            fprintf(file3,"%u %u %u %u\n",LEExtractV(cL_[v][i]),LEExtractD(cL_[v][i]),\
                LEExtractC(cL_[v][i]), rank_[LEExtractV(cL_[v][i])]);
        }
        fclose(file3);
        for (int i = 0; i < rev_cL_[v + n_/2].size(); ++i) {
            fprintf(file4,"%u %u %u %u\n",LEExtractV(rev_cL_[v + n_/2][i]),LEExtractD(rev_cL_[v + n_/2][i]),\
                LEExtractC(rev_cL_[v + n_/2][i]), rank_[LEExtractV(rev_cL_[v + n_/2][i])]);
        }
        fclose(file4);
    }

    for (int i = 0; i < rev_cL_[v].size(); ++i) {
        fprintf(file1,"%u %u %u %u\n",LEExtractV(rev_cL_[v][i]),LEExtractD(rev_cL_[v][i]),\
            LEExtractC(rev_cL_[v][i]), rank_[LEExtractV(rev_cL_[v][i])]);
    }
    fclose(file1);
    for (int i = 0; i < cL_[v + n_/2].size(); ++i) {
        fprintf(file2,"%u %u %u %u\n",LEExtractV(cL_[v + n_/2][i]),LEExtractD(cL_[v + n_/2][i]),\
            LEExtractC(cL_[v + n_/2][i]), rank_[LEExtractV(cL_[v + n_/2][i])]);
    }

    fclose(file2);
}

// clean redundant index
/*
uint32_t USPCQuery::clean_index(int fg, const uint32_t v){
    // in
    uint32_t cnt = 0;

    if (fg == 0) {
        // v's in label
        for(std::vector<LabelEntry>::iterator it    = cL_[v].begin(); it != cL_[v].end();)
        {
            uint32_t lv = LEExtractV(*(it));
            uint32_t ld = LEExtractD(*(it));

            if (lv == v) break;
            if (ld > Count(lv,v).first) {
                it = cL_[v].erase(it);
                for(std::vector<uint32_t>::iterator it1    = inv_in_[lv].begin(); it1 != inv_in_[lv].end();) {
                    if (*it1 == v) {
                        inv_in_[lv].erase(it1);
                        break;
                    } else {
                        it1++;
                    }
                }
                cnt++;
            }
            else {
                it++;
            }
        }

        // out-labels with v as hub
        for(std::vector<uint32_t>::iterator it1    = inv_out_[v].begin(); it1 != inv_out_[v].end();) {
            uint32_t pos = Find_Can_Noncan_DC(v,*it1,1);

            if (LEExtractD(rev_cL_[*it1][pos]) > Count(*it1,v).first) {
                std::vector<LabelEntry>::iterator it    = rev_cL_[*it1].begin();
                rev_cL_[*it1].erase(it + pos);
                it1 = inv_out_[v].erase(it1);
                cnt++;
            } else {
                it1++;
            }
        }
    } else {
        // v's out-labels
        for(std::vector<LabelEntry>::iterator it    = rev_cL_[v].begin(); it != rev_cL_[v].end();)
        {
            uint32_t lv = LEExtractV(*(it));
            uint32_t ld = LEExtractD(*(it));

            if (lv == v) break;

            if (ld > Count(v,lv).first) {
                it = rev_cL_[v].erase(it);

                for(std::vector<uint32_t>::iterator it1    = inv_out_[lv].begin(); it1 != inv_out_[lv].end();) {
                    if (*it1 == v) {
                        inv_out_[lv].erase(it1);
                        break;
                    } else {
                        it1++;
                    }
                }
                cnt++;
            }
            else {
                it++;
            }
        }

        // in-labels with v as hub
        for(std::vector<uint32_t>::iterator it1    = inv_in_[v].begin(); it1 != inv_in_[v].end();) {
            uint32_t pos = Find_Can_Noncan_DC(v,*it1,0);

            if (LEExtractD(cL_[*it1][pos]) > Count(v,*it1).first) {
                std::vector<LabelEntry>::iterator it    = cL_[*it1].begin();
                cL_[*it1].erase(it + pos);
                it1 = inv_in_[v].erase(it1);
                cnt++;
            } else {
                it1++;
            }
        }
    }

    return cnt;
}
*/

// check whether hub exist and update label
std::pair<bool,uint32_t> USPCQuery::CheckHubUpdate(uint32_t v, int in_out, uint32_t hub, uint32_t dis, uint32_t cnt, int fg) {
    uint32_t rmv_cnt = 0;

    // 0 for in, others for out
    if (in_out == 0) {
        for (uint32_t i = 0; i < cL_[v].size(); ++i) {
            uint32_t current_hub = LEExtractV(cL_[v][i]);

            // if found the hub
            if (current_hub == hub) {

                // if dis unchange, add cnt
                if (dis == LEExtractD(cL_[v][i])) {
                    uint32_t new_cnt = cnt + LEExtractC(cL_[v][i]);
                    cL_[v][i] = LEMerge(hub, dis, new_cnt);
                } else {
                    cL_[v][i] = LEMerge(hub, dis, cnt);
                    // if (fg == 1) rmv_cnt += clean_index(0,v);
                }

                return std::make_pair(true,rmv_cnt);
            }
        } // end for: each hub
    } else {
        for (uint32_t i = 0; i < rev_cL_[v].size(); ++i) {
            uint32_t current_hub = LEExtractV(rev_cL_[v][i]);

            // if found the hub
            if (current_hub == hub) {
                // if dis unchange, add cnt; else, update with new cnt
                if (dis == LEExtractD(rev_cL_[v][i])) {
                    uint32_t new_cnt = cnt + LEExtractC(rev_cL_[v][i]);
                    rev_cL_[v][i] = LEMerge(hub, dis, new_cnt);
                } else {
                    rev_cL_[v][i] = LEMerge(hub, dis, cnt);
                    // if (fg == 1) rmv_cnt += clean_index(1,v);
                }

                return std::make_pair(true,rmv_cnt);
            }
        } // end for each hub    
    }

    return std::make_pair(false,0);
}

// update index with edge insertion
std::tuple<uint32_t, uint32_t,uint32_t> USPCQuery::IncEdge(uint32_t s, uint32_t t, int fg) {
    G_[s].push_back(t);
    rev_G_[t].push_back(s);

    uint32_t cnt_add = 0;
    uint32_t cnt_all = 0;
    uint32_t rmv = 0;
    // find affected nodes, the hubs of union of in-label of s and out-label of t
    uint32_t i1 = 0, i2 = 0;
    std::set<uint32_t> aff_s;
    std::set<uint32_t> aff_t;
    std::vector<uint32_t> affect_node;

    while (i1 < cL_[s].size() && i2 < rev_cL_[t].size()) {
        if (rank_[LEExtractV(cL_[s][i1])] == rank_[LEExtractV(rev_cL_[t][i2])]) {
            affect_node.push_back(LEExtractV(cL_[s][i1]));
            aff_s.insert(LEExtractV(cL_[s][i1]));
            aff_t.insert(LEExtractV(cL_[s][i1]));
            i1++; i2++;
        } else if (rank_[LEExtractV(cL_[s][i1])] < rank_[LEExtractV(rev_cL_[t][i2])]) {
            affect_node.push_back(LEExtractV(cL_[s][i1]));
            aff_s.insert(LEExtractV(cL_[s][i1]));
            i1++;
        } else {
            affect_node.push_back(LEExtractV(rev_cL_[t][i2]));
            aff_t.insert(LEExtractV(rev_cL_[t][i2]));
            i2++;
        }
    }

    while (i1 < cL_[s].size()) {
        affect_node.push_back(LEExtractV(cL_[s][i1]));
        aff_s.insert(LEExtractV(cL_[s][i1]));
        i1++;
    }

    while (i2 < rev_cL_[t].size()) {
        affect_node.push_back(LEExtractV(rev_cL_[t][i2]));
        aff_t.insert(LEExtractV(rev_cL_[t][i2]));
        i2++;
    }

    // update process 
    std::vector<uint32_t> D(n_, UINT32_MAX);
    std::vector<uint32_t> C(n_, 0);

    // update for each hub 
    for (uint32_t i = 0; i < affect_node.size(); ++i) {
        uint32_t hub = affect_node[i];

        //*************************************** FORWARD **************************************
        // the dis and cnt from hub to s under old index (Only its own, can or non-can)
        if (aff_s.count(hub) > 0 && rank_[hub] < rank_[t]) {
            uint32_t dc_index = Find_Can_Noncan_DC(hub,s,0);
            D[t] = LEExtractD(cL_[s][dc_index]) + 1;
            C[t] = LEExtractC(cL_[s][dc_index]);
            std::vector<uint32_t> reset({t});
            std::queue<std::tuple<uint32_t,uint32_t,uint64_t>> Q;
            Q.push(std::make_tuple(t, D[t], C[t]));

            // BFS
            while (!Q.empty()) {
                std::tuple<uint32_t,uint32_t,uint64_t> cur_entry = Q.front();
                Q.pop();
                uint32_t v = std::get<0>(cur_entry);
                uint32_t D_ = std::get<1>(cur_entry);
                uint32_t C_ = std::get<2>(cur_entry);
                uint32_t d_old = Count(hub, v).first;

                if (D_ > d_old) continue;
                else {
                    std::pair<bool,uint32_t> check_cnt = CheckHubUpdate(v,0,hub,D_,C_,fg);
                    rmv += check_cnt.second;

                    if (!check_cnt.first) {
                        // add index if hub not exist, find the location
                        uint32_t find;
                        for (find = 0; find < cL_[v].size(); ++find) {
                            if (rank_[LEExtractV(cL_[v][find])] > rank_[hub])
                                break;
                        }

                        cL_[v].emplace(cL_[v].begin()+find, LEMerge(hub,D_,C_));
                        cnt_add++;
                        // inv_in_[hub].push_back(v);

                        // clean v's in-labels and those out-labels with v as hub
                        // if (fg == 1) rmv += clean_index(0,v);
                    }
                    cnt_all++;
                }
                
                // bfs next round
                for (const uint32_t w : G_[v]) {
                    if (rank_[hub] >= rank_[w] ) continue;

                    if (D[w] > D[v] + 1) {
                        D[w] = D[v] + 1; C[w] = C[v];
                        Q.push(std::make_tuple(w, D[w],C[w]));
                        reset.push_back(w);
                    } else if (D[w] == D[v] + 1) {
                        C[w] = C[v];
                        Q.push(std::make_tuple(w, D[w],C[w]));
                    }
                } // end for: bfs
            } // end while: Q

            for (const uint32_t v : reset) {
                D[v] = UINT32_MAX; C[v] = 0;
            }
        } // end if: aff_s

        //************************************ REVERSE ****************************************
        // the dis and cnt from t to hub under old index 
        if (aff_t.count(hub) > 0 && rank_[hub] < rank_[s]) {
            uint32_t rev_dc_index = Find_Can_Noncan_DC(hub,t,1);
            D[s] = LEExtractD(rev_cL_[t][rev_dc_index]) + 1;
            C[s] = LEExtractC(rev_cL_[t][rev_dc_index]);
            std::vector<uint32_t> rev_reset({s});
            std::queue<std::tuple<uint32_t,uint32_t,uint64_t>> rev_Q;
            rev_Q.push(std::make_tuple(s, D[s], C[s]));

            // BFS 
            while (!rev_Q.empty()) {
                std::tuple<uint32_t,uint32_t,uint64_t> cur_entry = rev_Q.front();
                rev_Q.pop();
                uint32_t v = std::get<0>(cur_entry);
                uint32_t D_ = std::get<1>(cur_entry);
                uint32_t C_ = std::get<2>(cur_entry);
                uint32_t d_old = Count(v, hub).first;

                if (D_ > d_old) continue;
                else {
                    std::pair<bool,uint32_t> check_cnt = CheckHubUpdate(v,1,hub,D_,C_,fg);
                    rmv += check_cnt.second;

                    if (!check_cnt.first) {
                        // add index if hub not exist, find the location
                        uint32_t find;
                        for (find = 0; find < rev_cL_[v].size(); ++find) {
                            if (rank_[LEExtractV(rev_cL_[v][find])] > rank_[hub]) 
                                break;
                        }

                        rev_cL_[v].emplace(rev_cL_[v].begin()+find, LEMerge(hub,D_,C_));
                        cnt_add++;
                        // inv_out_[hub].push_back(v);

                        // clean v's out-labels and those in-labels with v as hub
                        // if (fg == 1) rmv += clean_index(1,v);
                    }
                    cnt_all++;
                }
                
                // bfs next round 
                for (const uint32_t w : rev_G_[v]) {
                    if (rank_[hub] >= rank_[w]) continue;

                    if (D[w] > D[v] + 1)    {
                        D[w] = D[v] + 1; C[w] = C[v];
                        rev_Q.push(std::make_tuple(w, D[w],C[w]));
                        rev_reset.push_back(w);
                    } else if (D[w] == D[v] + 1) {
                        C[w] = C[v];
                        rev_Q.push(std::make_tuple(w, D[w],C[w]));
                    }
                } // end for: bfs
            } // end while: Q

            for (const uint32_t v : rev_reset) {
                D[v] = UINT32_MAX; C[v] = 0;
            }
        } //end if: aff_t hub
    } // end for: each hub

    return std::make_tuple(cnt_add, cnt_all-cnt_add, rmv);
}

// update index with edge deletion (affS,affT,rmv,add)
std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> USPCQuery::DecEdge(uint32_t s, uint32_t t) {

    // use set to store affS and affT
    std::set<uint32_t> affS_set;
    std::set<uint32_t> affT_set;

    // find AFF(s)
    std::queue<uint32_t> QS;
    std::set<uint32_t> visS;
    uint32_t cp_s = s+n_/2;
    QS.push(cp_s);
    visS.insert(cp_s);

    while (!QS.empty()) {
        uint32_t v = QS.front();
        QS.pop();

        if (Count_dis(v,cp_s) + 2 == Count_dis(v,t)) {
            affS_set.insert(v);

            for (auto nbrS: rev_G_[v]) {
                if (visS.count(nbrS+n_/2) < 1) {
                    QS.push(nbrS+n_/2);
                    visS.insert(nbrS+n_/2);
                }
            }
        }
    }

    // find AFF(t)
    std::queue<uint32_t> QT;
    std::set<uint32_t> visT;
    QT.push(t);
    visT.insert(t);

    while (!QT.empty()) {
        uint32_t v = QT.front();
        QT.pop();

        if (v == cp_s && Count_dis(s,cp_s) == Count_dis(t,cp_s) + 1) {
            visT.insert(v);
            affT_set.insert(v);
        }

        if (Count_dis(cp_s,v) == Count_dis(t,v) + 2) {
            affT_set.insert(v);

            uint32_t cp_t = v-n_/2;
            for (auto nbrT: G_[cp_t]) {
                if (visT.count(nbrT) < 1) {
                    QT.push(nbrT);
                    visT.insert(nbrT);
                }
            }
        }
    }

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ affS and affT are found ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    uint32_t cnt_rmv = 0;

    // remove labels from cL_ and rev_cL_
    for (auto rm_s: affS_set) {
        bool affect = false;
        for (auto rm_t: affT_set) {
            // do not remove own (h,0,1) label
            if (rm_s != rm_t) {

                if (rank_[rm_s] < rank_[rm_t]) {
                    uint32_t posIN = Find_Can_Noncan_DC(rm_s, rm_t, 0);
                    if (posIN != -1) {
                        // remove (rm_s,d,c) from Lin(rm_t); remove rm_t from inv_in(rm_s)
                        cL_[rm_t].erase(cL_[rm_t].begin() + posIN);
                        cL_[rm_t-n_/2].erase(cL_[rm_t-n_/2].begin() + posIN);
                        cnt_rmv++;
                    }
                } else {
                    uint32_t posOUT = Find_Can_Noncan_DC(rm_t, rm_s, 1);
                    if (posOUT != -1) {
                        // remove (rm_t,d,c) from Lout(rm_s); remove rm_s from inv_out(rm_t)
                        rev_cL_[rm_s].erase(rev_cL_[rm_s].begin() + posOUT);
                        rev_cL_[rm_s-n_/2].erase(rev_cL_[rm_s-n_/2].begin() + posOUT);
                        cnt_rmv++;
                    }
                }
            } else {
                // remove couple_in from couple_out's out label
                if (rev_cL_[rm_s-n_/2].size() >= 2) {
                    if (LEExtractV(rev_cL_[rm_s-n_/2][rev_cL_[rm_s-n_/2].size()-2]) == rm_t) {
                        auto last = rev_cL_[rm_s-n_/2].back();
                        rev_cL_[rm_s-n_/2].pop_back();
                        rev_cL_[rm_s-n_/2].pop_back();
                        rev_cL_[rm_s-n_/2].push_back(last);
                    }
                }
            }
        }
    }

    uint32_t rmv_edge = 0;
    for (rmv_edge = 0; rmv_edge < G_[s].size(); ++rmv_edge) {
        if (G_[s][rmv_edge] == t) break;
    }
    G_[s].erase(G_[s].begin()+rmv_edge);

    uint32_t rmv_edge_rev = 0;
    for (rmv_edge_rev = 0; rmv_edge_rev < rev_G_[t].size(); ++rmv_edge_rev) {
        if (rev_G_[t][rmv_edge_rev] == s) break;
    }
    rev_G_[t].erase(rev_G_[t].begin()+rmv_edge_rev);

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ remove possible redundant labels ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    std::set<uint32_t> all_aff_set; // all affected vertices in rank order
    all_aff_set.insert(affS_set.begin(),affS_set.end());
    all_aff_set.insert(affT_set.begin(),affT_set.end());

    std::vector<uint32_t> all_aff;
    all_aff.assign(all_aff_set.begin(),all_aff_set.end());

    std::vector<uint32_t> rank_all_aff(n_);
    for (uint32_t i = 0; i < all_aff.size(); ++i) {
        rank_all_aff[all_aff[i]] = rank_[all_aff[i]];
    }
    std::stable_sort(all_aff.begin(), all_aff.end(),
        [rank_all_aff](const uint32_t v1, const uint32_t v2) {
            return rank_all_aff[v1] < rank_all_aff[v2];
        });

    std::vector<uint32_t> dLu(n_,UINT32_MAX);
    // add new labels
    uint32_t add_cnt = 0;

    for (auto aff:all_aff) {

        // aff in affS, push in forward direction
        if (affS_set.find(aff) != affS_set.end()) {
            std::vector<uint32_t> reset({aff, aff-n_/2});
            std::vector<uint32_t> D(n_, UINT32_MAX);
            std::vector<uint32_t> C(n_, 0);
            D[aff] = 0; C[aff] = 1;
            D[aff-n_/2] = 1; C[aff-n_/2] = 1;

            std::queue<uint32_t> Q({aff-n_/2});

            for (auto e: rev_cL_[aff]) dLu[LEExtractV(e)] = LEExtractD(e);

            while (!Q.empty()) {
                uint32_t v = Q.front(); Q.pop();

                if ((v != aff || v != aff-n_/2) && (affT_set.find(v) != affT_set.end())){
                    const uint32_t dSoFar = Distance(dLu,cL_[v]);

                    if (dSoFar < D[v]) continue;
                    else {
                        // add to v's in-label
                        uint32_t i = 0;
                        int32_t low = 0; int32_t high = cL_[v].size() - 1; int32_t mid;

                        while (low <= high) {
                            mid = (low + high) / 2;
                            if (rank_[aff] >= rank_[LEExtractV(cL_[v][mid])]) low = mid + 1;
                            else high = mid - 1;
                        }
                        i = low > (cL_[v].size()-1) ? cL_[v].size()-1 : low;

                        cL_[v].insert(cL_[v].begin()+i,LEMerge(aff,D[v],C[v]));
                        
                        D[v-n_/2] = D[v] + 1;
                        C[v-n_/2] = C[v];

                        // add to v's couple's in-label
                        reset.push_back(v-n_/2);
                        cL_[v-n_/2].insert(cL_[v-n_/2].begin()+i,LEMerge(aff,D[v-n_/2],C[v-n_/2]));
                        v = v-n_/2;
                        add_cnt++;
                    }
                }

                // if meeting Vin and it is not in affT, move directly to its couple
                if (v >= n_/2) {
                    D[v-n_/2] = D[v]+1;
                    C[v-n_/2] = C[v];
                    v = v-n_/2;
                    reset.push_back(v);
                }

                for (auto nbr:G_[v]) {
                    if (rank_[nbr] <= rank_[aff]) continue;
                    if (D[nbr] == D[v] + 1) {
                        C[nbr] += C[v];
                    } else if (D[nbr] > D[v] + 1){
                        D[nbr] = D[v] + 1;
                        C[nbr] = C[v];
                        Q.push(nbr);
                        reset.push_back(nbr);
                    }
                }
            }

            for (auto v:reset) {
                D[v] = UINT32_MAX; C[v] = 0;
            }
            for (auto e:rev_cL_[aff]) {
                dLu[LEExtractV(e)] = UINT32_MAX;
            }
        } // if G

        // aff in affT, push in reverse direction
        if (affT_set.find(aff) != affT_set.end()) {
            std::vector<uint32_t> reset({aff});
            std::vector<uint32_t> D(n_, UINT32_MAX);
            std::vector<uint32_t> C(n_, 0);
            D[aff] = 0; C[aff] = 1;
            
            std::queue<uint32_t> Q({aff});

            for (auto e: cL_[aff]) dLu[LEExtractV(e)] = LEExtractD(e);

            while (!Q.empty()) {
                uint32_t v = Q.front(); Q.pop();

                if ((v != aff) && (affS_set.count(v+n_/2) > 0)) {

                    const uint32_t dSoFar = Distance(dLu, rev_cL_[v]);

                    if (dSoFar < D[v]) continue;
                    else {
                        // add to v's out-label
                        uint32_t i = 0;
                        int32_t low = 0; int32_t high = rev_cL_[v].size() - 1; int32_t mid;

                        while (low <= high) {
                            mid = (low + high) / 2;
                            if (rank_[aff] >= rank_[LEExtractV(rev_cL_[v][mid])]) low = mid + 1;
                            else high = mid - 1;
                        }
                        i = low > (rev_cL_[v].size()-1) ? rev_cL_[v].size()-1 : low;

                        rev_cL_[v].insert(rev_cL_[v].begin()+i,LEMerge(aff,D[v],C[v]));

                        if (v+n_/2 != aff) {
                            D[v+n_/2] = D[v] + 1;
                            C[v+n_/2] = C[v];
                            reset.push_back(v+n_/2);
                            rev_cL_[v+n_/2].insert(rev_cL_[v+n_/2].begin()+i,LEMerge(aff,D[v+n_/2],C[v+n_/2]));

                            v = v+n_/2;
                            add_cnt++;
                        } else {
                            continue;
                        }
                    }
                }

                if (v < n_/2) {
                    D[v+n_/2] = D[v]+1;
                    C[v+n_/2] = C[v];
                    v = v+n_/2;
                    reset.push_back(v);
                }

                for (auto nbr:rev_G_[v]) {
                    if (rank_[nbr] < rank_[aff]) continue;
                    if (D[nbr] == D[v] + 1) {
                        C[nbr] += C[v];
                    } else if (D[nbr] > D[v] + 1) {
                        D[nbr] = D[v] + 1;
                        C[nbr] = C[v];
                        Q.push(nbr);
                        reset.push_back(nbr);
                    }
                }
            }

            for (auto v:reset) {
                D[v] = UINT32_MAX; C[v] = 0;
            }
            for (auto e:cL_[aff]) {
                dLu[LEExtractV(e)] = UINT32_MAX;            
            }
        } // if G_bar
    }

    return std::make_tuple(affS_set.size(), affT_set.size(), cnt_rmv, add_cnt);
}
} // namespace spc
