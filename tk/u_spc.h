#ifndef SPC_U_SPC_H_
#define SPC_U_SPC_H_

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include <set>

#include "macros.h"
#include "u_label.h"

namespace spc {
  class USPC {
   public:
    // ctors
    USPC() = default;
    USPC(const USPC&) = delete;
    USPC& operator=(const USPC&) = delete;

   protected:
    void OrderRank() {
      for (uint32_t i = 0; i < n_; ++i) {
        rank_[order_[i]] = i;
      }
    }
    // core data
    uint32_t n_;
    // graph and reverse graph
    Graph G_, rev_G_;
    // only cL_ and rev_cL_ are in use
    Label dL_, cL_, rev_dL_, rev_cL_;
    // for vertex v: its inlabel is can/noncan-L_[v], then it is a map, with key is hub's rank

    std::vector<uint32_t> deg_G_;
    // order and rank
    std::vector<uint32_t> order_;
    std::vector<uint32_t> rank_;
    // inverted index
    std::vector<std::vector<uint32_t>> inv_in_;
    std::vector<std::vector<uint32_t>> inv_out_;
    // index reduction
    std::vector<uint32_t> eqm_;
    std::vector<uint32_t> rev_eqm_;
  };

  class USPCIndex final: private USPC {
   public:
    enum class OrderScheme {
      kDegree,   // degree-based ordering,
      kBDegree,
      kDFS,
      kInvalid
    };
    // ctors
    USPCIndex() = default;
    USPCIndex(const USPCIndex&) = delete;
    USPCIndex& operator=(const USPCIndex&) = delete;
    // construct an index
    void BuildIndex_Ori(const Graph& const_graph, const Graph& rev_graph);
    void BuildIndex(const Graph& const_graph, const Graph& rev_graph);
    // merge index
    void MergeIndexAll();
    void MergeIndex();
    // write the index to disk
    uint64_t IndexWrite_Ori(const std::string& filename) const;
    uint64_t IndexWrite(const std::string& filename) const;
    // BFS count
    std::tuple<uint32_t, uint32_t> BFS_Count(const Graph& const_graph, uint32_t node);
    void set_os(const OrderScheme os) { os_ = os; }

   private:
    void find_dfs_order(const Graph& graph, uint32_t node, std::vector<bool>& check_nodes, std::vector<uint32_t>& orders);
    uint32_t Distance(const std::vector<uint32_t>& dLu,
                      const std::vector<LabelEntry>& dLv) const;
    // ordering functions
    void DegreeOrder(const std::vector<uint32_t>& graph_deg);
    void BiDegreeOrder(const std::vector<uint32_t>& graph_deg);
    void InvalidOrder(const std::vector<uint32_t>&) {
      ASSERT_INFO(false, "invalid ordering");
    }
    std::map<OrderScheme, void (USPCIndex::*)(const std::vector<uint32_t>&)> of_ = {
      {OrderScheme::kDegree,  &USPCIndex::DegreeOrder},
      {OrderScheme::kBDegree,  &USPCIndex::BiDegreeOrder},
      {OrderScheme::kInvalid, &USPCIndex::InvalidOrder}
    };
    OrderScheme os_ = OrderScheme::kInvalid;
  };

  class USPCQuery final: private USPC {
   public:
    USPCQuery() = default;
    USPCQuery(const USPCQuery&) = delete;
    USPCQuery& operator=(const USPCQuery&) = delete;
    // read the index from disk
    void IndexRead_Ori(const std::string& filename);
    void IndexRead(const std::string& filename);
    uint64_t IndexWrite(const std::string& filename) const;
    //uint64_t IndexWrite_q(const std::string& filename) const;
    // add edge
    std::tuple<uint32_t, uint32_t,uint32_t> IncEdge(uint32_t s, uint32_t t, int fg);
    // delete edge
    std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> DecEdge(uint32_t s, uint32_t t);
    // find the particular label entry L(v)-(hub,d,c)
    uint32_t Find_Can_Noncan_DC(uint32_t hub, uint32_t v, int in_out);
    // check hub exist or not and update
    std::pair<bool,uint32_t> CheckHubUpdate(uint32_t v, int in_out, uint32_t hub, uint32_t dis, uint32_t cnt, int fg);
    // remove redundant labels
    // uint32_t clean_index(int fg, const uint32_t v);

    // query
    std::uint32_t Count_dis(uint32_t v1, uint32_t v2) const;
    std::pair<uint32_t, uint64_t> Count(uint32_t v1, uint32_t v2) const;
    std::pair<uint32_t, uint64_t> Count_map(uint32_t v1, uint32_t v2);
    std::pair<uint32_t, uint64_t> Count_Bi_map(uint32_t node);
    std::tuple<uint32_t, uint32_t> Count_Ori_Cycle(uint32_t node) const;
    std::tuple<uint32_t, uint32_t> Count_Bigraph(uint32_t node) const;
    std::tuple<uint32_t, uint64_t, bool> Count_self(uint32_t v1, uint32_t v2) const ;
    std::tuple<uint32_t, uint64_t, std::set<uint32_t>,bool> Count_hub(uint32_t v1, uint32_t v2) const;
    void printLabel(const uint32_t v, int flag);
   private:
    uint32_t Distance(const std::vector<uint32_t>& dLu,
                      const std::vector<LabelEntry>& dLv) const;
    // used for IS-Counting
    mutable std::vector<uint32_t> dH_;
    mutable std::vector<uint64_t> cH_;
  };
}

#endif
