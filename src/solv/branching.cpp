#include "../util/defs.hpp"
#include "../util/graphs.hpp"
#include "branching.hpp"

#include <algorithm> // for sort
#include <unordered_map>
#include <unordered_set>


namespace vc{


  vertex_p find_max_deg_vertex(graph& g){
    vertex_p result = g.vertices.begin();
    for(vertex_p v = g.vertices.begin(); v != g.vertices.end(); ++v)
      if(v->degree() > result->degree()) result = v;
    return result;
  }

  vertex_p find_min_deg_vertex(graph& g){
    vertex_p result = g.vertices.begin();
    for(vertex_p v = g.vertices.begin(); v != g.vertices.end(); ++v)
      if(v->degree() < result->degree()) result = v;
    return result;
  }

  inline void deg0_reduct(graph& g, vertex_p& v, solution_t& s){
    g.delete_vertex(v);
  }

  inline void deg1_reduct(graph& g, vertex_p& v, solution_t& s){
    const vertex_p neighbor(v->adj_list.front().head);
    select_vertex(g, neighbor, s);
  }


  inline void deg2_reduct(graph& g, vertex_p& v, solution_t& s){
    const vertex_p n1(v->adj_list.front().head);
    const vertex_p n2(v->adj_list.back().head);
    // that's all we needed from v
    g.delete_vertex(v);
    // if the two neighbors are adjacent, then take them (triangle->needs at least 2)
    if(find_edge(n1, n2) != n1->adj_list.end()){
      select_vertex(g, n1, s);
      select_vertex(g, n2, s);
    } else {
      // otherwise, contract both edges and decrease k by 1
      // 1. hang all neighbors of n1 to n2
      for(edge_p e = n1->adj_list.begin(); e != n1->adj_list.end(); ++e)
        g.add_edge_secure(n2, e->head);
      // 2. delete n1
      g.delete_vertex(n1);
      // 3. increase k
      s += ( n1->name + "/" + n2->name );
    }

  }

  solution_t run_branching_algo(graph& g){
    DEBUG2(cout << "running branching for graph with vertices: "<<g.vertices<<endl);
    if(g.vertices.size() <= 1) return solution_t();
    if(g.vertices.size() == 2){
      if(g.vertices.front().adj_list.empty()) return solution_t();
        else return solution_t(1, g.vertices.front().name);
    }
    vertex_p min_deg = find_min_deg_vertex(g);
    DEBUG2(cout << "min degree vertex: "<<min_deg<<endl);
    if(min_deg->degree() > 2){
      // min-deg > 2
      solution_t s1, s2;

      const vertex_p max_deg(find_max_deg_vertex(g));
      unordered_map<uint, vertex_p> id_to_vertex;
      graph gprime(g, &id_to_vertex);
      const vertex_p max_deg_prime(id_to_vertex[max_deg->id]);

      // either take him...
      select_vertex(g, max_deg_prime, s1);
      s1 += run_branching_algo(gprime);
      DEBUG2(cout << " selecting "<<s1.front()<<" yielded size-"<<s1.size()<<" solution "<<s1<<endl);
      // or take all his neighbors
      for(edge_p e = max_deg->adj_list.begin(); e != max_deg->adj_list.end();){
        vertex_p v = e->head;
        ++e;
        select_vertex(g, v, s2);
      }
      s2 += run_branching_algo(g);
      // return the smaller solution
      if(s1.size() < s2.size()) return s1; else return s2;
    } else {
      solution_t s;
      // min-deg < 2
      switch(min_deg->degree()){
        case 0: // degree-0, just delete it
          deg0_reduct(g, min_deg, s); break;
        case 1: // degree-1, take its neighbor
          deg1_reduct(g, min_deg, s); break;
        case 2: // degree-2,
          deg2_reduct(g, min_deg, s); break;
      }
      s += run_branching_algo(g);
      return s;
    }

  }



}; // end namespace


