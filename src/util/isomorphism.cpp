#include "isomorphism.hpp"
#include <algorithm>


namespace vc{

  // take a partial isomorphism (that's relating everyone with degree less than 'degree' and extend it)
  bool isomorphic_recursive(const degree_map_t& deg_map1,
                            const degree_map_t& deg_map2,
                            const degree_map_t::const_iterator current_deg,
                            isomorphism_t& partial_iso){
    if(current_deg == deg_map1.end()){
      // if all degrees are included in partial_iso, check if it is an isomorphism
      for(auto i = deg_map1.begin(); i != deg_map1.end(); ++i)
        for(vertex_p v : i->second)
          for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
            if(!adjacent(partial_iso[v], partial_iso[e->head])) return false;
      // if all edges of deg_map1 are represented in deg_map2 and the degrees match up, then its isomorphic
      return true;
    } else {
      // else run thru the possible permutations for the current degree
      vector<vertex_p> vlist2(deg_map2.at(current_deg->first));
      vector<vertex_p> vlist1(current_deg->second);
      // get the first 
      std::sort(vlist2.begin(), vlist2.end());
      do {
        vector<vertex_p>::iterator i = vlist1.begin();
        vector<vertex_p>::iterator j = vlist2.begin();
        // assign one of the permutations to the current degree
        while(i != vlist1.end()){
          partial_iso[*i] = *j;
          ++i; ++j;
        }
        // recurse for the next degree
        degree_map_t::const_iterator next(current_deg); ++next;
        if(isomorphic_recursive(deg_map1, deg_map2, next, partial_iso))
          return true; // if an isomorphism is found return success
      } while(std::next_permutation(vlist2.begin(), vlist2.end()));
      return false;
    }
  }

  bool isomorphic_brute_force(const degree_map_t& deg_map1, const degree_map_t& deg_map2){
    isomorphism_t iso;
    return isomorphic_recursive(deg_map1, deg_map2, deg_map1.begin(), iso);
  }

  bool isomorphic(graph& g1, graph& g2){
    // check the vertex numbers
    if(g1.vertices.size() != g2.vertices.size()) return false;

    // get the degree maps to speed up isomorphism checking
    degree_map_t deg_map1, deg_map2;
    for(vertex_p v = g1.vertices.begin(); v != g1.vertices.end(); ++v)
      deg_map1[v->degree()].push_back(v);
    for(vertex_p v = g2.vertices.begin(); v != g2.vertices.end(); ++v)
      deg_map2[v->degree()].push_back(v);

    // check the degree sequences
    for(auto deg_pair = deg_map1.begin(); deg_pair != deg_map1.end(); ++deg_pair)
      if(deg_map2[deg_pair->first].size() != deg_pair->second.size())
        return false; // number of vertices of degree deg_pair->first in g1 and g2 differ

    DEBUG3(cout << g2 << " passed easy tests, going on to brute force"<<endl);

    bool result = isomorphic_brute_force(deg_map1, deg_map2);
    DEBUG3(cout << "isomorphic? "<<result<<endl);
    return result;
  }

}
