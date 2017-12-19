#ifndef ISOMORPHISM_HPP
#define ISOMORPHISM_HPP

#include "defs.hpp"
#include "graphs.hpp"

namespace vc{

  typedef map<uint, vector<vertex_p> >  degree_map_t;
  typedef map<vertex_p, vertex_p> isomorphism_t;


  // take a partial isomorphism (that's relating everyone with degree less than 'degree' and extend it)
  bool isomorphic_recursive(const degree_map_t& deg_map1,
                            const degree_map_t& deg_map2,
                            const degree_map_t::const_iterator current_deg,
                            isomorphism_t& partial_iso);

  
  bool isomorphic_brute_force(const degree_map_t& deg_map1, const degree_map_t& deg_map2);

  bool isomorphic(graph& g1, graph& g2);


}

#endif
