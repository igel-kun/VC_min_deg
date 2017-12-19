#ifndef BRANCHING_HPP
#define BRANCHING_HPP


namespace vc{
  // run the complete branching recursively and return the number of operation it took
  solution_t run_branching_algo(graph& g);

  inline void select_vertex(graph& g, const vertex_p& v, solution_t& sol){
    sol += v->name;
    g.delete_vertex(v);
  }
}

#endif
