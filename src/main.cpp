
#include "util/graphs.hpp"
#include "util/profile.hpp"
#include "util/isomorphism.hpp"
#include "solv/branching.hpp"
#include "math.h"
#include <algorithm>

#define num_edges(x) ((x*(x-1))/2)

string profile_names[] = { "A", "B", "C", "D", "E", "F", "G", "H" };
string internal_names[] = { "0", "1", "2", "3" , "4", "5", "6", "7", "8" };


/******** argument parsing *************/

typedef map<string, vector<string> >  arg_map;
const std::pair<string, int> _requires_params[] = {
  { "graph", 1 },
  { "enum", 0 },
  { "profile",  1 },
  { "all", 0 },
  { "-n", 1 }, // number of internal vertices
  { "-p", 1 } // number of profile vertices (heaps)
};

void usage(const char* progname, std::ostream& o){
  o << "usage: " << progname << " graph <file to read> [more opts]" << std::endl;
  o << "       " << progname << " profile <file to read> [more opts] "<< std::endl;
  o << "       " << progname << " all [more opts] "<< std::endl;
  o << "       " << progname << " enum [more opts] "<< std::endl;
  o << "more opts: " << " -n x\t <int>\t search for graphs with x internal vertices (default: 4, max: 8)"<< std::endl;
  o << "           " << " -p x\t <int>\t size of the profile (default: 4, max: 8)"<< std::endl;
  exit(1);
}

arg_map parse_args(int argc, char** argv){
  map<string, int>  requires_params(begin(_requires_params), end(_requires_params));
  arg_map result;
  int arg_ptr;
  arg_ptr = 1;
  while(arg_ptr < argc){
    const std::string arg(argv[arg_ptr++]);
    // if the argument is not registered in requires_args, then exit with usage
    if(requires_params.find(arg) == requires_params.end()) usage(argv[0], std::cerr);
    // if there are not enough parameters for this argument
    if(argc < arg_ptr + requires_params[arg]) usage(argv[0], std::cerr);
    // otherwise fill the argument map
    std::vector<string> params(requires_params[arg]);
    for(int i = 0; i < requires_params[arg]; ++i)
      params[i] = argv[arg_ptr++];
    result.insert(make_pair(arg, params));
  }
  return result;
}


/********* serious business ************/

using namespace vc;

// first 'profile_vertices' ids get the profile_names, rest gets the internal names
string get_name_by_id(const uint id, const uint profile_vertices){
  if(id < profile_vertices)
    return profile_names[id];
  else
    return internal_names[id - profile_vertices];
}

// construct the graph with given edges
graph get_graph(const vector<vector<bool> >& edges, const uint profile_vertices){
  const uint num_verts(edges.size());
  graph g;
  vector<vertex_p> vertices(num_verts);
  // vertices are A B C D v 0 1 2 3 ...
  for(uint id = 0; id < num_verts; ++id)
    vertices[id] = g.add_vertex_fast(get_name_by_id(id, profile_vertices));

  for(uint i1 = 0; i1 < num_verts; ++i1)
    for(uint i2 = max(i1, profile_vertices); i2 < num_verts; ++i2)
      if(edges[i1][i2]) g.add_edge_fast(vertices[i1], vertices[i2]);
  
  DEBUG3(cout << "done constructing new graph"<<endl);
  return g;
}

void add_profile_to_internal(graph& g, const vector<vector<bool> >& edges){
  const uint internal_vertices(edges.size());
  const uint profile_vertices(edges[0].size());

  vector<vertex_p> internals(internal_vertices);
  vector<vertex_p> profiles(profile_vertices);

  // get the internal vertices by name
  for(uint i = 0; i < internal_vertices; ++i)
    internals[i] = g.find_vertex_by_name(internal_names[i]);
  // get the profile vertices by adding them
  for(uint i = 0; i < profile_vertices; ++i)
    profiles[i] = g.add_vertex_fast(profile_names[i]);

  for(uint i1 = 0; i1 < internal_vertices; ++i1)
    for(uint i2 = 0; i2 < profile_vertices; ++i2)
      if(edges[i1][i2]) g.add_edge_fast(internals[i1], profiles[i2]);
}
// advance_to_next_graph to the next graph
bool advance_to_next_graph(vector<vector<bool> >& edges, const uint profile_vertices){
  const uint num_verts(edges.size());
  uint i1 = 0, i2 = 0;
  for(i1 = 0; i1 < num_verts; ++i1){
    for(i2 = max(i1 + 1, profile_vertices); i2 < num_verts; ++i2)
      if(!edges[i1][i2]) {
        edges[i1][i2] = true; // if a 0 was found, make it 1
        break;
      } else edges[i1][i2] = false; // all trailing 1's become 0
    if(i2 < num_verts) break;
  }
 
  DEBUG3(cout << "adj-matrix now:"<<endl;
      for(uint i = 0; i < num_verts; ++i) cout<<edges[i]<<endl; );
 
  // if all digits are 1, then return false, there is no successor
  return (i2 < num_verts);
}

// advance_to_next_graph to the next bipartite graph
bool advance_to_next_bipartite_graph(vector<vector<bool> >& edges){
  const uint part1(edges.size());
  const uint part2(edges[0].size());
  uint i1 = 0, i2 = 0;
  for(i1 = 0; i1 < part1; ++i1){
    for(i2 = 0; i2 < part2; ++i2)
      if(!edges[i1][i2]) {
        edges[i1][i2] = true; // if a 0 was found, make it 1
        break;
      } else edges[i1][i2] = false; // all trailing 1's become 0
    if(i2 < part2) break;
  }
 
  DEBUG2(cout << "adj-matrix now:"<<endl;
      for(uint i = 0; i < part1; ++i) cout<<edges[i]<<endl; );
 
  // if all digits are 1, then return false, there is no successor
  return (i2 != part2);
}

bool advance_to_next_border(vector<bool>& border){
  uint i;
  for(i = 0; i < border.size(); ++i)
    if(border[i]) border[i] = false; else {
      border[i] = true;
      break;
    }
  return (i < border.size());
}


uint vc_counter = 0;

// check if the profile of g matches (+/- offset) the p
bool profile_equal(const graph& g, const profile_t& p, const uint profile_vertices, const uint vc_num){
  // profile border: 1 = 'all neighbors are in the VC'
  vector<bool> profile_border(profile_vertices);
  // get the offset using the vc_num of g
  int offset = (int)vc_num - p.back();
  uint index = 0;

  DEBUG3(cout << "computing profile"<<endl);
  ++vc_counter;
  if(vc_counter % 500000 == 0) DEBUG5(cerr<<"crunched "<<vc_counter/1000<<"k graphs"<<endl);

  do{
    DEBUG2(cout << "profile containment in VC: "<<profile_border<<endl);
    graph gprime(g);
    solution_t s;
    for(uint i = 0; i < profile_vertices; ++i){
      // get i'th vertex (it's a meta-vertex)
      vertex_p X = gprime.find_vertex_by_name(get_name_by_id(i, profile_vertices));
      // if all of X are in the VC, delete v
      if(profile_border[i]) gprime.delete_vertex(X); else // else select all of N(X)
        for(edge_p e = X->adj_list.begin(); e != X->adj_list.end();){
          vertex_p u(e->head);
          ++e;
          select_vertex(gprime, u, s);
        }
    }
    // solve the rest of g
    s += run_branching_algo(gprime);
    DEBUG2(cout << "got size-"<<s.size()<<" solution: " << s<< endl);
    // if the solution size (offset by 'offset') does not match the profile, return failure
    if(s.size() != p[index++] + offset) return false;
  } while(advance_to_next_border(profile_border));
  return true;
}



profile_t get_profile(const graph& g, const uint profile_vertices){
  profile_t result(pow(2, profile_vertices));
  // profile border: 1 = 'all neighbors are in the VC'
  vector<bool> profile_border(profile_vertices);
  uint index = 0;

  DEBUG3(cout << "computing profile"<<endl);
  ++vc_counter;
  if(vc_counter % 100000 == 0) DEBUG5(cerr<<"crunched "<<vc_counter/1000<<"k graphs"<<endl);

  do{
    DEBUG2(cout << "profile containment in VC: "<<profile_border<<endl);
    graph gprime(g);
    solution_t s;
    for(uint i = 0; i < profile_vertices; ++i){
      // get i'th vertex
      vertex_p v = gprime.find_vertex_by_name(get_name_by_id(i, profile_vertices));
      // if not all of v's neighbors are in the VC, use v (but don't put it into a solution)
      if(profile_border[i]) gprime.delete_vertex(v); else // else select all
        for(edge_p e = v->adj_list.begin(); e != v->adj_list.end();){
          vertex_p u(e->head);
          ++e;
          select_vertex(gprime, u, s);
        }
    }
    s += run_branching_algo(gprime);
    DEBUG2(cout << "got size-"<<s.size()<<" solution: " << s<< endl);
    // save the optimal solution size in result[index]
    result[index++] = s.size();
  } while(advance_to_next_border(profile_border));
  DEBUG3(cout << "profile for graph "<<g<<":"<<endl);
  DEBUG3(cout << result<<endl);
  return result;
}

void output_all_profiles(const uint internal_vertices, const uint profile_vertices){
  unordered_map<profile_t, list<graph>, profile_hasher > equiv_class;
  // for each graph with n vertices, get its profile (that it, 16 solution sizes, depending on whether the neighbors of the first 4 vertices are selected or not)

  // forbit edges between A, B, C, D
  uint num_verts = internal_vertices + profile_vertices;
  vector<vector<bool> > edges(num_verts, vector<bool>(num_verts)); // wastes space, but simplifies the program
  for(uint i1 = 0; i1 < num_verts; ++i1)
    for(uint i2 = max(i1, profile_vertices); i2 < num_verts; ++i2) edges[i1][i2] = false;
  

  DEBUG2(cout << "done initializing edges"<<endl);
  do {
    // get the graph based on 'edges'
    graph g(get_graph(edges, profile_vertices));
    DEBUG3(cout << "got new graph"<< endl);
    DEBUG1(cout << "created graph "<< g<< endl);
    // get the profile of 'g'
    profile_t p(get_profile(g, profile_vertices));
    // add 'g' to the equivalence class of this profile
    equiv_class[p].push_back(g);
  } while(advance_to_next_graph(edges, profile_vertices));

  // output the equivalence classes
  cout << "EQUIVALENCE CLASSES:"<<endl;
  for(auto m = equiv_class.begin(); m != equiv_class.end(); ++m){
    cout << "================ "<< m->first << " ============================= "<<endl;
    // go through the list of graphs
    for(auto l = m->second.begin(); l != m->second.end(); ++l)
      l->print_edges(cout);
  }

}


// generate graphs that might have profile 'target' from an internal graph by adding profile vertices
void generate_and_apply(const graph& internal,
                        const profile_t& target,
                        const uint profile_vertices,
                        const int offset,
                        list<pair<vertex_p, vertex_p> >& allowed_edges,
                        void* apply(list<graph>&, const graph&, const profile_t, const uint)){
}

void print_if_equal(list<graph>& eq_class,
                    const graph& g,
                    const profile_t& p,
                    const uint profile_vertices,
                    const uint vc_num){
  if(profile_equal(g, p, profile_vertices, vc_num)){
    DEBUG5(cerr<<"found "; g.print_edges(cerr));
    eq_class.push_back(g);
  }
}



// add all graphs of the equivalence class of the target profile, agreeing on a fixed internal graph
void equiv_class_fixed_internal(const graph& g, 
                                const profile_t& target,
                                const uint profile_vertices,
                                list<graph>& equiv_class,
                                const uint vc_num){
  const uint internal_vertices(g.vertices.size());
  // forbit edges between A, B, C, D
  vector<vector<bool> > edges(internal_vertices, vector<bool>(profile_vertices)); // wastes space, but simplifies the program
  for(uint i1 = 0; i1 < internal_vertices; ++i1)
    for(uint i2 = 0; i2 < profile_vertices; ++i2)
      edges[i1][i2] = false;

  DEBUG5(cerr << "internal graph: "<<endl; g.print_edges(cerr););
  do {
    // get the graph based on 'edges'
    graph gprime(g);
    add_profile_to_internal(gprime, edges);
    DEBUG1(cout << "added profile vertices, finished graph is: "<< gprime<< endl);
  
    print_if_equal(equiv_class, gprime, target, profile_vertices, vc_num);
  } while(advance_to_next_bipartite_graph(edges));
}


bool push_back_if_not_isomorphic(list<graph>& created_graphs, graph& g){
  DEBUG3(cout << "checking "<<g<<" against "<<created_graphs.size() << " saved graphs..."<<endl);
  for(list<graph>::iterator i = created_graphs.begin(); i != created_graphs.end(); ++i)
    if(isomorphic(g, *i)) return false;
  DEBUG3(cout << "passed"<<endl);
  // g is not isomorphic to anything in the list, so add it
  created_graphs.push_back(g);
  return true;
}

void output_equivalence_class(const profile_t& target, const uint internal_vertices, const uint profile_vertices){
  list<graph> equiv_class;
  const uint last_profile_entry(target[pow(2, profile_vertices) - 1]);

  DEBUG3(cout << "generating all "<<internal_vertices<<"-vertex graphs of VC num "<<last_profile_entry<<endl);
  // STEP 1. generate all internal graph whose vertex cover is at most the profile's last entry
  // (when all profile vertices are in)
  vector<vector<bool> > internal_edges(internal_vertices, vector<bool>(internal_vertices)); // wastes space, but simplifies the program
  list<graph> created_graphs;

  for(uint i1 = 0; i1 < internal_vertices; ++i1)
    for(uint i2 = i1; i2 < internal_vertices; ++i2)
      internal_edges[i1][i2] = false;  
  do {
    // get the graph based on 'internal_edges', no profile
    graph g(get_graph(internal_edges, 0));
    DEBUG1(cout << "created graph "<< g<< endl);
    // compute its vertex cover s
    solution_t s(run_branching_algo(g));
    if(s.size() <= last_profile_entry){ 
      graph gprime(get_graph(internal_edges, 0));
      if(push_back_if_not_isomorphic(created_graphs, gprime))
        equiv_class_fixed_internal(created_graphs.back(), target, profile_vertices, equiv_class, s.size());
    }
  } while(advance_to_next_graph(internal_edges, 0));

  // output the equivalence classes
  cout << "EQUIVALENCE CLASSES:"<<endl;
  cout << "================ "<< target << " ============================= "<<endl;
  // go through the list of graphs
  for(auto l = equiv_class.begin(); l != equiv_class.end(); ++l)
    l->print_edges(cout);
}

void output_all_non_isomorphic(const uint num_verts){
  DEBUG3(cout << "generating all "<<num_verts<<"-vertex graphs"<<endl);
  // STEP 1. generate all internal graph whose vertex cover is at most the profile's last entry
  // (when all profile vertices are in)
  vector<vector<bool> > edges(num_verts, vector<bool>(num_verts)); // wastes space, but simplifies the program
  list<graph> created_graphs;

  for(uint i1 = 0; i1 < num_verts; ++i1)
    for(uint i2 = i1; i2 < num_verts; ++i2)
      edges[i1][i2] = false;  
  do {
    // get the graph based on 'internal_edges', no profile
    graph g(get_graph(edges, 0));
    // compute its vertex cover
    push_back_if_not_isomorphic(created_graphs, g);
  } while(advance_to_next_graph(edges, 0));
 
  // output the equivalence classes
  cout << "non-isomorphic "<<num_verts<<"-vertex graphs:"<<endl;
  // go through the list of graphs
  for(auto l = created_graphs.begin(); l != created_graphs.end(); ++l)
    l->print_edges(cout);
}

uint internal_vertices = 4;
uint profile_vertices = 4;

int main(int argc, char** argv){
  // parse arguments
  arg_map arguments(parse_args(argc, argv));
  // act on arguments
  // first: parse options
  if(arguments.find("-n") != arguments.end()) internal_vertices = atoi(arguments["-n"][0].c_str());
  if(arguments.find("-p") != arguments.end()) profile_vertices = atoi(arguments["-p"][0].c_str());
  // then: parse actions
  if(arguments.find("graph") != arguments.end()){
    // read profile from graph and output equivalent graphs
    graph g;
    g.read_from_file(arguments["graph"][0].c_str());
    profile_t target(get_profile(g, profile_vertices));
    DEBUG5(cout << "found profile: "<<target<<" now looking for equivalent profiles..."<<endl);
    output_equivalence_class(target, internal_vertices, profile_vertices);

  } else if(arguments.find("profile") != arguments.end()){
// TODO: implement me
//    profile_t target_profile = read_profile_from_file(arguments["profile"][0].c_str());
//    output_equivalence_class(target_profile);

  } else if(arguments.find("all") != arguments.end()){
    output_all_profiles(internal_vertices, profile_vertices);
  } else if(arguments.find("enum") != arguments.end()){
    output_all_non_isomorphic(internal_vertices);
  } else usage(argv[0], std::cerr);
}
