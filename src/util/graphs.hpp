#ifndef GRAPHS_HPP
#define GRAPHS_HPP

#include <cstdlib>
#include <set>
#include <list>
#include <vector>
#include <map>
#include <unordered_map>

#include <fstream>
#include <iostream>
#include <sstream>

#include "defs.hpp"


using namespace std;


 
namespace vc {

  // my graph will be a list of adjacency lists
  //  for each vertex, its adjacency list contains pointers to an edge
  //  edges contain pointers to both its endpoints _and_ to the appropriate positions in the respective adjacency lists

  class vertex;
  class edge;
  class graph;


  typedef list<vertex>::iterator vertex_p;
  typedef list<vertex>::const_iterator vertex_pc;
  typedef list<edge>::iterator edge_p;
  typedef list<edge>::const_iterator edge_pc;

  typedef list<vertex_p> vertexlist;
  typedef list<edge_p> edgelist;
  
  typedef vertexlist::iterator vertex_pp;
  typedef vertexlist::const_iterator vertex_ppc;
  typedef edgelist::iterator edge_pp;
  typedef edgelist::const_iterator edge_ppc;

  // hash an edge using the IDs of its vertices
  class hash_edge {
    public:
      uint operator()(const edge_p& e) const;
  };
  typedef unordered_map<edge_p, uint, hash_edge> weighted_edges;
  typedef pair<edge_p, uint> weighted_edge;

  // hash a graph using the ids of its vertices and their degrees
  class hash_graph {
    public:
      uint operator()(const graph& e) const;
  };

    // can't take a list of edges here because some reduction rules replace parts without knowing the exact edge to delete; they will just register "some edge incident to v"
  typedef list<string> solution_t;

};


namespace vc{
  class vertex {
    friend class edge;
    friend class graph;
  public:
    uint id;
    string name;
    // the adjacency list of the vertex
    list<edge> adj_list;

    /****************
     * constructors
     ****************/

    vertex(const uint new_id):id(new_id){}
    vertex(const vertex& v);
    
    // empty destructor for testing purposes, TODO: remove
    ~vertex(){};

    inline bool operator<(const vertex& v) const {return id < v.id;}
    inline bool operator==(const vertex& v) const {return id == v.id;}
    inline operator string() const{return name;}

    // This is preferred to is_on_cyclic_core, because if we only have a tree, the root doesn't have a parent!
    inline uint degree() const {
      return adj_list.size();
    }

    void init(const string& _name){
      name = _name;
    }
    void init_from(const vertex_p& v){
      init(v->name);
    }
    vertex& operator=(const vertex& v) {
      init(v.name);
      return *this;
    }

  };

  class vertex_hasher{
  public:
    uint operator()(const vertex_p& x) const{
      return x->id;
    }
  };
  typedef unordered_set<vertex_p, vertex_hasher> vertexset;


  class edge{
    friend class vertex;
    friend class graph;
  public:
    const vertex_p head;
    edge_p head_adj_pos;


    /******************
     * constructors
     ******************/
    edge(const vertex_p& _head):head(_head){}
    // empty destructor for testing purposes, TODO: remove
    ~edge(){};
        
    edge_p& get_reversed(){return head_adj_pos;}
    const edge_p& get_reversed() const{return head_adj_pos;}

    const vertex_p& get_head() const {return head;}
    const vertex_p& get_tail() const {return get_reversed()->head;}


    bool operator==(const vc::edge& e) const{
      return (head == e.head) && (head_adj_pos->head == e.head_adj_pos->head);
    }
    operator string() const {return (string)*(head_adj_pos->head) + "->" + (string)*(head);}
  };

  class edge_hasher{
  public:
    uint operator()(const edge_p& x) const{
      const vertex_hasher h;
      return h(x->head) + h(x->get_tail());
    }
  };
  typedef unordered_set<edge_p, edge_hasher> edgeset;

  class graph {
  private:
    void compute_bridges(edgelist& bridgelist, list<uint>& split_off_sizes);
  public:
    uint current_id;
    list<vertex> vertices;

    /****************************
     * constructors
     ***************************/
    graph():current_id(0),vertices(){}
    // initialize while translating the edge list el
    // note that the edgelist may change, but the pointer wont
    graph(const graph& g, edgelist * const el);
    graph(const graph& g, unordered_map<uint, vertex_p>* id_to_vertex = NULL);
    graph(const graph& g, edge_p& e){
      edgelist el; el.push_back(e);
      graph(g, &el);
    }

    /**************************
     * read-only informative functions
     **************************/

    uint num_vertices() const;

    // find a vertex by specifying its id, return vertices.end() if its just not there
    vertex_p find_vertex_by_id(const uint id);
    vertex_p find_vertex_by_name(const string& s);

    // simple output,
    // Prints the edgelist of g (plus the number of vertices/edges in verbose mode)
    // out: output stream.
    void write_to_stream(ostream& out, const bool verbose = true) const;

    void print_edges(ostream& out) const {
      for(vertex_pc v = vertices.begin(); v != vertices.end(); ++v)
        for(edge_pc e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
          if(e->head->id > v->id)
            out << "("<<v->name<<","<<e->head->name<<") ";
      out << endl;
    }

    // test whether two graphs are equal (including vertex id's and _the_order_in_the_vertex_list)
    // TODO: the second condition can be dropped if we keep a map:id->vertex in the graph instead of a list<vertex>
    bool operator==(const graph& g) const;

    /************************
     * graph modifications
     ************************/

    // clear the graph (remove all vertices and edges)
    void clear();

    // add a vertex given as id to the graph and return a fresh iterator to it
    vertex_p add_vertex_fast();
    vertex_p add_vertex_fast(const string& s);
    vertex_p add_vertex_fast(const uint id);
    vertex_p add_vertex_fast(const uint id, const string& s);
    vertex_p add_vertex_secure(const uint id);

    // add an edge to the graph - modify adjacency lists
    // this is the _fast_ variant: no check is done whether this edge already exists!
    edge_p add_edge_fast(const vertex_p& u, const vertex_p& v);
    edge_p add_edge_fast(const vertex_p& u, const vertex_p& v, const edge_pc& copy_from);

    // add an edge to the graph - modify adjacency lists
    // this is the _secure_ variant: all sanity checks are performed
    edge_p add_edge_secure(const vertex_p& u, const vertex_p& v);

    // delete a vertex and all incident edge
    void delete_vertex(const vertex_p& v);
    void delete_vertices(list<vertex_p>& vl);

    // delete an edge, return next edge_p in "current" adjacency list
    edge_p delete_edge(const edge_p& e);
    void delete_edges(const edgelist& l);

    void copy_graph(const vertex_p& v, graph& gto, unordered_map<uint,vertex_p>* id_to_vertex = NULL);

    // simple input
    // reads the edgelist into g
    // in: input stream.
    void read_from_stream(istream& in);

    // simple input
    // reads the edgelist into f
    // infile: input file namf
    void read_from_file(const char* infile);

    // add a graph to this one
    void add_disjointly(const graph& Gfrom, unordered_map<uint, vertex_p>* id_to_vertex = NULL);

//    void update_subtree_NH();
    // get all degree-one vertices
    list<vertex_p> get_leaves();
  };


  inline bool operator<(const vertex_p& v1, const vertex_p& v2){
    return v1->id < v2->id;
  }

  inline edge_p find_edge(const vertex_p& u, const vertex_p& v){
    for(edge_p e = u->adj_list.begin(); e != u->adj_list.end(); ++e)
      if(e->head == v) return e;
    return u->adj_list.end();
  }

  inline bool adjacent(const vertex_p& u, const vertex_p& v){
    return (find_edge(u, v) != u->adj_list.end());
  }


};

//ostream& operator<<(ostream& os, const vc::solution_t& s) {return os << s.i; }

inline ostream& operator<<(ostream& os, const vc::vertex& v) {return os << (string)v; }
inline ostream& operator<<(ostream& os, const vc::vertex_p& v) {return os << (const vc::vertex&)(*v); }
inline ostream& operator<<(ostream& os, const vc::vertex_pc& v) {return os << (const vc::vertex&)(*v); }

inline ostream& operator<<(ostream& os, const vc::edge& e) {return os << (string)e;}
inline ostream& operator<<(ostream& os, const vc::edge_p& e) {return os << (vc::edge)*e;}
inline ostream& operator<<(ostream& os, const vc::edge_pc& e) {return os << (vc::edge)*e;}

inline ostream& operator<<(ostream& os, const vc::graph& g){
  g.write_to_stream(os);
  return os;
}




#endif
