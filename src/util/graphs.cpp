#include "graphs.hpp"
#include <unordered_map>
#include <sstream>

namespace vc{

  uint hash_edge::operator()(const edge_p& e) const{
    uint id1(e->head->id);
    uint id2(e->get_tail()->id);
    return id1 * max(id1, id2) + id2;
  }

  // for all j in {0..sizeof(result)-1}:
  //    let V[j] be the set of vertices whose ID is j modulo sizeof(result)
  //    let D[j] be the set of last bits of the degrees of all vertices in V[j]
  //    set j'th most significant bit of result to the parity of D[j]
  uint hash_graph::operator()(const graph& g) const{
    uint result = 0;
    for(auto v = g.vertices.begin(); v != g.vertices.end(); ++v)
      result ^= ( (v->degree() & 1) << (v->id % (8*sizeof(result))) );
  // since the sum over all degrees is always even, the parity of result is always even
  // to compensate for the lost image space, shift it to the left and set the least significant
  // bit of result to the parity of the vertex set
    result = ( result << 1 ) | ( g.vertices.size() & 1 );
    return result;
  }


    
  bool graph::operator==(const graph& g) const{
    // make sure *this and g have the same # of vertices
    if(vertices.size() != g.vertices.size()) return false;
    auto v(vertices.begin());
    auto vprime(g.vertices.begin());
    while(v != vertices.end()){
      if(v->id != vprime->id) return false;
      if(v->degree() != vprime->degree()) return false;
      // test the sets of incident edges by the ids of the heads
      set<uint> el;
      for(auto e = v->adj_list.begin(); e != v->adj_list.end(); ++v)
        el.insert(e->head->id);
      for(auto e = v->adj_list.begin(); e != v->adj_list.end(); ++v)
        if(el.find(e->head->id) == el.end()) return false;

      ++v; ++vprime;
    }
    return true;
  }


  graph::graph(const graph& g, unordered_map<uint, vertex_p>* id_to_vertex):
    // copy graph infos
    current_id(-1)
  {
    DEBUG5(cout << "copy constructing a new graph with "<<g.vertices.size()<<" vertices"<<endl);
    add_disjointly(g, id_to_vertex);
  }

  // copy constructor - NOTE THAT trr_infos ARE NOT up to date for the copy
  graph::graph(const graph& g, edgelist* const el):
    // copy graph infos
    current_id(-1)
  {
    unordered_map<uint, vertex_p> id_to_vertex;
    add_disjointly(g, &id_to_vertex);
    // if we are also tasked with translating the edgelist el, then do so using id_to_vertex
    if(!el->empty()){
      DEBUG4(cout << "translating edgelist "<< *el << " using "<<id_to_vertex<<endl);
      edgelist::iterator eprime = el->begin();
      bool reached_end;
      do{
        // advance the iterator, since we're going to delete e (and add the corresponding edge to the front of el)
        edgelist::iterator e(eprime++);
        reached_end = (eprime == el->end());

        // push the translated edge to the front of el
        DEBUG4(cout << "translating ("<<id_to_vertex[(*e)->get_tail()->id]<<","<<id_to_vertex[(*e)->head->id]<<")"<<" with adjacency list "<<id_to_vertex[(*e)->head->id]->adj_list<<endl);
        el->push_front(find_edge(id_to_vertex[(*e)->get_tail()->id], id_to_vertex[(*e)->head->id]));
        // and delete e
        el->erase(e);
      } while(!reached_end);
      DEBUG4(cout << "done"<<endl);
    }
  }


    vertex::vertex(const vertex& v): id(v.id), name(v.name){}


  // return number of vertices present in the graph
  uint graph::num_vertices() const{
    return vertices.size();
  }

  // clear the graph (remove all vertices and edges)
  void graph::clear(){
    vertices.clear();
  }

  // find a vertex by specifying its id, return vertices.end() if its just not there
  vertex_p graph::find_vertex_by_id(const uint id){
    for(vertex_p i = vertices.begin(); i != vertices.end(); ++i)
      if(i->id == id) return i;
    return vertices.end();
  }
  vertex_p graph::find_vertex_by_name(const string& s){
    for(vertex_p i = vertices.begin(); i != vertices.end(); ++i)
      if(i->name == s) return i;
    return vertices.end();
  }

  // add a vertex with a brand new id to the graph and return a fresh iterator to it
  // for id, just take the last vertex in the list and invcease his id
  vertex_p graph::add_vertex_fast(){
    return add_vertex_fast(++current_id);
  }

  vertex_p graph::add_vertex_fast(const string& s){
    vertex_p v(add_vertex_fast());
    v->name = s;
    return v;
  }

  // add a vertex given as id to the graph and return a fresh iterator to it
  vertex_p graph::add_vertex_fast(const uint id){
    return vertices.insert(vertices.end(), vertex(id));
  }

  vertex_p graph::add_vertex_fast(const uint id, const string& s){
    vertex_p v(add_vertex_fast(id));
    v->name = s;
    return v;
  }


  // add a vertex given as id to the graph and return a fresh iterator to it
  // if it was there already, return an iterator to this one!
  vertex_p graph::add_vertex_secure(const uint id){
    vertex_p duplicate = find_vertex_by_id(id);
    if(duplicate == vertices.end())
      return add_vertex_fast(id);
    else
      return duplicate;
  }


  // add an edge to the graph - modify adjacency lists
  // this is the _fast_ variant: no check is done whether this edge already exists!
  edge_p graph::add_edge_fast(const vertex_p& u, const vertex_p& w){
    edge eu(u);
    edge ew(w);

    DEBUG4(cout << "adding edge "<<*u<<"-"<<*w<<endl);
    edge_p uadj_pos = u->adj_list.insert(u->adj_list.end(), ew);
    edge_p wadj_pos = w->adj_list.insert(w->adj_list.end(), eu);

    uadj_pos->head_adj_pos = wadj_pos;
    wadj_pos->head_adj_pos = uadj_pos;

    return uadj_pos;
  }

  edge_p graph::add_edge_fast(const vertex_p& u, const vertex_p& v, const edge_pc& copy_from){
    edge_p result(add_edge_fast(u, v));
    return result;
  }
  // add an edge to the graph - modify adjacency lists
  // this is the _secure_ variant: all sanity checks are performed
  edge_p graph::add_edge_secure(const vertex_p& u, const vertex_p& v){
    // first, disallow loops
    if(u == v) return u->adj_list.end();

    // next, disallow parallel edges
    if(adjacent(u, v)) return u->adj_list.end();

    // if the edge is not yet there, insert it
    return add_edge_fast(u,v);
  }

    // delete a vertex and all incident edge
    void graph::delete_vertex(const vertex_p& v){
      DEBUG5(cout << "deleting "<< *v << endl);
      // delete incident edges
      while(!v->adj_list.empty()) delete_edge(v->adj_list.begin());
      // and remove it from the vertex list
      vertices.erase(v);
    }
    void graph::delete_vertices(list<vertex_p>& vl){
      for(list<vertex_p>::iterator v = vl.begin(); v != vl.end(); ++v)
        delete_vertex(*v);
    }

    // delete an edge
    edge_p graph::delete_edge(const edge_p& e){
      const edge_p& mirror_e(e->head_adj_pos);
      
      DEBUG5(cout << "deleting edge "<< *e << endl);

      vertex_p w(e->head);
      vertex_p u(mirror_e->head);

      list<edge>& wadj(w->adj_list);
      list<edge>& uadj(u->adj_list);

      // perform the delete and return the next edge_p "in line"
      uadj.erase(mirror_e);
      return wadj.erase(e);
    }

    void graph::delete_edges(const edgelist& l){
      for(edgelist::const_iterator e = l.begin(); e != l.end(); ++e)
        delete_edge(*e);
    }

  // simple output,
  // Prints the edgelist of g (plus the number of vertices/edges in verbose mode)
  // g: a graph
  // P: Property map to access the information that should be print for the vertices 
  // out: output stream.
  void graph::write_to_stream(ostream& out, const bool verbose) const{
    
    if(verbose){
      out<<"number of vertices: "<<  num_vertices()  <<endl;
    }
  
    set<vertex> seen;
    for(list<vertex>::const_iterator v = vertices.begin(); v != vertices.end(); ++v){
      // mark v 'seen' so no edges involving v are printed later
      //std::pair<set<vertex>::iterator, bool> seen_entry = 
        seen.insert(*v);
      for(list<edge>::const_iterator e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
        // if v's current head has not been 'seen' yet, print the edge
        if(seen.find(*(e->head)) == seen.end()) out<< *e << endl;
    }
  } // end write_graph


  // simple input
  // reads the edgelist into g
  // in: input stream. 
  void graph::read_from_stream(istream& in){
    typedef unordered_map<string, vertex_p>::iterator map_iter;
    
    unordered_map<string, vertex_p> id2vertex;
    string name0,name1;
    pair<map_iter, bool> i,j;

    // clear the graph
    clear();

    while(in){
      // read the two endpoints
      in >> name0;
      in >> name1;

      // get their respective vertices (or vceate if they don't exist yet)
      i = id2vertex.insert(pair<string, vertex_p>(name0,vertex_p()));
      if(i.second) { // new vertex
        i.first->second = add_vertex_fast();
        i.first->second->name = name0;
      }
      j = id2vertex.insert(make_pair(name1,vertex_p()));
      if(j.second) { // new vertex
        j.first->second = add_vertex_fast();
        j.first->second->name = name1;
      }
      add_edge_secure(i.first->second,j.first->second);
    }
  } // end of read_graph


  // simple input
  // reads the edgelist into f
  // infile: input file namf
  void graph::read_from_file(const char* infile)
  {
    ifstream f(infile);
    read_from_stream(f);
  }


  // get all degree-one vertices
  list<vertex_p> graph::get_leaves(){
    list<vertex_p> tmp;
    for(vertex_p v = vertices.begin(); v != vertices.end(); ++v)
      if(v->degree() == 1) tmp.push_back(v);
    return tmp;
  }


  void graph::add_disjointly(const graph& Gfrom, unordered_map<uint, vertex_p>* id_to_vertex){
    if(Gfrom.vertices.empty()) return;
    const bool destroy_map(id_to_vertex == NULL);
    if(destroy_map) id_to_vertex = new unordered_map<uint, vertex_p>();
    // copy vertices
    for(vertex_pc x = Gfrom.vertices.begin(); x != Gfrom.vertices.end(); ++x){
      vertex_p y(add_vertex_fast(x->name));
      (*id_to_vertex)[x->id] = y;
      // copy edges involving v
      for(edge_pc e = x->adj_list.begin(); e != x->adj_list.end(); ++e){
        const unordered_map<uint, vertex_p>::const_iterator i(id_to_vertex->find(e->head->id));
        if(i != id_to_vertex->end()) add_edge_fast(i->second, y, e);
      }
    }
    if(destroy_map) delete id_to_vertex;
  }


};





