#ifndef PROFILE_H
#define PROFILE_H


#include "graphs.hpp"

namespace vc {
  
  typedef vector<uint> profile_t;

  class profile_hasher{
    public:
    uint operator()(const profile_t& p) const{
      uint result = 0;
      const uint profile_size(p.size());
      const uint bits_per_record = 8 * sizeof(uint) / profile_size;
      const uint lower_mask = (1 << bits_per_record) - 1;
      for(uint i = 0; i < profile_size; ++i)
        result |= ( p[i] && lower_mask  ) << (i * bits_per_record);
      return result;
    }
  };


  bool operator==(const profile_t& p1, const profile_t& p2){
    const uint profile_size(p1.size());
    if(profile_size != p2.size()) return false;
    // get the additive offset between p1 and p2
    int offset = p2[0] - p1[0];
    // check whether all entries of the profiles differ by exactly that offset
    for(uint i = 1; i < profile_size; ++i)
      if(offset + p1[i] != p2[i]) return false;
    return true;
  }

  int get_profile_offset(const profile_t& p1, const profile_t& p2){
    return p2.back() - p1.back();
  }

  profile_t read_profile_from_stream(istream& in) {
    return profile_t();
  }

  profile_t read_profile_from_file(const char* infile) {
    ifstream f(infile);
    return read_profile_from_stream(f);
  }

  uint border_to_profile_index(vector<bool>& border){
    const uint border_size = border.size();
    uint multiplicator = 1;
    uint result = 0;
    for(uint i = 0; i < border_size; ++i)
      result += (multiplicator << i) * border[i];
    return result;
  }
}

#endif
