# VC_min_deg
code to discover reduction rules for VC on small-degree-graphs

# input (see example input file)
### "graph" mode
input files should be lists of pairs of vertex names, where names "A" -- "H" are interpreted as "border sets"

### "profile" mode
not yet implemented

# options
-p XP -- number of profile vertices (default: 4, max: 8)

-n XN -- number of internal vertices (default: 4, max: 8)

# output (for debuglevel 0)
### "graph" mode
the output is 2 lines of header (including the profile of the input graph)
followed by a list of graphs with equivalent profile and NX internal nodes, one per line, with edges in format "(X,Y)"
### "profile" mode
not yet implemented
### "all" mode
for each of the 2^XP profiles with XP profile nodes, list all graphs with XN internal nodes that have this profile
### "enum" mode
enumerate all graphs with XN nodes
