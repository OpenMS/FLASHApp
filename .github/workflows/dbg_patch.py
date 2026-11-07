# Temporary: adds prints for IHFA hits to FLASHExtenderAlgorithm.cpp (verification only)
import re, sys
p = sys.argv[1]
s = open(p).read()
key = "auto mod_count = getBlindModNumber_(*iter) + var_count; // blind and variable modifications so far\n"
assert key in s
s = s.replace(key, key + '              if (hit.getAccession().find("IHFA") != std::string::npos) std::cout << "DBGPATH pre " << hi.calculated_precursor_mass_ << " m " << m << " ni " << node_index << " pi " << pro_index << " mz " << t_node_spec[node_index].getMZ() << " shift " << mass_shift << " mods " << mod_count << std::endl;\n', 1)
m = re.search(r"\n( *)const int to = shiftRangeEnd_\([^\n]*\n", s)
assert m
ins = ('{0}if (hit.getAccession().find("IHFA") != std::string::npos) {{ std::cout << "DBGRANGE m " << m << " pre_node " << pre_node_index << " pre_pro " << pre_pro_index << " node " << node_index << " pro " << pro_index << " shift " << mass_shift << " to " << to << std::endl;'
       ' for (int q = pre_node_index + 1; q < node_index; q++) std::cout << "DBGMID " << q << " " << t_node_spec[q].getMZ() << " tol " << hi.tol_spec_map_[m][q].getIntensity() << std::endl; }}\n').format(m.group(1))
s = s[:m.end()] + ins + s[m.end():]
open(p, 'w').write(s)
print("patched")
# per path: the merged masses, sources and tolerances of each node, and whether it is left out
key2 = "            for (auto iter = best_path.rbegin(); iter != best_path.rend(); iter++)\n            {\n              if (reread["
if key2 in s:
    s = open(p).read()
    add = ('            if (hit.getAccession().find("IHFA") != std::string::npos) { const auto& rg = hi.node_range_map_[m]; const auto& sr = hi.node_source_map_[m];'
           ' std::cout << "DBGSRC m " << m << " ranges " << rg.size() << " spec " << t_node_spec.size() << " tols " << hi.tol_spec_map_[m].size() << " sources " << sr.size() << std::endl;'
           ' for (Size a = 0; a < best_path.size(); a++) { int na = getNodeIndex_(best_path[a], t_pro_masses.size()); std::cout << "DBGSRC node " << na << " mz " << std::to_string(t_node_spec[na].getMZ()) << " tol " << hi.tol_spec_map_[m][na].getIntensity() << " reread " << (int)reread[a] << " src";'
           ' if (na < (int)rg.size()) for (int k = rg[na].first; k <= rg[na].second; k++) std::cout << " " << sr[k]; std::cout << std::endl; } }\n')
    s = s.replace(key2, add + key2, 1)
    open(p, 'w').write(s)
    print("patched sources")
