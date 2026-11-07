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
