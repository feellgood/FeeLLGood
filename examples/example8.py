#!/usr/bin/env python3
# Initial magnetization depending on the volume region.
#
# The sample is a bilayer: two magnetic layers, "bottom" and "top", stacked along z and meshed
# together (they share the nodes of their interface). The initial magnetization is given by a
# JavaScript function of the position (x, y, z) and of the list of the names of the volume
# regions the node belongs to:
#  - nodes of the top layer only:    along x
#  - nodes of the bottom layer only: at 45 degrees in the plane (x + y)
#  - nodes of the interface (in both regions): in between (22.5 degrees)
#
# The function must take four parameters to receive the region names; it returns the
# magnetization [mx, my, mz], which is normalized by feeLLGood. x, y, z are in metres.
# Requires the gmsh python module (pip install gmsh).

import json
import gmsh

# Names of generated files.
file_basename = "bilayer"
mesh_filename = file_basename + ".msh"
json_filename = file_basename + ".json"

# Dimensions in nm.
length = 60
width = 30
thickness = {"bottom": 5, "top": 5}
mesh_size = 3

# Mesh: two boxes fragmented into a single conforming mesh.
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.model.add(file_basename)
bottom = gmsh.model.occ.addBox(0, 0, 0, length, width, thickness["bottom"])
top = gmsh.model.occ.addBox(0, 0, thickness["bottom"], length, width, thickness["top"])
gmsh.model.occ.fragment([(3, bottom)], [(3, top)])
gmsh.model.occ.synchronize()

# Physical volumes, identified by the height of their center.
for dim, tag in gmsh.model.getEntities(3):
    z_center = gmsh.model.occ.getCenterOfMass(dim, tag)[2]
    name = "bottom" if z_center < thickness["bottom"] else "top"
    gmsh.model.addPhysicalGroup(3, [tag], 300 if name == "bottom" else 301)
    gmsh.model.setPhysicalName(3, 300 if name == "bottom" else 301, name)

gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
gmsh.model.mesh.generate(3)
gmsh.write(mesh_filename)
gmsh.finalize()

# Initial magnetization: JavaScript function of (x, y, z, regions), 'regions' being the array of
# the names of the volume regions sharing the node.
initial_magnetization = """function(x, y, z, regions) {
    if (regions.length > 1)                            // interface: node of both layers
        return [cos(PI/8), sin(PI/8), 0];
    if (regions.includes("top")) return [1, 0, 0];
    return [1, 1, 0];                                  // bottom layer (normalized by feeLLGood)
}"""

# Simulation settings.
settings = {
    "outputs": {
        "file_basename": file_basename,
        "evol_time_step": 1e-12,
        "final_time": 2e-11,
        "evol_columns": [
            "t",
            "bottom:<Mx>", "bottom:<My>", "bottom:<Mz>",
            "top:<Mx>", "top:<My>", "top:<Mz>",
            "E_tot"
        ],
        "mag_config_every": 5
    },
    "mesh": {
        "filename": mesh_filename,
        "length_unit": 1e-9,
        "volume_regions": {
            "bottom": {"Ms": 800e3, "Ae": 1e-11, "alpha_LLG": 0.1},
            "top": {"Ms": 400e3, "Ae": 1e-11, "alpha_LLG": 0.1}
        }
    },
    "initial_magnetization": initial_magnetization,
    "time_integration": {
        "min(dt)": 1e-16,
        "max(dt)": 1e-12
    }
}
with open(json_filename, "w") as outfile:
    json.dump(settings, outfile, indent=4)
    outfile.write("\n")

print(f"Prepared simulation of a bilayer: {mesh_filename} and {json_filename}")
print(f"Simulate with: feellgood {json_filename}")
