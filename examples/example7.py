# Prepare a simulation for a pair on interacting nanoparticles:
#  - each nanoparticle is a tiny icosahedron (13-node mesh)
#  - the nanoparticles interact via their stray fields.

from math import sqrt
import numpy as np
from feellgood.settingsMaker import Settings
import json

# Names of generated files.
file_basename = "twinIco"
mesh_filename = file_basename + ".msh"
json_filename = file_basename + ".json"

# Single icosahedron mesh: nodes and faces.
phi = (1 + sqrt(5)) / 2
ico_nodes = np.array([
    [phi, 1, 0], [phi, -1, 0], [-phi, 1, 0], [-phi, -1, 0],
    [1, 0, phi], [1, 0, -phi], [-1, 0, phi], [-1, 0, -phi],
    [0, phi, 1], [0, phi, -1], [0, -phi, 1], [0, -phi, -1], [0, 0, 0]
])
ico_faces = np.array([
    [0, 8, 4], [0, 4, 1], [0, 1, 5], [0, 5, 9], [0, 9, 8],
    [8, 9, 2], [8, 2, 6], [8, 6, 4], [6, 10, 4], [6, 3, 10],
    [6, 2, 3], [10, 3, 11], [10, 11, 1], [10, 1, 4], [1, 11, 5],
    [11, 3, 7], [11, 7, 5], [3, 2, 7], [7, 9, 2], [5, 7, 9]
])

# Meshed pair of nanoparticles.
particles = [
    {
        "center": [0, 0, 0],
        "volume": { "name": "particle1", "tag": 300 },
        "surface":{ "name": "surface(particle1)", "tag": 200 }
    },
    {
        "center": [4, 0, 0],
        "volume": { "name": "particle2", "tag": 301 },
        "surface":{ "name": "surface(particle2)", "tag": 201 }
    }
]

# Return a list of mesh elements as a 2D array of numbers.
# Arguments:
#  - first_tag: tag of first element in the list
#  - extra_cols: extra constant columns as a 1D array of numbers
#  - offset_idx: offset to add to the node indices from ico_faces
def elements(first_tag, extra_cols, offset_idx):
    column_shape = (len(ico_faces), 1)
    tags = np.arange(len(ico_faces)).reshape(column_shape) + first_tag
    cst_cols = np.ones(column_shape, dtype=int) * extra_cols
    return np.hstack((tags, cst_cols, ico_faces + offset_idx))

# Prepare the mesh data.
nodes = np.vstack((
    ico_nodes + particles[0]["center"],
    ico_nodes + particles[1]["center"]
))
tetrahedrons = np.vstack((
    elements( 1, [4, 2, particles[0]["volume"]["tag"], 1, 13], 1),
    elements(21, [4, 2, particles[1]["volume"]["tag"], 2, 26], 14)
))
triangles = np.vstack((
    elements(41, [2, 2, particles[0]["surface"]["tag"], 1], 1),
    elements(61, [2, 2, particles[1]["surface"]["tag"], 2], 14)
))

# Write out the mesh.
meshFile = open(mesh_filename, "w")
meshFile.write(f"""$MeshFormat
2.2\t0\t8
$EndMeshFormat
$PhysicalNames
4
2 {particles[0]["surface"]["tag"]} "{particles[0]["surface"]["name"]}"
2 {particles[1]["surface"]["tag"]} "{particles[1]["surface"]["name"]}"
3 {particles[0]["volume"]["tag"]} "{particles[0]["volume"]["name"]}"
3 {particles[1]["volume"]["tag"]} "{particles[1]["volume"]["name"]}"
$EndPhysicalNames
$Nodes
{len(nodes)}
""")
for i in range(0, len(nodes)):
    meshFile.write("\t".join(map(str, [i+1]+list(nodes[i]))) + "\n")
meshFile.write(f"""$EndNodes
$Elements
{4*len(ico_faces)}
""")
for tetrahedron in tetrahedrons:
    meshFile.write("\t".join(map(str, tetrahedron)) + "\n")
for triangle in triangles:
    meshFile.write("\t".join(map(str, triangle)) + "\n")
meshFile.write("$EndElements\n")
meshFile.close()

# Simulation settings.
settings = {
    "outputs": {
        "evol_time_step": 1e-11,
        "final_time": 5e-09,
        "evol_columns": [
            "t",
            "particle1:<Mx>", "particle1:<My>", "particle1:<Mz>",
            "particle2:<Mx>", "particle2:<My>", "particle2:<Mz>",
            "E_tot"
        ],
        "mag_config_every": False
    },
    "mesh": {
        "filename": mesh_filename,
        "volume_regions": {
            particles[0]["volume"]["name"]: {
                "Ms": 800e3,
                "Ae": 1e-11,
                "alpha_LLG": 0.1
            },
            particles[1]["volume"]["name"]: {
                "Ms": 400e3,
                "Ae": 1e-11,
                "alpha_LLG": 0.1
            }
        },
        "surface_regions": {
            particles[0]["surface"]["name"]: {},
            particles[1]["surface"]["name"]: {}
        }
    },
    "initial_magnetization": [1, 0, 1],
    "time_integration": {
        "min(dt)": 1e-13,
        "max(dt)": 1e-11
    }
}
with open(json_filename, "w") as outfile:
    json.dump(settings, outfile, indent=4)
    outfile.write("\n")

print("Prepared simulation of two icosahedrons: "
        f"{mesh_filename} and {json_filename}")
print(f"Simulate with: feellgood {json_filename}")
