# Optimize the parameter finite_element_solver.max(iter).
#
# This script is an example of how a solver parameter might be
# optimized. In order to find an optimal value for max(iter), we
# repeatedly perform a short simulation with a selection of values,
# and we record the computation time as a function of max(iter).
#
# The problem addressed by the simulation is a beefed-up version of the
# “Quick start guide”: a cylinder with height = diameter, with the
# initial magnetization at 45° from the axis, and no applied field.

from feellgood.meshMaker import Cylinder
import subprocess
import json

# Set this to True in order to get a very quick test that only checks
# that this script is working.
quick_test = False

if quick_test:
    # 3312-node mesh, very short simulated time.
    height = 40
    elt_size = 2.5
    final_time = 2e-12
    iter_values = (100, 315, 1000, 3150, 10000)
else:
    # 101889-node mesh, longer simulated time.
    height = 136
    elt_size = 2.5
    final_time = 1e-10

    # Selected values from the Renard series R10: geometric progression
    # 10**(n/10), conveniently rounded.
    iter_values = (
          100,  125,  160,  200,  250,  315,  400,  500,  630,  800,
         1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,
        10000)

# Create the mesh and save it in the file "cylinder.msh".
mesh = Cylinder(height / 2, height, elt_size, "surface", "volume")
mesh.make("cylinder.msh")

# Simulation settings.
#
# max(iter) plays a role in selecting the time step:
#
#  - if max(iter) is small, the solver will fail to converge on large
#    time steps, which will force feeLLGood to keep the time steps
#    suitably small
#
#  - if max(iter) is large, feeLLGood will use larger time steps, but
#    each one will take many iterations to be solved for.
#
# In order to highlight this role of max(iter), we make sure that no
# other parameter can affect the time step. We thus set the following
# parameters to very large values:
#
#  - outputs.evol_time_step
#  - time_integration.max(du)
#  - time_integration.max(dt)
settings = {
    "outputs": {
        "evol_time_step": final_time,
        "final_time": final_time,
        "mag_config_every": False
    },
    "mesh": {
        "filename": "cylinder.msh",
        "volume_regions": { "volume": { "alpha_LLG": 0.05 } },
        "surface_regions": { "surface": {} }
    },
    "initial_magnetization": [1, 0, 1],
    "demagnetizing_field_solver": { "nb_threads": 10 },
    "finite_element_solver": {
        "max(iter)": None  # this will be set in the main loop
    },
    "time_integration": {
        "max(du)": 2,
        "min(dt)": 1e-14,
        "max(dt)": 1e-10
    }
}

# Print a table of  (max(iter), time) pairs.
print("# max(iter)\ttime (s)")

# Loop over the selected max(iter) values.
for max_iter in iter_values:
    settings["finite_element_solver"]["max(iter)"] = max_iter

    # Run feeLLGood.
    process = subprocess.run(["feellgood", "-"],
        input=json.dumps(settings), capture_output=True, text=True)

    # The output of feeLLGood ends with lines that look like:
    #
    #    Computing time:
    #        total: ... s
    #        per time step: ... s
    #
    # Grab the number right after the string "total:".
    for line in process.stdout.splitlines():
        if "total: " in line:
            time = line.split()[1]
            print(f"{max_iter}\t{time}")
