import os
import sys
import json
import subprocess
from math import sqrt, pi, atan2, remainder

# Main simulation parameters.
A = 4e-3                 # field amplitude [T]
start_frequency = 5.3e9  # [Hz]
stop_frequency = 7.3e9   # [Hz]
nb_steps_frequency = 21

# Output files.
output_directory = "ferro_resonance"
output_path = f"{output_directory}/resonance.tsv"

# Using many threads is not useful with such a small problem.
thread_count = 4

# Simulation settings. Those set to None depend on the frequency. They
# will be set to their actual values within the frequency scan loop.
simulation_settings = {
    "outputs": {
        "directory": output_directory,
        "file_basename": None,
        "evol_time_step": 5e-12,
        "final_time": 1e-9,
        "evol_columns": [
            "t", "<Mx>", "<My>", "<Mz>", "Hx", "Hy", "Hz",
            "E_ex", "E_demag", "E_zeeman", "E_tot"
        ],
        "mag_config_every": False
    },
    "mesh": {
        "filename": "ellipsoid.msh",
        "length_unit": 1e-9,
        "volume_regions": { "ellipsoid_volume": { "alpha_LLG": 0.05 } },
        "surface_regions": { "ellipsoid_surface": {} }
    },
    "initial_magnetization": [0, 0, 1],
    "Bext": None,
    "demagnetizing_field_solver": { "nb_threads": thread_count },
    "finite_element_solver": { "max(iter)": 2000 },
    "time_integration": {
        "max(du)": 0.01,
        "min(dt)": 2e-14,
        "max(dt)": 2e-12
    }
}

# Create the output directory if needed.
if (os.path.isdir(output_directory)):
    print(f"Directory {output_directory} already exists.")
else:
    os.system(f"mkdir {output_directory}")

# Frequency scan.
output_file = open(output_path, "w")
output_file.write("#f\tamplitude\tphase\n")
for i in range(0, nb_steps_frequency):

    # Frequency-specific settings.
    freq = start_frequency \
        + i * (stop_frequency-start_frequency) / (nb_steps_frequency-1)
    omega = 2 * pi * freq
    setting_path = f"{output_directory}/settings_{freq/1e9:.1f}GHz.json"
    file_basename = f"response_{freq/1e9:.1f}GHz"  # for the .evol file
    evol_path = f"{output_directory}/{file_basename}.evol"
    simulation_settings["outputs"]["file_basename"] = file_basename
    simulation_settings["Bext"] = [
        f"{A}*cos({omega}*t)", f"{A}*sin({omega}*t)", "0"
    ]

    # Run the simulation.
    with open(setting_path, "w") as settings_file:
        json.dump(simulation_settings, settings_file, indent=4)
    print(f"\nJSON file {setting_path} generated.")
    sys.stdout.flush()
    process = subprocess.run(["../feellgood", setting_path])
    if process.returncode != 0:
        print("FeeLLGood failed -- aborting.")
        sys.exit(1)

    # Get amplitude and phase from the last line of the evolution file.
    with open(evol_path, "r") as evol_file:
        for line in evol_file:
            pass
    data = line.split()
    mx = float(data[1])
    my = float(data[2])
    Hx = float(data[4])
    Hy = float(data[5])
    amplitude = sqrt(mx**2 + my**2)
    angle_m = atan2(my, mx)
    angle_H = atan2(Hy, Hx)
    phase = remainder(angle_m - angle_H, 2 * pi)
    output_file.write(f"{freq}\t{amplitude}\t{phase}\n")
    output_file.flush()
