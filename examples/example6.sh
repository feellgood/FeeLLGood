#!/bin/bash

# This example shows how multiple feeLLGood simulations can be run from
# a shell script. The goal is to simulate ferromagnetic resonance on an
# elongated nanoparticle. The sample is a prolate ellipsoid of 8×40 nm,
# with a 1 T/µ₀ magnetization initially along the long (Z) axis. There
# is no DC field. A 4 mT rotating field is applied on the XY plane. The
# frequency of the field is varied in order to find the resonance.

# This must be executed first, before $_ gets clobbered by subsequent
# commands.
full_script_path="$_"

# Main simulation parameters.
amplitude=0.004    # magnitude of the rotating field [T]
duration=1e-9      # duration of the excitation [s]
f_start=9.3        # start frequency [GHz]
f_step=0.2         # frequency step [GHz]
f_stop=13.3        # stop frequency [GHz]
f_format='%04.1f'  # format for the frequency in file names

# Create a directory named "resonance.d" under the directory holding
# this script, and use it as a working directory.
cd $(dirname $full_script_path)
mkdir -p resonance.d
cd resonance.d
echo "Storing all results in $PWD"

# Template configuration file. {FREQUENCY} is a placeholder that will be
# replaced with the frequency, in GHz, of each individual simulation.
# Note that, at the end of the simulation, the applied field is along X.
cat > settings-template.yml <<EOF
outputs:
  file_basename: response_{FREQUENCY}GHz
  evol_time_step: 5e-12
  final_time: $duration
  evol_columns: [t, <Mx>, <My>, <Mz>]
  evol_header: true
  take_photo: false
mesh:
  filename: ../ellipsoid.msh
  scaling_factor: 1e-9
  volume_regions:
    ellipsoid_volume:
      Ae: 1e-11
      Js: 1
      alpha_LLG: 0.05
  surface_regions:
    ellipsoid_surface:  # use defaults
initial_magnetization: [0, 0, 1]
Bext: [
    $amplitude * cos(2 * pi * {FREQUENCY}e9 * (t - $duration)),
    $amplitude * sin(2 * pi * {FREQUENCY}e9 * (t - $duration)),
    0
]
time_integration:
  max(du): 0.01
  min(dt): 2e-14
  max(dt): 2e-12
EOF

# Prepare the output file. The columns are:
# f: frequency in GHz
# m_x: magnetization along the excitation field
# m_y: magnetization perpendicular to the excitation field
echo -e "#f(GHz)\tm_x\tm_y" > resonance.tsv

# Frequency scan.
echo -n "Scanning from $f_start GHz to $f_stop GHz "
echo "in steps of $f_step GHz:"
for f in $(seq -f "$f_format" $f_start $f_step $f_stop); do
    echo "    f = $f GHz"
    output=output_${f}GHz.log   # stdout and stderr of feellgood
    evol=response_${f}GHz.evol  # evolution file

    # Run the simulation.
    sed "s/{FREQUENCY}/$f/" settings-template.yml | \
    ../../feellgood - > $output 2>&1

    # Extract response from last line of the evolution.
    last_line=$(tail -n 1 $evol)
    magnetization=$(cut -f2,3 <<< $last_line)
    echo -e "$f\t$magnetization" >> resonance.tsv
done
echo "Done."

# Make sure gnuplot is available.
if ! command -v gnuplot > /dev/null; then
    echo "gnuplot command not found: not plotting the results."
    exit
fi

# Plot the frequency response.
gnuplot -persist <<EOF
set title "Ferromagnetic resonance (4 mT rotating field)"
set xlabel "f (GHz)"
set ylabel "response (reduced magnetization)"
set xrange [$f_start:$f_stop]
set grid
set style data linespoints
plot 'resonance.tsv' using 1:(sqrt(\$2**2+\$3**2)) title "amplitude", \
    '' using 1:2 title "in phase", \
    '' using 1:3 title "in quadrature"
EOF
