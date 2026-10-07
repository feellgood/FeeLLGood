#!/usr/bin/env python3

"""
Benchmark of the computation of the demagnetizing field with scalfmm 3.

For each mesh, and each set of fast multipole parameters (order, tree height, group size), this
script measures:
- the accuracy: relative L2 error on the magnetic scalar potential phi at the initial time, versus a
  reference computation (by default a higher order),
- for each number of threads: the median, min and max durations of one computation of the
  magnetostatics, as reported by feellgood in verbose mode ("magnetostatics done in ... ms").

The output is a tab separated text file, one line per (mesh, parameters, number of threads).
"""

import os
import sys
import json
import socket
import statistics
import subprocess
import tempfile
from datetime import datetime
from math import log2, floor, sqrt

__version__ = '1.0.0'

def makeSettings(mesh, volume_name, surface_name, out_dir, nbThreads, final_time, order, height,
                 group_size):
    """ returns a dictionary of settings for feellgood input """
    settings = {
        "outputs": {
            "directory": out_dir,
            "file_basename": "bench",
            "evol_time_step": 1e-12,
            "final_time": final_time,
            "evol_columns": ["t", "<Mx>", "<My>", "<Mz>", "E_demag"],
            "mag_config_every": False
        },
        "mesh": {
            "filename": mesh,
            "length_unit": 1e-9, # we use nanometers
            "volume_regions": { volume_name: {} }
        },
        "initial_magnetization": [0, 0.6, 0.8],
        "Bext": [1, 0, 1],
        "demagnetizing_field_solver": {
            "nb_threads": nbThreads,
            "order": order,
            "tree_height": height,
            "group_size": group_size
        },
        "time_integration": {
            "min(dt)": 5e-18,
            "max(dt)": 1e-12,
            "max(du)": 0.1
        }
    }
    if surface_name:
        settings["mesh"]["surface_regions"] = { surface_name: {} }
    return settings

def nbPhysicalCores():
    """ number of physical cores (hyperthreads not counted), from lscpu, else os.cpu_count() """
    try:
        out = subprocess.run(["lscpu", "-p=Core,Socket"], text=True, capture_output=True).stdout
        cores = {line for line in out.splitlines() if line and not line.startswith('#')}
        if cores:
            return len(cores)
    except OSError:
        pass
    return os.cpu_count()

def makeListNbThreads():
    """ powers of two from ncores/4 to ncores, ncores being the number of physical cores, so that
    with OMP_PLACES=cores two threads never share a core """
    maxNbThreads = 2**floor(log2(nbPhysicalCores()))
    listNbThreads = []
    nb = maxNbThreads
    while nb >= max(1, maxNbThreads//4):
        listNbThreads.append(nb)
        nb = nb // 2
    return listNbThreads

def ompEnv():
    """ environment of the runs: OMP_MAX_TASK_PRIORITY, OMP_PROC_BIND and OMP_PLACES get default
    values (one thread per physical core) unless they are already set """
    env = dict(os.environ)
    env.setdefault("OMP_MAX_TASK_PRIORITY", "11") # needed by the task priorities of scalfmm 3
    env.setdefault("OMP_PROC_BIND", "close")
    env.setdefault("OMP_PLACES", "cores")
    return env

def run(executable, settings):
    """ runs feellgood in verbose mode with seed=2 for being deterministic, returns stdout """
    env = ompEnv()
    val = subprocess.run([executable, "-v", "--seed", "2", "-"], input=json.dumps(settings),
                         text=True, capture_output=True, env=env)
    if val.returncode != 0:
        sys.exit("feellgood failed:\n" + val.stdout + val.stderr)
    return val.stdout

def parse_log(log):
    """ returns the durations of the magnetostatics (ms), the tree height and the number of nodes """
    durations = []
    height = None
    nb_nodes = None
    for line in log.splitlines():
        words = line.split()
        if line.startswith("magnetostatics done in"):
            durations.append(float(words[3]))
        elif line.startswith("Magnetostatics: order"):
            height = int(words[words.index("height") + 1].rstrip(','))
        elif line.strip().startswith("nodes:"):
            nb_nodes = int(words[1])
    return durations, height, nb_nodes

def read_phi(sol_file):
    """ returns the list of phi values (last column) of a .sol file """
    phi = []
    with open(sol_file) as f:
        for line in f:
            if not line.startswith('#') and line.strip():
                phi.append(float(line.split()[-1]))
    return phi

def compute_phi(executable, mesh, volume_name, surface_name, nbThreads, order, height, group_size):
    """ single magnetostatics computation (final_time = 0), returns phi """
    with tempfile.TemporaryDirectory() as out_dir:
        settings = makeSettings(mesh, volume_name, surface_name, out_dir, nbThreads, 0, order,
                                height, group_size)
        settings["outputs"]["mag_config_every"] = 1 # writes bench_iter0.sol
        run(executable, settings)
        return read_phi(os.path.join(out_dir, "bench_iter0.sol"))

def relative_error(phi, phi_ref):
    """ relative L2 error """
    num = sum((a - b)**2 for a, b in zip(phi, phi_ref))
    den = sum(b**2 for b in phi_ref)
    return sqrt(num/den)

def bench_mesh(f, args, mesh, volume_name, surface_name, mesh_label):
    """ loop over the fmm parameters and the numbers of threads for a single mesh """
    maxThreads = max(args.nbThreads)
    phi_ref = compute_phi(args.executable, mesh, volume_name, surface_name, maxThreads,
                          args.ref_order, args.ref_height, args.ref_group_size)
    for order in args.orders:
        for height in args.heights:
            for group_size in args.group_sizes:
                phi = compute_phi(args.executable, mesh, volume_name, surface_name, maxThreads,
                                  order, height, group_size)
                err = relative_error(phi, phi_ref)
                for nbThreads in args.nbThreads:
                    with tempfile.TemporaryDirectory() as out_dir:
                        settings = makeSettings(mesh, volume_name, surface_name, out_dir, nbThreads,
                                                args.final_time, order, height, group_size)
                        log = run(args.executable, settings)
                    durations, used_height, nb_nodes = parse_log(log)
                    # the first computation builds the interaction lists, it is not counted
                    timed = durations[1:] if len(durations) > 1 else durations
                    line = [mesh_label, nb_nodes, order, height, used_height, group_size,
                            nbThreads, len(timed), "{:.2f}".format(statistics.median(timed)),
                            "{:.2f}".format(min(timed)), "{:.2f}".format(max(timed)),
                            "{:.3e}".format(err)]
                    f.write('\t'.join(str(x) for x in line) + '\n')
                    f.flush()
                    print('\t'.join(str(x) for x in line))

def metadata(args):
    """ header of the output file """
    def cmd(c):
        try:
            return subprocess.run(c, text=True, capture_output=True).stdout.strip()
        except OSError:
            return ""
    cpu = ""
    for line in cmd(["lscpu"]).splitlines():
        if line.startswith("Model name:"):
            cpu = line.split(':', 1)[1].strip()
    version = cmd([args.executable, "--version"]).splitlines()
    lines = ["feellgood: " + (version[0] if version else "?"),
             "date: " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
             "host: " + socket.gethostname(),
             "cpu: " + cpu + " (nproc = " + str(os.cpu_count()) + ", physical cores = "
                 + str(nbPhysicalCores()) + ")",
             "environment: " + ' '.join(k + '=' + v for k, v in sorted(ompEnv().items())
                                        if k.startswith(("OMP_", "GOMP_", "KMP_"))),
             "reference: order " + str(args.ref_order) + ", tree height " + str(args.ref_height)
                 + ", group size " + str(args.ref_group_size),
             "final_time: " + str(args.final_time)]
    columns = ["mesh", "nodes", "order", "height", "used_height", "group_size", "threads",
               "n_calls", "median_ms", "min_ms", "max_ms", "err_phi"]
    return ''.join("# " + l + '\n' for l in lines) + '\t'.join(columns) + '\n'

def get_params():
    """ command line parser """
    import argparse
    description = 'feellgood benchmark of the demagnetizing field computation (scalfmm 3)'
    epilogue = '''
    examples:
    ./benchmark-fmm.py -s 4 3 2.5 -o 5 7 -H 0 4 5 6
        cylinders (radius 96 nm, height 32 nm) of average element sizes 4, 3 and 2.5 nm, orders 5
        and 7, automatic tree height and heights 4, 5, 6
    ./benchmark-fmm.py -m ../examples/ellipsoid.msh --volume ellipsoid_volume \\
            --surface ellipsoid_surface -n 1 8
        a given mesh, on 1 and 8 threads

    tree height 0 means automatic. OMP_MAX_TASK_PRIORITY, OMP_PROC_BIND and OMP_PLACES default to
    11, close and cores if they are not set: with at most as many threads as physical cores, two
    threads never share a core. All OMP_*, GOMP_* and KMP_* variables are written in the header of
    the output file.
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-e', '--executable', default='../src/feellgood',
                        help='feellgood executable (default: ../src/feellgood)')
    parser.add_argument('-m', '--mesh', help='mesh file, instead of generated cylinders')
    parser.add_argument('--volume', default='volume', help='volume region name of --mesh')
    parser.add_argument('--surface', help='surface region name of --mesh, if any')
    parser.add_argument('-s', '--sizes', type=float, nargs='+', default=[4.0, 3.0, 2.5],
                        help='average element sizes of the generated cylinders (nm)')
    parser.add_argument('-o', '--orders', type=int, nargs='+', default=[5, 7],
                        help='interpolation orders')
    parser.add_argument('-H', '--heights', type=int, nargs='+', default=[0],
                        help='tree heights, 0 means automatic')
    parser.add_argument('-g', '--group_sizes', type=int, nargs='+', default=[64],
                        help='group sizes')
    parser.add_argument('-n', '--nbThreads', type=int, nargs='+', default=makeListNbThreads(),
                        help='numbers of threads')
    parser.add_argument('-t', '--final_time', type=float, default=1e-11,
                        help='final physical simulation time of the timing runs (s)')
    parser.add_argument('--ref_order', type=int, default=10, help='order of the reference')
    parser.add_argument('--ref_height', type=int, default=0,
                        help='tree height of the reference, 0 means automatic')
    parser.add_argument('--ref_group_size', type=int, default=64,
                        help='group size of the reference')
    parser.add_argument('-f', '--output', help='output file (default: benchmark-fmm-<host>.txt)')
    parser.add_argument('--version', action='version', version=__version__,
                        help='show the version number')
    return parser.parse_args()

if __name__ == '__main__':
    args = get_params()
    args.executable = os.path.abspath(args.executable)
    outputFileName = args.output or 'benchmark-fmm-' + socket.gethostname() + '.txt'
    try:
        with open(outputFileName, 'w') as f:
            f.write(metadata(args))
            if args.mesh:
                bench_mesh(f, args, os.path.abspath(args.mesh), args.volume, args.surface,
                           os.path.basename(args.mesh))
            else:
                from feellgood.meshMaker import Cylinder
                height = 32
                radius = 3.0 * height
                with tempfile.TemporaryDirectory() as mesh_dir:
                    meshFileName = os.path.join(mesh_dir, "cylinder.msh")
                    for elt_size in args.sizes:
                        Cylinder(radius, height, elt_size, "surface", "volume").make(meshFileName)
                        bench_mesh(f, args, meshFileName, "volume", "surface",
                                   "cylinder-" + str(elt_size))
        print("results written in", outputFileName)
    except KeyboardInterrupt:
        print(" benchmark interrupted")
        sys.exit()
