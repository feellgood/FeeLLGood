import json

mySettings = {}

outputs = {}
outputs["directory"] = "data_out/"
outputs["file basename"] = "tom"
outputs["vtk file"] = False
outputs["save_energies"] = 3
outputs["take_photo"] = 5

mySettings["outputs"] = outputs

mesh = {}
mesh["filename"] = "wire_d70_L1000.msh"
mesh["scaling factor"] = 1e-9
mesh["epsilon"] = 1e-40
vol300 = {"Ae":1e-11, "Js":1.0, "Ka":0.0, "a" : [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]], 
        "alpha":0.05
        }

mesh["volume_regions"] = {"300" : vol300}

surf200 = {"Js" :"-1", "Ks" : 0.0, "uk" : [0.0, 0.0 ,1.0]}
mesh["surface_regions"] = {"200" : surf200}

mySettings["mesh"] = mesh

mySettings["restore"] = False
mySettings["restore from file"] = "sol.in"
mySettings["initial magnetization"] = {"Mx" : "1-z^2", "My" : "0" , "Mz": "tanh(z*20)"}

recentering = {}
recentering["recenter"] = False
recentering["direction"] = 'Z'
recentering["threshold"] = 0.1
mySettings["recentering"] = recentering

mySettings["Bext"] = [0.0,0.0,0.007]

current = {}
current["Uz"] = 0.0
current["beta"] = 0.0
mySettings["spin polarized current"] = current


FEMsolver = {}
FEMsolver["nbThreads"] = 16
FEMsolver["max(iter)"] = 500
FEMsolver["refresh preconditionner every"] = 20

mySettings["finite element solver"] = FEMsolver

DEMAGsolver = {}
DEMAGsolver["nbThreads"] = 16
DEMAGsolver["analytic corrections"] = True

mySettings["demagnetization field solver"] = DEMAGsolver

timing = {}
timing["final_time"] = 1e-7
timing["theta"] = 0.5
timing["min(du)"] = 1e-9
timing["max(du)"] = 0.02
timing["min(dt)"] = 1e-14
timing["max(dt)"] = 1e-7
timing["initial dt"] = 1e-9
mySettings["time integration"] = timing

with open('mySettings.json','w') as outfile:
    json.dump(mySettings,outfile)
