import json

## \class Settings
# this class gives an efficient way to build and write a json file from a dictionary to be able to script easily feellgood executable

class Settings(object):

    ## \brief creator
    # this is creator to build a nested dictionary to mimmic json settings
    def __init__(self):
        ## \brief dictionary
        # mySets is the strict equivalent to the json file needed by executable feellgood to run a simulation, any (key,value) not corresponding
        # to what is expected by feellgood is just ignored.
        self.mySets = {}

        outputs = {}
        outputs["directory"] = "data_out/"
        outputs["file basename"] = "tom"
        outputs["vtk file"] = False
        outputs["evol time step"] = 1e-7
        outputs["take_photo"] = 5
        outputs["verbose"] = True
        outputs["evol columns"] = ["iter","t","dt","max dm","<mx>","<my>","<mz>","E_ex","E_aniso","E_demag","E_zeeman","E_tot","DW_z","DW_dz"]
        outputs["evol header"] = True

        self.mySets["outputs"] = outputs

        mesh = {}
        mesh["filename"] = "wire_d70_L1000.msh"
        mesh["scaling factor"] = 1e-9
        mesh["epsilon"] = 1e-40
        vol300 = {"Ae":1e-11, "Js":1.0, "Ka":0.0, "a" : [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]], "alpha":0.05}

        mesh["volume_regions"] = {"300" : vol300}

        surf200 = {"Js" : 1.0, "Ks" : 0.0, "uk" : [0.0, 0.0 ,1.0]}
        mesh["surface_regions"] = {"200" : surf200}

        self.mySets["mesh"] = mesh

        self.mySets["restore"] = False
        self.mySets["restore from file"] = "sol.in"
        self.mySets["initial magnetization"] = {"Mx" : "1-z^2", "My" : "0" , "Mz": "tanh(z*20)"}

        recentering = {}
        recentering["recenter"] = False
        recentering["direction"] = 'Z'
        recentering["threshold"] = 0.1
        self.mySets["recentering"] = recentering

        self.mySets["Bext"] = {"Bx" : "cos(2*Pi*0.000000001*t)", "By" : "sin(2*Pi*0.000000001*t)" , "Bz": "0"}

        current = {}
        current["Uz"] = 0.0
        current["beta"] = 0.0
        self.mySets["spin polarized current"] = current

        FEMsolver = {}
        FEMsolver["nbThreads"] = 16
        FEMsolver["max(iter)"] = 500
        FEMsolver["refresh preconditionner every"] = 20

        self.mySets["finite element solver"] = FEMsolver

        DEMAGsolver = {}
        DEMAGsolver["nbThreads"] = 16
        
        self.mySets["demagnetization field solver"] = DEMAGsolver

        timing = {}
        timing["final_time"] = 1e-6
        timing["min(du)"] = 1e-9
        timing["max(du)"] = 0.02
        timing["min(dt)"] = 1e-14
        timing["max(dt)"] = 1e-7
        timing["initial dt"] = 1e-9
        
        self.mySets["time integration"] = timing
    
    ## \brief getter
    # standard getter
    #
    def __getitem__(self,key):
        return self.mySets[key]

    ## \brief setter
    # (key,value) setter
    def __setitem__(self,key,value):
        self.mySets[key] = value
    
    ## \brief write json file
    # this method write a json file from the dictionary built by the creator
    def write(self,fileName):
        with open(fileName,'w') as outfile:
            json.dump(self.mySets,outfile,indent = 4)
