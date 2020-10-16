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
        outputs["file_basename"] = "tom"
        outputs["vtk_file"] = False
        outputs["evol_time_step"] = 1e-7
        outputs["take_photo"] = 5
        outputs["verbose"] = True
        outputs["evol_columns"] = ["iter","t","dt","max_dm","<Mx>","<My>","<Mz>","E_ex","E_aniso","E_demag","E_zeeman","E_tot","DW_z","DW_dz"]
        outputs["evol_header"] = True

        self.mySets["outputs"] = outputs

        mesh = {}
        mesh["filename"] = "wire_d70_L1000.msh"
        mesh["scaling_factor"] = 1e-9
        mesh["epsilon"] = 1e-40
        mesh["volume_regions"] = {}

        mesh["surface_regions"] = {}

        self.mySets["mesh"] = mesh

        self.mySets["restore"] = False
        self.mySets["restore_from_file"] = "sol.in"
        self.mySets["initial_magnetization"] = {"Mx" : "1-z^2", "My" : "0" , "Mz": "tanh(z*20)"}

        recentering = {}
        recentering["recenter"] = False
        recentering["direction"] = 'Z'
        recentering["threshold"] = 0.1
        self.mySets["recentering"] = recentering

        self.mySets["Bext"] = {"Bx" : "cos(2*Pi*42*t)", "By" : "sin(2*Pi*42*t)" , "Bz": "0"}

        self.mySets["spin_polarized_current"] = {"Uz" : 0.0, "beta" : 0.0}

        FEMsolver = {}
        FEMsolver["nb_threads"] = 16
        FEMsolver["max(iter)"] = 500
        FEMsolver["refresh_preconditionner_every"] = 20

        self.mySets["finite_element_solver"] = FEMsolver

        self.mySets["demagnetization_field_solver"] = { "nb_threads" : 16}

        timing = {}
        timing["final_time"] = 1e-6
        timing["min(du)"] = 1e-9
        timing["max(du)"] = 0.02
        timing["min(dt)"] = 1e-14
        timing["max(dt)"] = 1e-7
        timing["initial_dt"] = 1e-9
        
        self.mySets["time_integration"] = timing
    
    ## \brief getter
    # standard getter
    #
    def __getitem__(self,key):
        return self.mySets[key]

    ## \brief setter
    # (key,value) setter
    def __setitem__(self,key,value):
        self.mySets[key] = value

    ## \brief create a new volume region
    # (key) , initialized with def values
    def createVolRegion(self,key):
        if key not in self.mySets["mesh"]["volume_regions"]:
            self.mySets["mesh"]["volume_regions"].update( { key : {"Ae":1e-11, "Js":1.0, "K":0.0, "uk":[0.0,0.0,0.0],"K3":0.0, "uk3":[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], "alpha":0.05} } )

    ## \brief create a new surface region
    # (key) , initialized with def values
    def createSurfRegion(self,key):
        if key not in self.mySets["mesh"]["surface_regions"]:
            self.mySets["mesh"]["surface_regions"].update( { key : {"Js" : 1.0, "Ks" : 0.0, "uk" : [0.0, 0.0 ,1.0]} } )

    ## \brief write json file
    # this method write a json file from the dictionary built by the creator
    def write(self,fileName):
        with open(fileName,'w') as outfile:
            json.dump(self.mySets,outfile,indent = 4)
