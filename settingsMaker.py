import json

## \class Settings
# this class gives an efficient way to build and write a json file from a dictionary to be able to script easily feellgood executable

class Settings(object):

    ## \brief creator
    # this is creator to build a nested dictionary to mimmic json settings
    def __init__(self,mshFileName):
        ## \brief dictionary
        # mySets is the strict equivalent to the json file needed by executable feellgood to run a simulation, any (key,value) not corresponding
        # to what is expected by feellgood is just ignored.
        self.mySets = {}

        self.mySets["outputs"] = {"directory" : "data_out/",
        "file_basename" : "tom",
        "vtk_file" : False,
        "evol_time_step" : 1e-7,
        "take_photo" : 5,
        "evol_columns" : ["iter","t","dt","max_dm","<Mx>","<My>","<Mz>","E_ex","E_aniso","E_demag","E_zeeman","E_tot","DW_z","DW_dz"],
        "evol_header" : True }

        self.mySets["mesh"] = {"filename" : mshFileName, "scaling_factor" : 1e-9, "epsilon" : 1e-40, "volume_regions" : {}, "surface_regions" : {} }

        self.mySets["initial_magnetization"] = {"Mx" : "1-z^2", "My" : "0" , "Mz": "tanh(z*20)"}

        self.mySets["Bext"] = {"Bx" : "0", "By" : "0" , "Bz": "0"}

        self.mySets["spin_polarized_current"] = {"Uz" : 0.0, "beta" : 0.0}

        self.mySets["finite_element_solver"] = {"nb_threads" : 16, "max(iter)" : 500, "refresh_preconditionner_every": 20}

        self.mySets["demagnetization_field_solver"] = {"nb_threads" : 16}

        self.mySets["time_integration"] = {"final_time" : 1e-6, "min(du)" : 1e-9, "max(du)" : 0.02, "min(dt)" : 1e-11, "max(dt)" : 1e-7, "initial_dt" : 1e-9}
    
    ## \brief getter
    # standard getter
    #
    def __getitem__(self,key):
        return self.mySets[key]

    ## \brief setter
    # (key,value) setter
    def __setitem__(self,key,value):
        self.mySets[key] = value

    ## \brief define recentering
    # direction and threshold input paramaters defining recentering
    def recentering(self,direction,threshold):
        self.mySets["recentering"] = {"direction" : direction, "threshold" : threshold}

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
