import os.path
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

        self.mySets["mesh"] = {"filename" : mshFileName, "scaling_factor" : 1e-9, "volume_regions" : {}, "surface_regions" : {} }

        self.mySets["initial_magnetization"] = {"Mx" : "1-z^2", "My" : "0" , "Mz": "tanh(z*20)"}

        self.mySets["Bext"] = {"Bx" : "0", "By" : "0" , "Bz": "0"}

        self.mySets["finite_element_solver"] = {"nb_threads" : 16, "max(iter)" : 500, "refresh_preconditioner_every": 20}

        self.mySets["demagnetizing_field_solver"] = {"nb_threads" : 16}

        self.mySets["time_integration"] = {"final_time" : 1e-6, "max(du)" : 0.02, "min(dt)" : 1e-11, "max(dt)" : 1e-7}
    
    ## \brief getter
    # standard getter
    #
    def __getitem__(self,key):
        return self.mySets[key]

    ## \brief setter
    # (key,value) setter
    def __setitem__(self,key,value):
        self.mySets[key] = value

    ## \brief define spin_polarized_current
    # uk and beta input parameters defining spin transfert torque
    def createSTT(self,volRegRef,gamma0,sigma,N0,beta,l_J,l_sf,bc_surf_reg_ref1,bc_val_name1,bc_val1,bc_surf_reg_ref2,bc_val_name2,bc_val2):
        self.mySets["spin_transfer_torque"] = {"volume_region_reference" : str(volRegRef), "gamma0" : gamma0, "sigma" : sigma,"dens_state" : N0, "beta" : beta,"l_J" : l_J,"l_sf" : l_sf,"boundary_conditions" : {  str(bc_surf_reg_ref1) : { bc_val_name1 : bc_val1 } ,  str(bc_surf_reg_ref2) : { bc_val_name2: bc_val2 } } }

    ## \brief define recentering
    # direction and threshold input parameters defining recentering
    def recentering(self,direction,threshold):
        self.mySets["recentering"] = {"direction" : direction, "threshold" : threshold}

    ## \brief create a new volume region
    # (key) , initialized with def values
    def createVolRegion(self,key):
        if key not in self.mySets["mesh"]["volume_regions"]:
            self.mySets["mesh"]["volume_regions"].update( { key : {"Ae":1e-11, "Js":1, "K":0.0, "uk":[0.0,0.0,1.0],"K3":0.0, "ex":[1.0,0.0,0.0], "ey":[0.0,1.0,0.0], "ez":[0.0,0.0,1.0], "alpha_LLG":0.05} } )

    ## \brief set material constants exchange , magnetization at saturation and alpha_LLG for region 'key', create if if key does not exists
    # (key) , initialized with def values
    def setMaterialRegion(self,key,cst_exch,mag_sat,cst_alpha):
        if key not in self.mySets["mesh"]["volume_regions"]:
            self.mySets["mesh"]["volume_regions"].update( { key : {"Ae":cst_exch, "Js":mag_sat, "K":0.0, "uk":[0.0,0.0,1.0],"K3":0.0, "ex":[1.0,0.0,0.0], "ey":[0.0,1.0,0.0], "ez":[0.0,0.0,1.0], "alpha_LLG":cst_alpha} } )
        else:
            self.mySets["mesh"]["volume_regions"][key]["Ae"] = cst_exch
            self.mySets["mesh"]["volume_regions"][key]["Js"] = mag_sat
            self.mySets["mesh"]["volume_regions"][key]["alpha_LLG"] = cst_alpha
        
    ## \brief create a new surface region
    # (key) , initialized with def values
    def createSurfRegion(self,key):
        if key not in self.mySets["mesh"]["surface_regions"]:
            self.mySets["mesh"]["surface_regions"].update( { key : {"suppress_charges" : False, "Ks" : 0.0, "uk" : [0.0, 0.0 ,1.0]} } )

    ## \brief some test to check the validity of the parameters
    # this method is called before writing the json file, it might also be called from another script for some checking before calling the write method
    def check(self):
        if "spin_transfer_torque" in self.mySets:
            print("spin transfer torque subsection detected.")
            if "boundary_conditions" not in self.mySets["spin_transfer_torque"]:
                print("undefined boundary conditions.")
            elif len(self.mySets["spin_transfer_torque"]["boundary_conditions"]) == 2 :
                bc = self.mySets["spin_transfer_torque"]["boundary_conditions"]
                nbV=0
                nbJ=0
                for i in bc.values():
                    if isinstance(i,dict) and "V" in i.keys():
                        nbV += 1
                    if isinstance(i,dict) and "J" in i.keys():
                        nbJ += 1
                if nbV == 2:
                    print("Dirichlet boundary conditions.")
                elif nbV ==1 and nbJ==1:
                    print("mixt boundary conditions.")
                elif nbJ==2:
                    print("Neumann conditions.")
                else:
                    print("not supported boundary conditions.")
            else :
                print("wrong number of boundary conditions, should be 2.")
    
    ## \brief write json file
    # this method write a json file from the dictionary built by the creator
    def write(self,fileName):
        self.check()
        with open(fileName,'w') as outfile:
            json.dump(self.mySets,outfile,indent = 4)
        print("json file " + fileName + " generated.")
