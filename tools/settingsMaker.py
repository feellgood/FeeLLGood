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
        self.mySets = {
            "outputs": {},
            "mesh": {
                "filename": mshFileName,
                "volume_regions": {},
                "surface_regions": {}
            },
            "finite_element_solver": {},
            "demagnetizing_field_solver": {},
            "time_integration": {}
        }

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
    def createSTT(self,volRegRef,gamma0,sigma,N0,beta,l_J,l_sf,bc_surf_reg_ref1,bc_val1,bc_surf_reg_ref2,bc,bc_val2):
        self.mySets["spin_transfer_torque"] = {"enable": True, "volume_region_reference" : str(volRegRef), "gamma0" : gamma0, "sigma" : sigma,"dens_state" : N0, "beta" : beta,"l_J" : l_J,"l_sf" : l_sf,"boundary_conditions" : {  str(bc_surf_reg_ref1) : { "V" : bc_val1 } ,  str(bc_surf_reg_ref2) : { "V" : bc_val2 } } }

    ## \brief define recentering
    # direction and threshold input parameters defining recentering
    def recentering(self,direction,threshold):
        self.mySets["recentering"] = {"enable": True, "direction": direction, "threshold": threshold}

    ## \brief create a new volume region
    # (key) , initialized with def values
    def createVolRegion(self,key):
        regions = self.mySets["mesh"]["volume_regions"]
        if key not in regions:
            regions[key] = {}
        return regions[key]

    ## \brief set material constants exchange , magnetization at saturation and alpha_LLG for region 'key', create if if key does not exists
    # (key) , initialized with def values
    def setMaterialRegion(self,key,cst_exch,mag_sat,cst_alpha):
        region = self.createVolRegion(key)
        region["Ae"] = cst_exch
        region["Js"] = mag_sat
        region["alpha_LLG"] = cst_alpha

    ## \brief create a new surface region
    # (key) , initialized with def values
    def createSurfRegion(self,key):
        regions = self.mySets["mesh"]["surface_regions"]
        if key not in regions:
            regions[key] = {}
        return regions[key]

    ## \brief some test to check the validity of the STT parameters
    # this method is called before writing the json file, it might also be called from another script for some checking before calling the write method
    def check(self):
        if "spin_transfer_torque" in self.mySets:
            print("spin transfer torque subsection detected.")
            if "boundary_conditions" not in self.mySets["spin_transfer_torque"]:
                print("undefined boundary conditions.")
            elif len(self.mySets["spin_transfer_torque"]["boundary_conditions"]) == 2 :
                bc = self.mySets["spin_transfer_torque"]["boundary_conditions"]
                nbV=0
                for i in bc.values():
                    if isinstance(i,dict) and "V" in i.keys():
                        nbV += 1
                if nbV == 2:
                    print("Dirichlet boundary conditions.")
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
