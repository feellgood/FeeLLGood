import os
import sys
import subprocess
from math import sqrt,pi

from feellgood.settingsMaker import Settings

mySettings = Settings("ellipsoid.msh")
mySettings.createVolRegion( "ellipsoid_volume" )
mySettings.createSurfRegion( "ellipsoid_surface" )

MaxNbThreads = int(subprocess.check_output(["getconf","_NPROCESSORS_ONLN"]))

mySettings["finite_element_solver"]["nb_threads"] = 4 # MaxNbThreads
mySettings["finite_element_solver"]["max(iter)"] = 700

mySettings["demagnetizing_field_solver"]["nb_threads"] = 4 # MaxNbThreads

mySettings["outputs"]["evol_columns"] = ["t","<Mx>","<My>","<Mz>","<dMx/dt>","<dMy/dt>","<dMz/dt>","Hx","Hy","Hz","E_ex","E_demag","E_zeeman","E_tot"]
mySettings["outputs"]["mag_config_every"] = 500
mySettings["outputs"]["directory"] = "test_data_out/"
mySettings["outputs"]["file_basename"] = "ellipsoid"

mySettings["outputs"]["evol_time_step"] = 5e-16
mySettings["outputs"]["final_time"] = 1e-13

mySettings["mesh"]["length_unit"] = 1e-10

mySettings["mesh"]["volume_regions"]["ellipsoid_volume"]["alpha_LLG"] = 0.05

mySettings["time_integration"]["min(dt)"] = 5e-18
mySettings["time_integration"]["max(dt)"] = 2.5e-16

mySettings["time_integration"]["max(du)"] = 0.1

mySettings["initial_magnetization"] = [0.01, 0, 0.99]

if(os.path.exists(mySettings["outputs"]["directory"]) and os.path.isdir(mySettings["outputs"]["directory"]) ):
    print("directory " + mySettings["outputs"]["directory"] + " already exists.")
else:
    os.system("mkdir " + mySettings["outputs"]["directory"])


A = 1.0

startFrequency = 4.5e9
stopFrequency = 5.5e9
nbStepsFrequency = 200
currentPath = os.getcwd()
os.chdir("test_data_out")
myFile = open("FerroRes.txt","w")
myFile.write("#f	MinMax\n")
os.chdir(currentPath)
for i in range(0,nbStepsFrequency) :
    freq = startFrequency + i*(stopFrequency - startFrequency)/(nbStepsFrequency-1)
    omega = 2*pi*freq
    mySettings["Bext"] = [f"{A}*cos({omega}*t)", f"{A}*sin({omega}*t)", "0"]

    mySettings.write('mySettings.json')
    sys.stdout.flush()
    val = subprocess.run(["../feellgood","-v","mySettings.json"])

    amplitudeX = 0
    amplitudeY = 0
    if(val.returncode==0):
        print("FeeLLGood terminated correctly")
        minMx = 1
        maxMx = -1
        minMy = 1
        maxMy = -1
        with open(mySettings["outputs"]["directory"] + mySettings["outputs"]["file_basename"] + ".evol","r") as f:
            for line in f:
                lastLine = line
                data = lastLine.split()
                if (data[0] != "#"):
                    t = float(data[0])
                    mx = float(data[1])
                    my = float(data[2])
                    #print("values= " + str(t) +";"+ str(mx) +";"+ str(my))
                    if(t>1e-8):
                        if(minMx>mx):
                            minMx = mx
                        if(maxMx<mx):
                            maxMx = mx

                        if(minMy>my):
                            minMy = my
                        if(maxMy<my):
                            maxMy = my
        f.close()
        amplitudeX = maxMx - minMx
        amplitudeY = maxMy - minMy
        print("f=" +str(freq) + " amplitudeX= " + str(amplitudeX) + " amplitudeY= " + str(amplitudeY))
    myFile.write(str(freq) + "\t" + str(sqrt(amplitudeX**2 + amplitudeY**2)) + "\n")
    myFile.flush()
myFile.close()
sys.exit(0)

