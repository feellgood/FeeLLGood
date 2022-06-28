import os
import sys
import subprocess
from math import sqrt,pi

sys.path.insert(0,'../tools')
from settingsMaker import Settings

mySettings = Settings("../ellipsoid.msh")
mySettings.createVolRegion( "300" )
mySettings.createSurfRegion( "200" )

MaxNbThreads = int(subprocess.check_output(["getconf","_NPROCESSORS_ONLN"]))

mySettings["finite_element_solver"]["nb_threads"] = 4 # MaxNbThreads
mySettings["finite_element_solver"]["max(iter)"] = 700

mySettings["demagnetizing_field_solver"]["nb_threads"] = 4 # MaxNbThreads

mySettings["outputs"]["evol_columns"] = ["t","<Mx>","<My>","<Mz>","<dMx/dt>","<dMy/dt>","<dMz/dt>","Hx","Hy","Hz","E_ex","E_demag","E_zeeman","E_tot"]
mySettings["outputs"]["take_photo"] = 500
mySettings["outputs"]["directory"] = "test_data_out/"

mySettings["outputs"]["evol_time_step"] = 0.1e-9

mySettings["mesh"]["scaling_factor"] = 1e-10

mySettings["mesh"]["volume_regions"]["300"]["alpha_LLG"] = 0.05

mySettings["time_integration"]["final_time"] = 0.2e-7
mySettings["time_integration"]["min(dt)"] = 0.1e-11
mySettings["time_integration"]["max(dt)"] = 0.5e-10

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
	mySettings["Bext"] = {"Bx" : str(A) + "*cos(" + str(omega) + "*t)", "By" : str(A) + "*sin("  + str(omega) +  "*t)" , "Bz": "0"}
	
	mySettings.write('mySettings.json')
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

