import matplotlib.pyplot as plt
import numpy as np

def Plot_Line(fileName):
	with open(fileName, "r") as f:
		axisLabel = f.readline().split()
		lineLabel = f.readline().split()
		xy = np.loadtxt(f, usecols = np.arange(0, len(lineLabel) ))
	
	for i in range(1,len(lineLabel)):
		plt.plot(xy[:,0], xy[:,i], "-")

	plt.xlabel(axisLabel[0])
	plt.ylabel(axisLabel[1])
	plt.legend(lineLabel[1:])

	plt.savefig("./figures/"+ fileName[0:-4] + ".png")
	plt.close()

def Plot_LineSlope(fileName):
	with open(fileName, "r") as f1:
		axisLabel = f1.readline().split()
		lineLabel = f1.readline().split()
		xy = np.loadtxt(f1, usecols = np.arange(0, len(lineLabel) ))
	
	xySlope = np.ndarray([np.size(xy,0)-1, np.size(xy,1)])
	xySlope[:,0] = xy[0:-1,0]	
	
	n = np.size(xySlope, 0)

	for i in range(1,len(lineLabel)):
		xySlope[:,i] = np.diff(xy[:,i])/np.diff(xy[:,0])
		plt.plot(xySlope[0:n,0], xySlope[0:n,i], ".-")
	
	plt.xlabel(axisLabel[0])
	plt.ylabel(r"$\Delta$" + axisLabel[1] + r"/$\Delta$" + axisLabel[0])
	plt.legend(lineLabel[1:])
	
	plt.savefig("./figures/"+ fileName[0:-4] + "_slope.png")
	plt.close()

	with open(fileName[0:-4] + "_slope.dat", "w") as f2:
		f2.write(np.array2string(xySlope))

def Plot_Line_2(fileName_1, fileName_2):
	with open(fileName_1, "r") as f1:
		axisLabel_1 = f1.readline().split()
		lineLabel_1 = f1.readline().split()
		xy_1 = np.loadtxt(f1, usecols = np.arange(0, len(lineLabel_1) ))
	
	with open(fileName_2, "r") as f2:
		lineLabel_2 = f2.readline().split()
		xy_2 = np.loadtxt(f2, usecols = np.arange(0, len(lineLabel_2) ))
	
	for i in range(1,len(lineLabel_1)):
		plt.plot(xy_1[:,0], xy_1[:,i], "-")
		plt.plot(xy_2[:,0], xy_2[:,i], ":")

	plt.xlabel(axisLabel_1[0])
	plt.ylabel(axisLabel_1[1])
	plt.legend(["appSol_f", "exSol_f", "appSol_s", "exSol_s"])

	plt.savefig("./figures/CompEx.png")
	plt.close()

"""
#OVS & PIM
Plot_Line("Temp_f.dat")
Plot_Line("Temp_s.dat")
Plot_Line("discErr_f.dat")
Plot_Line("discErr_s.dat")
Plot_LineSlope("discErr_f.dat")
Plot_LineSlope("discErr_s.dat")
Plot_Line("discLoc_f.dat")
Plot_Line("discLoc_s.dat")
"""
#VisualMotion of thermocline
Plot_Line("Motion_f.dat")
Plot_Line("Motion_s.dat")

#Compare with Exact Solutions
Plot_Line_2("charge.dat", "sol-exact.dat")

#StorageCycle
Plot_Line("charge.dat")
Plot_Line("dischg.dat")

