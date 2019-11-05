import matplotlib.pyplot as plt
import numpy as np

def Plot_Line(fileName):
	with open(fileName, "r") as f:
		axisLabel = f.readline().split()
		xy = np.loadtxt(f, usecols=(0,1))
		
	plt.plot(xy[:,0], xy[:,1], ".-")
	plt.xlabel(axisLabel[0])
	plt.ylabel(axisLabel[1])
	plt.title(fileName)
	plt.savefig("./figures/"+ fileName[0:-4] + ".png")
	plt.close()

def Plot_LineSlope(fileName):
	with open(fileName, "r") as f:
		axisLabel = f.readline().split()
		xy = np.loadtxt(f, usecols=(0,1))
	
	xySlope = np.ndarray([np.size(xy,0)-1, np.size(xy,1)])
	xySlope[:,0] = xy[0:-1,0]
	xySlope[:,1] = np.diff(xy[:,1])/np.diff(xy[:,0])
	
	n = np.size(xySlope, 0)
	
	plt.plot(xySlope[1:n,0], xySlope[1:n,1], ".-")
	plt.xlabel(axisLabel[0])
	plt.ylabel(r"$\Delta$" + axisLabel[1] + r"/$\Delta$" + axisLabel[0])
	plt.title(fileName[0:-4] + "_slope")
	plt.savefig("./figures/"+ fileName[0:-4] + "_slope.png")
	plt.close()

	with open(fileName[0:-4] + "_slope.dat", "w") as f2:
		f2.write(np.array2string(xySlope))

Plot_Line("storstate.dat")
Plot_Line("discL1.dat")
Plot_LineSlope("discL1.dat")
Plot_Line("discL2.dat")
Plot_LineSlope("discL2.dat")
Plot_Line("discInf.dat")
Plot_LineSlope("discInf.dat")
Plot_Line("discInfLoc.dat")

