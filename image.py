import pylab


data=pylab.loadtxt("out.txt")
pylab.plot(data[:,0],data[:,1])
pylab.legend()
pylab.title("Avg x over time")
pylab.xlabel("Time (s)")
pylab.ylabel("x")
pylab.savefig("output.png")
