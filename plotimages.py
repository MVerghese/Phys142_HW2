import pylab

x0 = -4
dx = 8/600

xPot = []
yPot = []

for i in range (0,600):
    x = x0 + dx*i
    xPot.append(x)
    y = (x**2 - 2**.5)**2
    yPot.append(y)

for i in range(0,351):
    data=pylab.loadtxt("out" + str(i) + ".txt")

    print(data)

    x = []
    y = []

    for datapoint in data:
        x.append(datapoint[0])
        y.append(datapoint[1])



    pylab.clf()
    fig, axes = pylab.subplots()
    axes.set_xlim([-4,4])
    axes.set_ylim([-0,1.25])
    #axes.set_aspect('equal')

    #circle = pylab.Circle((data[i,2], 0), radius=0.5, fc='b')
    #xes.add_patch(circle)

    pylab.plot(x,y);
    pylab.plot(xPot,yPot)

    pylab.legend()
    pylab.title("Quantam Tunneling at time step, t= %.2f " %  i)
    pylab.xlabel("Position (m)")
    #pylab.grid(True)
        
    #pylab.savefig("output.png")

    pylab.savefig("out"+("%03d" % i)+".png")
