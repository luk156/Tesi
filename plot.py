#!/usr/bin/env python

import numpy as N
import pylab

data = N.loadtxt("roots-wavevectors.csv", delimiter=",")

pylab.figure(figsize=(3.46,2.14), frameon="False")
#            xLL   yLL   W         H
pylab.axes([0.17, 0.20, 0.9-0.17, 0.75])
pylab.plot(data[:,0], data[:,1], '-k')
pylab.xlim([-15.0, 15.0])
pylab.ylim([0.0, 1.0])
pylab.xlabel(r"$k R$")
pylab.ylabel(r"$|f(k R)|$")
pylab.savefig("plot.png", dpi=300)
pylab.savefig("plot.pdf")
