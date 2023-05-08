import shutil
from decimal import *
getcontext().prec = 8 # memorize 8 digits after the comma

import numpy as np
import myTools


if __name__ == '__main__':

#	param2D = myTools.Create2DListFromTXT(NameTXTinput='inf_param.dat')
#	myTools.CreateNewDir(dirName='out', isLocal=True)
#	shutil.copy("../inf_param.dat", "./")
#	
#	if len(param2D) > 0 :
#		print "Number of lines in parameter file: ", len(param2D)
#		ln = 0
#		for line in param2D :
#			
#			m = line[0] # potential parameter
#			l = line[1] # potential parameter
#			lp = line[2] # potential parameter
#			xi = line[3] # potential parameter
#			lpp = line[4] # baryogenesis parameter
#			x0 = line[5] # initial value of x
#			pix0 = 0.0 # initial value of pix
#			y0 = line[6] # initial value of y
#			piy0 = 0.0
#			psi0 = 0.0
#			Nstart = 0.0 # initial number of e-folding (for odeint)
#			Nsteps = int(1000+1) # number of steps in NN list (its length; +1 is for nice purposes)
#			IntMethod = 'RK45'
#			#IntMethod = 'LSODA'
#
#			mllpxilpp = (m, l, lp, xi, lpp)
#
#			data = myTools.inflate(xpixypiypsi0=[x0, pix0, y0, piy0, psi0], mllpxilpp=mllpxilpp, Nstart=Nstart, Nsteps=Nsteps, eps=1e-10, method='RK45', rtol=1.e-4, atol=1.e-7)
#			if data == [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1] :
#				continue
#			else :
#				Nend, Nstar, As_star, ns_star, r_star, etaBlate, T_RR_stop, T_RS_stop, T_SR_stop, T_SS_stop, betaiso_stop, cosD_stop, cosD0_stop = data
#			print "\nNend = ", Nend
#			print "Nstar = ", Nstar
#			print "As_star = ", As_star
#			print "ns_star = ", ns_star
#			print "r_star = ", r_star
#			print "etaBlate = ", etaBlate
#			print "T_RR_stop = ", T_RR_stop
#			print "T_RS_stop = ", T_RS_stop
#			print "T_SR_stop = ", T_SR_stop
#			print "T_SS_stop = ", T_SS_stop
#
#			print "Done line %d!" % ln
#			ln += 1
	

	m = 6.43E-07
	l = 8.73E-12
	lp = 6.89E-13
	xi = 0.0596
	lpp = 8.67E-05
	x0 = 18.4
	y0 = 6.63
	param = (m, l, lp, xi, lpp, x0, y0)
	fname = 'chainMod1'

	matry = 1e1
	print "len(initial chain) = %d\n" % int(matry)

	print "Start to build the chain..."
	myTools.MCMC(fname=fname, param=param, matry=matry, parchange=0.05, Nstart=0.0, Nsteps=1e3, eps=1e-10, method='RK45', rtol=1.e-4, atol=1.e-7, isrescale=False, isconvfact=True)
	print "Done!"

