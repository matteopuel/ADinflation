import os, sys, shutil, datetime
import math, bigfloat
import numpy as np
from scipy import integrate
from scipy import interpolate

from decimal import *
getcontext().prec = 8 # memorize 8 digits after the comma

import random


def Create2DListFromTXT(NameTXTinput):

	"""
	This function copies the content of a txt file into a 2D python list, except for the first line.
	Input:
	* NameTXTinput = name of the input txt file
	Output: 2D list
	"""

	with open(NameTXTinput,'r') as txtinput:
	
		vector2D = [] # initialize a 2D list
	
		txtinput.next() # skip the first line in the txt file
	
		for line in txtinput:
			row = [] # initialize a 1D list corresponding to each line
	
			i = 0
			while i < len(line):
				if line[0] == '#': # skip the last line and all the lines that are comments
					break
				elif line[i] != ' ' and line[i] != '\t' and line[i] != '\n':
					string = line[i] # memorize in 'string' the fist number that it encounters
					for j in xrange(1, len(line)):
						if i+j < len(line) and line[i+j] != ' ' and line[i+j] != '\t':
							string += line[i+j] # combine the numbers that are consecutive, e.g. '9','2','1' becomes '921'
						else:
							i += j+1 # start in the counting after the last number
							break
	
					#row.append(Decimal(string)) # transform the string into a float number with highly precision and add to 'row'
					row.append(float(string))
	
				else:
					i += 1
	
			if row :
				vector2D.append(row) # add the 'row' to the 2D list
	
	return vector2D


def CreateNewDir(dirName, isLocal=True) :
	"""
	Create a new directory in the current folder.
	"""

	if isLocal == True :
		if not os.path.exists(dirName):
			os.mkdir(dirName)
			print("Directory %s created" % dirName)
		else:    
			print("Directory %s already exists" % dirName)
			yn = raw_input("Do you want to overwrite it? ")
			if yn == "yes" or yn == "y" or yn == "Y" or yn == "YES" :
				shutil.rmtree(dirName)
				os.mkdir(dirName)
			else :
				sys.exit(1)
		nameDir = dirName
	else :
		date_time = datetime.datetime.now().strftime("%Y%m%d_%I%M%p")
		nameDir = dirName + "_" + date_time
		if os.path.exists(nameDir):
			shutil.rmtree(nameDir)
		os.mkdir(nameDir)

	path = os.getcwd() + "/"
	os.chdir(path + nameDir)

	return


def InvMat(A) :
	"""
	Invert the 2x2 matrix A.
	Input:
	* A = 2x2 matrix
	Output: 2x2 matrix
	"""

	det = A[0][0]*A[1][1] - A[1][0]*A[0][1]
	Ai_00 = A[1][1]/det
	Ai_01 = -A[0][1]/det
	Ai_10 = -A[1][0]/det
	Ai_11 = A[0][0]/det

	return np.array([[Ai_00, Ai_01], [Ai_10, Ai_11]])


def Dot(A, v) :
	"""
	Dot product between a 2x2 matrix A and a 1x2 vector v.
	Input:
	* A = matrix
	* v = vector
	Output: 1x2 vector
	"""

	vec_0 = A[0][0]*v[0] + A[0][1]*v[1]
	vec_1 = A[1][0]*v[0] + A[1][1]*v[1] 
	vec = [vec_0, vec_1]

	return np.array(vec)


def MatMul(A, B) :
	"""
	Multiply two matrices A and B.
	Input:
	* A = first matrix
	* B = second matrix
	Output: 2x2 matrix
	"""

	C_00 = A[0][0]*B[0][0] + A[0][1]*B[1][0]
	C_01 = A[0][0]*B[0][1] + A[0][1]*B[1][1]
	C_10 = A[1][0]*B[0][0] + A[1][1]*B[1][0]
	C_11 = A[1][0]*B[0][1] + A[1][1]*B[1][1] 

	return np.array([[C_00, C_01], [C_10, C_11]])


def secant(f, x0, x1, eps):
	"""
	Solve nonlinear equation f=0 by secant method.
	x0 is the final point of integration, while x1 is the inital point.
	The iteration continues until ||f|| < eps.
	"""

	f_x0 = f(x0)
	f_x1 = f(x1)
	iteration_counter = 0
	while np.abs(f_x1) > eps and iteration_counter < 100:
		try:
			denominator = float(f_x1 - f_x0)/(x1 - x0)
			x = x1 - float(f_x1)/denominator
		except ZeroDivisionError:
			print "Error! - denominator zero for x = ", x
			sys.exit(1)     # Abort with error
		x0 = x1
		x1 = x
		f_x0 = f_x1
		f_x1 = f(x1)
		iteration_counter += 1
		# Here, either a solution is found, or too many iterations
	if np.abs(f_x1) > eps:
		iteration_counter = -1

	return x, iteration_counter


def Newton(f, dfdx, x, eps):
	"""
	Solve nonlinear equation f=0 by Newton's method.
	dfdx is the first derivative of f with respect to x. Both f and dfdx must be functions of x.
	At input, x holds the start value. The iteration continues until ||f|| < eps.
	"""

	f_value = f(x)
	iteration_counter = 0
	while np.abs(f_value) > eps and iteration_counter < 100:
		try:
			x = x - float(f_value)/dfdx(x)
		except ZeroDivisionError:
			print "Error! - derivative zero for x = ", x
			sys.exit(1)     # Abort with error

		f_value = f(x)
		iteration_counter += 1

	# Here, either a solution is found, or too many iterations
	if np.abs(f_value) > eps:
		iteration_counter = -1
		print "More than 100 iterations...initial condition given!"

	return x, iteration_counter


def Newton_system(F, J, x, eps):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    """
    F_value = F(x)
    F_norm = np.linalg.norm(F_value, ord=2)  # l2 norm of vector
    iteration_counter = 0
    while np.abs(F_norm) > eps and iteration_counter < 100:
    	nF_value = [i * (-1) for i in F_value]
        delta = np.linalg.solve(J(x), nF_value)
        x = x + delta
        F_value = F(x)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1

    # Here, either a solution is found, or too many iterations
    if np.abs(F_norm) > eps:
        iteration_counter = -1
        print "More than 100 iterations...initial condition given!"

    return x, iteration_counter


def Hubble(xpixypiy, mllpxilpp) :
	"""
	Hubble constant. From Wolfram Mathematica.
	Input:
	* xpixypiy = initial condition vector: [x0, pix0, y0, piy0]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: value of Hubble constant
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	H = 1./(np.sqrt(6.0)*np.sqrt((1 + xi*(x**2 + y**2))**2/(4*lp*x*y*(-x**2 + y**2) + m**2*(x**2 + y**2) + (l*(x**2 + y**2)**2)
		/2. + (6*xi**2*(1 + xi*(x**2 + y**2))**4*(x*pix + y*piy)**2)/(1 + xi*(1 + 6*xi)*(x**2 + y**2))**2 + ((1 + xi*
		(x**2 + y**2))**3*((1 + xi*x**2 + xi*(1 + 6*xi)*y**2)*pix - 6*xi**2*x*y*piy)**2)
		/(1 + xi*(1 + 6*xi)*(x**2 + y**2))**2 + ((1 + xi*(x**2 + y**2))**3*(piy + xi*(-6*xi*x*y*pix + ((1 + 6*xi)
		*x**2 + y**2)*piy))**2)/(1 + xi*(1 + 6*xi)*(x**2 + y**2))**2)))

	return H


def xdyd(xpixypiy, mllpxilpp) :
	"""
	Derivatives of X and Y fields with respect to time. From Wolfram Mathematica.
	Input:
	* xpixypiy = list of 4 elements [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: [xd, yd]
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	xd = ((pix + pix*x**2*xi - 6*piy*x*xi**2*y + pix*xi*(1 + 6*xi)*y**2)*(1 + xi*(x**2 + y**2)))/(1 + xi
		*(1 + 6*xi)*(x**2 + y**2))

	yd = ((1 + xi*(x**2 + y**2))*(piy + xi*(x**2*(piy + 6*piy*xi) - 6*pix*x*xi*y + piy*y**2)))/(1 + xi
		*(1 + 6*xi)*(x**2 + y**2))

	return np.array([xd, yd])


def xpyp(xpixypiy, mllpxilpp) :
	"""
	Derivatives of X and Y fields with respect to N. From Wolfram Mathematica.
	Input:
	* xpixypiy = initial condition vector: [x0, pix0, y0, piy0]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: [xp, yp]
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	H = Hubble(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
	xd, yd = xdyd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)

	xp = xd/H
	yp = yd/H

	return np.array([xp, yp])


def pip(pi, x, y, mllpxilpp) :
	"""
	Equations for the derivatives of the canonical momenta pix, piy with respect to N. From Wolfram Mathematica.
	Input:
	* pi = [pix, piy] variables
	* x = initial value of x
	* y = initial value of y
	* mllpxilpp = [m, l, lp, xi, lpp] parameters
	Output: array [pixp, piyp]
	"""

	pix, piy = pi
	m, l, lp, xi, lpp = mllpxilpp

	Hi = 1./Hubble(xpixypiy=[x, pix, y, piy], mllpxilpp=mllpxilpp)

	pixp = -3*pix + Hi*(-((m**2*x + l*x**3 - 6*lp*x**2*y + l*x*y**2 + 2*lp*y**3)/(1 + xi*(x**2 + y**2))**2) + (x*xi
		*(-8*lp*x*(x - y)*y*(x + y) + 2*m**2*(x**2 + y**2) + l*(x**2 + y**2)**2))/(1 + xi*(x**2 + y**2))**3 + (xi
		*(-(piy**2*x*(x**4*xi**2*(1 + 6*xi)**2 + 2*x**2*xi*(1 + 6*xi)*(1 + xi*(1 + 6*xi)*y**2) + (1 + xi*y**2)
		*(1 + xi*(y + 6*xi*y)**2))) + 6*pix*piy*xi*y*(1 + 2*xi*(x**2 + y**2) + 6*xi**3*(x**2 + y**2)**2 + xi**2
		*(x**4 + 6*y**2 + y**4 + 2*x**2*(-3 + y**2))) - pix**2*x*(1 + 2*xi*(-3 + x**2 + y**2) + xi**2*(x**2 + y**2)**2 + 6
		*xi**3*(-6*y**2 + (x**2 + y**2)**2))))/(1 + xi*(1 + 6*xi)*(x**2 + y**2))**2)

	piyp = -3*piy + Hi*((xi*y*(-8*lp*x*(x - y)*y*(x + y) + 2*m**2*(x**2 + y**2) + l*(x**2 + y**2)**2))/(1 + xi
		*(x**2 + y**2))**3 + (2*lp*x*(x**2 - 3*y**2) - y*(m**2 + l*(x**2 + y**2)))/(1 + xi*(x**2 + y**2))**2 + (xi
		*(-(piy**2*y*(1 + 2*xi*(-3 + x**2 + y**2) + xi**2*(x**2 + y**2)**2 + 6*xi**3*(x**4 + y**4 + 2*x**2*(-3 + y**2)))) + 6
		*pix*piy*x*xi*(1 + 2*xi*(x**2 + y**2) + 6*xi**3*(x**2 + y**2)**2 + xi**2*(x**4 - 6*y**2 + y**4 + 2*x**2*(3 + y**2))) - pix**2
		*y*(x**4*xi**2*(1 + 6*xi)**2 + (1 + xi*(1 + 6*xi)*y**2)**2 + 2*x**2*xi*(1 + xi*(6 + 18*xi + (y + 6*xi*y)**2)))))/(1 + xi
		*(1 + 6*xi)*(x**2 + y**2))**2)

	return [pixp, piyp]


def Jpip(pi, x, y, mllpxilpp) :
	"""
	Jacobian of the derivatives of the canonical momenta pix, piy with respect to N. 
	Useful for applying the Netwon method to compute initial pix, piy satisfying Slow-roll. From Wolfram Mathematica.
	Input:
	* pi = [pix, piy] variables
	* x = initial value of x
	* y = initial value of y
	* mllpxilpp = [m, l, lp, xi, lpp] parameters
	Output: matrix [[pixpDpix, pixpDpiy],[piypDpix, piypDpiy]]
	"""

	pix, piy = pi
	m, l, lp, xi, lpp = mllpxilpp

	Hi = 1./Hubble(xpixypiy=[x, pix, y, piy], mllpxilpp=mllpxilpp)

	pixpDpix = -3 + (2*Hi*xi*(-(pix*x**5*xi**2*(1 + 6*xi)) + 3*piy*x**4*xi**3*(1 + 6*xi)*y - 2*pix*x**3*xi*(1 + xi*(1 + 6*xi)*y**2) + 3
		*piy*xi*y*(1 + 2*xi*(1 + 3*xi)*y**2 + xi**2*(1 + 6*xi)*y**4) + 6*piy*x**2*xi**2*y*(1 + 6*xi**2*y**2 + xi*(-3 + y**2)) - pix*x
		*(1 + xi**2*y**4 + 6*xi**3*y**2*(-6 + y**2) + 2*xi*(-3 + y**2))))/(1 + x**2*xi*(1 + 6*xi) + xi*(1 + 6*xi)*y**2)**2

	pixpDpiy = (-2*Hi*xi*(piy*x*(1 + x**2*xi*(1 + 6*xi))**2 - 3*pix*xi*(1 + 2*x**2*(1 - 3*xi)*xi + x**4*xi**2*(1 + 6*xi))*y + 2*piy*x
		*xi*(1 + (6 + x**2)*xi + 6*(3 + 2*x**2)*xi**2 + 36*x**2*xi**3)*y**2 - 6*pix*xi**2*(1 + (3 + x**2)*xi + 6*x**2*xi**2)*y**3 + piy*x
		*xi**2*(1 + 6*xi)**2*y**4 - 3*pix*xi**3*(1 + 6*xi)*y**5))/(1 + x**2*xi*(1 + 6*xi) + xi*(1 + 6*xi)*y**2)**2

	piypDpix = (2*Hi*xi*(-(pix*x**5*xi**2*(1 + 6*xi)) + 3*piy*x**4*xi**3*(1 + 6*xi)*y - 2*pix*x**3*xi*(1 + xi*(1 + 6*xi)*y**2) + 3*piy*xi
		*y*(1 + 2*xi*(1 + 3*xi)*y**2 + xi**2*(1 + 6*xi)*y**4) + 6*piy*x**2*xi**2*y*(1 + 6*xi**2*y**2 + xi*(-3 + y**2)) - pix*x*(1 + xi**2
		*y**4 + 6*xi**3*y**2*(-6 + y**2) + 2*xi*(-3 + y**2))))/(1 + x**2*xi*(1 + 6*xi) + xi*(1 + 6*xi)*y**2)**2

	piypDpiy = -3 - (2*Hi*xi*(piy*x*(1 + x**2*xi*(1 + 6*xi))**2 - 3*pix*xi*(1 + 2*x**2*(1 - 3*xi)*xi + x**4*xi**2*(1 + 6*xi))*y + 2*piy
		*x*xi*(1 + (6 + x**2)*xi + 6*(3 + 2*x**2)*xi**2 + 36*x**2*xi**3)*y**2 - 6*pix*xi**2*(1 + (3 + x**2)*xi + 6*x**2*xi**2)*y**3 + piy
		*x*xi**2*(1 + 6*xi)**2*y**4 - 3*pix*xi**3*(1 + 6*xi)*y**5))/(1 + x**2*xi*(1 + 6*xi) + xi*(1 + 6*xi)*y**2)**2

	return [[pixpDpix, pixpDpiy], [piypDpix, piypDpiy]]


def Solve_pi(xpixypiy0, mllpxilpp, eps=1e-10) :
	"""
	Find the initial values of pix and piy satisfying the slow-roll condition.
	Input:
	* xpixypiy0 = initial condition vector: [x0, pix0, y0, piy0]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* eps = accuracy in Newton's method
	Output: [pix, piy]
	"""

	x0, pix0, y0, piy0 = xpixypiy0
	pi0, numIter = Newton_system(lambda pi: pip(pi, x0, y0, mllpxilpp), lambda pi: Jpip(pi, x0, y0, mllpxilpp), x=[pix0, piy0], eps=eps)
	pix00, piy00 = pi0

	return [pix00, piy00]


def psix_y(x, y, mllpxilpp) :
	"""
	Derivative of psi with repsect to X or Y fields. From Wolfram Mathematica.
	Input:
	* x = value of the field X
	* y = value of the field Y
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: [psix, psiy]
	"""

	m, l, lp, xi, lpp = mllpxilpp

	psix = -(y/((x**2 + y**2 + xi*(x**2 + y**2)**2)*np.sqrt(1 + (6*xi**2*(x**2 + y**2))/(1 + xi*(x**2 + y**2)))))
	psiy = x/((x**2 + y**2 + xi*(x**2 + y**2)**2)*np.sqrt(1 + (6*xi**2*(x**2 + y**2))/(1 + xi*(x**2 + y**2))))

	return [psix, psiy]


def psip(xpixypiy, mllpxilpp) :
	"""
	Differential equation for psi prime with respect to the number of e-folding NN. From Wolfram Mathematica.
	Input:
	* xpixypiy = initial condition vector: [x0, pix0, y0, piy0]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: 1D list of 1 element representing the differential equation for psi prime
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	xp, yp = xpyp(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
	psix, psiy = psix_y(x=x, y=y, mllpxilpp=mllpxilpp)

	psip = psix*xp + psiy*yp

	return psip


def EoMs(NN, xpixypiypsi, mllpxilpp) :
	"""
	Equations of motion for original field x, its conjugate momentum pix, the field y, its conjugate momentum piy 
	in terms of the number of the number of e-folding NN. From Wolfram Mathematica.
	Input:
	* NN = integration variable: number of e-foldings
	* xpixypiy = initial condition vector: [x0, pix0, y0, piy0]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: 1D list of 4 elements representing the system of differential equations: xp, pixp, yp, piyp
	"""

	x, pix, y, piy, psi = xpixypiypsi
	xpixypiy = [x, pix, y, piy]
	m, l, lp, xi, lpp = mllpxilpp

	xp, yp = xpyp(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
	pixp, piyp = pip(pi=[pix, piy], x=x, y=y, mllpxilpp=mllpxilpp)
	psip1 = psip(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)

	return [xp, pixp, yp, piyp, psip1]


def Solve_EoMs(NN0, xpixypiypsi0, mllpxilpp, method='RK45', rtol=1.e-4, atol=1.e-7, dense_output=False, ist_eval=True, stopEnd=False) :
	"""
	Find the solution of the equations of motion for x, pix, y, piy.
	Input:
	* NN0 = array corresponding to the number of e-foldings
	* xpixypiy0 = initial condition vector: [x0, pix0, y0, piy0]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* method = method used to solve the system of differential equations: RK45, RK23, Radau, BDF, LSODA (fixed: 'RK45')
	* dense_output = False (True if we want a callable function in output)
	* rtol = relative accuracy, i.e. number of correct digits (fixed: 1.e-4)
 	* atol = absolute tollerance with local error given by atol+rtol*sol (fixed: 1.e-7)
 	* ist_eval = True when you want the solution to be computed at the point given in an array
 	Output: [NN, x, pix, y, piy]
	"""

	def firstx0(t, xpixypiypsi) : return xpixypiypsi[0]
	
	if stopEnd :
		firstx0.terminal = True
	else :
		firstx0.terminal = False
	firstx0.direction = -1

	if ist_eval :
		sol = integrate.solve_ivp(lambda t, xpixypiypsi: EoMs(t, xpixypiypsi, mllpxilpp), [NN0[0], NN0[-1]], y0=xpixypiypsi0, method=method, t_eval=NN0, dense_output=dense_output, events=firstx0, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	else :
		sol = integrate.solve_ivp(lambda t, xpixypiypsi: EoMs(t, xpixypiypsi, mllpxilpp), [NN0[0], NN0[-1]], y0=xpixypiypsi0, method=method, t_eval=None, dense_output=dense_output, events=firstx0, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	## LSODA fortran method in order to solve the system of linear first order differential equations
	x = sol.y[0, :]
	pix = sol.y[1, :]
	y = sol.y[2, :]
	piy = sol.y[3, :]
	psi1 = sol.y[4, :]

	NN = sol.t # this is identical to NN0
	NNend = sol.t_events

	if dense_output :
		Fsol = sol.sol
		return [Fsol, NNend[0][0]]
	else :
		return [NN, x, pix, y, piy, psi1, NNend[0][0]]


def NB(xpixypiy, mllpxilpp) :
	"""
	Baryon number density. From Wolfram Mathematica.
	Input:
	* xpixypiy = list of 4 elements [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: NB
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	nb = 2*(1 + xi*(x**2 + y**2))*(x*piy - pix*y)
	
	return nb


def NI(xpixypiy, mllpxilpp) :
	"""
	Inflaton number density. From Wolfram Mathematica.
	Input:
	* xpixypiy = list of 4 elements [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: NI
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	H = Hubble(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
	rho = 3*H**2
	ni = rho/m
	
	return ni


def nrh(mllpxilpp) :
	"""
	Density of phi field at reheating in Matter-Dominated universe. Approximate formula.
	Input:
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: nrh
	"""

	m, l, lp, xi, lpp = mllpxilpp

	Gamma_rate = 3*lpp**2*m/(256*np.pi**3)
	n = np.abs(4.0*Gamma_rate**2/(3.0*m))

	return n


def etaB(xpixypiye, mllpxilpp, corr=True) :
	"""
	Baryon-to-entropy ratio at reheating, which should be conserved into the late universe.
	Input:
	* xpixypiye = final condition vector at Nend: [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* corr = boolean if we want to apply the correction factor (fixed to True)
	Output: etaB
	"""

	m, l, lp, xi, lpp = mllpxilpp

	gth = 106.75 + 18.0
	Gamma_rate = 3*lpp**2*m/(256*np.pi**3)
	Trh = np.abs((90.0/(np.pi**2*gth))**0.25*(Gamma_rate)**0.5)
	s = (2.0*np.pi**2/45)*gth*Trh**3

	nbrh = NB(xpixypiy=xpixypiye, mllpxilpp=mllpxilpp)
	nirh = NI(xpixypiy=xpixypiye, mllpxilpp=mllpxilpp)
			
	eta_e = np.abs(nbrh/nirh)
	etab = eta_e*nrh(mllpxilpp=mllpxilpp)/s

	if corr :
		etab = (33.0/111)*etab

	return etab


def Lkin(xpixypiy, mllpxilpp) :
	"""
	Kinetic energy Lkin. From Wolfram Mathematica.
	Input:
	* xpixypiy = [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: Lkin
	"""

	m, l, lp, xi, lpp = mllpxilpp
	x, pix, y, piy = xpixypiy

	xd, yd = xdyd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)

	Lk = 0.5*(xd**2 + yd**2) + 3.0*xi**2*(x*xd + y*yd)**2/(1 + xi*(x**2 + y**2))**2

	return Lk


def VE(x, y, mllpxilpp) :
	"""
	Potential energy VE. From Wolfram Mathematica.
	Input:
	* x = value of the field X
	* y = value of the field Y
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: VE
	"""
	m, l, lp, xi, lpp = mllpxilpp

	Ve = ((m**2*(x**2 + y**2))/2. + (l*(x**2 + y**2)**2)/4. - 2*lp*(x**3*y - x*y**3))/(1 + xi*(x**2 + y**2))**2
	return Ve


def VEx_y(x, y, mllpxilpp) :
	"""
	First derivative of the potential with respect to the field X or Y. From Wolfram Mathematica.
	Input:
	* x = value of the field X
	* y = value of the field Y
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: [VEx, VEy]
	"""
	m, l, lp, xi, lpp = mllpxilpp

	Vex = (l*x**3 + 2*lp*x**4*xi*y + l*x*y**2 + 2*lp*y**3*(1 + xi*y**2) - m**2*x*(-1 + x**2*xi + xi*y**2) - 6*lp
		*x**2*y*(1 + 2*xi*y**2))/(1 + x**2*xi + xi*y**2)**3

	Vey = (-2*lp*x**5*xi + l*x**2*y + l*y**3 - 2*lp*x*y**2*(-3 + xi*y**2) - m**2*y*(-1 + x**2*xi + xi*y**2) + 2
		*lp*x**3*(-1 + 6*xi*y**2))/(1 + x**2*xi + xi*y**2)**3
	
	return np.array([Vex, Vey])


def VExx_xy_yy(x, y, mllpxilpp) :
	"""
	Second derivative of the potential with respect to the field X or Y or mixed. From Wolfram Mathematica.
	Input:
	* x = value of the field X
	* y = value of the field Y
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: matrix [[VExx, VExy], [VEyx, VEyy]]
	"""
	m, l, lp, xi, lpp = mllpxilpp

	Vexx = (-3*l*x**4*xi - 4*lp*x**5*xi**2*y + l*x**2*(3 - 2*xi*y**2) + l*y**2*(1 + xi*y**2) + 8*lp*x**3*xi
		*y*(4 + 7*xi*y**2) -12*lp*x*y*(1 + 4*xi*y**2 + 3*xi**2*y**4) + m**2*(1 + 3*x**4*xi**2 - xi**2*y**4 + 2
		*x**2*xi*(-4 + xi*y**2)))/(1 + x**2*xi + xi*y**2)**4

	Vexy = (2*(lp*x**6*xi**2 + 2*x**3*xi*(-l + m**2*xi)*y - lp*x**4*xi*(2 + 23*xi*y**2) + lp*y**2*(3 + 2*xi
		*y**2 - xi**2*y**4) + lp*x**2*(-3 + 23*xi**2*y**4) + x*y*(l - 2*l*xi*y**2 + 2*m**2*xi*(-2 + xi
		*y**2))))/(1 + x**2*xi + xi*y**2)**4

	Veyx = Vexy

	Veyy = (l*x**4*xi + 36*lp*x**5*xi**2*y - 3*l*y**2*(-1 + xi*y**2) - 8*lp*x**3*xi*y*(-6 + 7*xi*y**2) + x**2
		*(l - 2*l*xi*y**2) + 4*lp*x*y*(3 - 8*xi*y**2 + xi**2*y**4) + m**2*(1 - x**4*xi**2 + 2*xi*(-4 + x**2*xi)
		*y**2 + 3*xi**2*y**4))/(1 + x**2*xi + xi*y**2)**4
	
	return [[Vexx, Vexy], [Veyx, Veyy]]


def Z0(x, y, mllpxilpp) :
	"""
	Z0 matrix. From Wolfram Mathematica.
	Input:
	* x = value of the field X
	* y = value of the field Y
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: 2x2 matrix
	"""

	m, l, lp, xi, lpp = mllpxilpp

	ct = x/np.sqrt(x**2 + y**2)
	st = y/np.sqrt(x**2 + y**2)
	Omega = np.sqrt(1.0/(1 + xi*(x**2 + y**2)))
	e1 = 1.0/(Omega*np.sqrt(1 + 6*Omega**2*xi**2*(x**2 + y**2)))
	e2 = 1.0/Omega

	return np.array([[ct*e1, -st*e2], [st*e1, ct*e2]])


def dZ0(x, y, mllpxilpp) :
	"""
	Derivative of Z0 matrix with respect to X or Y fields. From Wolfram Mathematica.
	Input:
	* x = value of the field X
	* y = value of the field Y
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	Output: [dZ0x, dZ0y]
	"""

	m, l, lp, xi, lpp = mllpxilpp

	dZ0x_00 = (np.sqrt(1/(1 + xi*(x**2 + y**2)))*(y**2 + 6*xi**3*(x**2 + y**2)**3 + xi*(x**2 + y**2)*(x**2 + 2
		*y**2) + xi**2*(x**2 + y**2)*(x**4 + 6*y**2 + y**4 + 2*x**2*(-3 + y**2))))/((x**2 + y**2)**1.5*(1 + xi
		*(1 + 6*xi)*(x**2 + y**2))*np.sqrt(1 + (6*xi**2*(x**2 + y**2))/(1 + xi*(x**2 + y**2))))
	dZ0x_01 = (x*y*np.sqrt(1/(1 + xi*(x**2 + y**2))))/(x**2 + y**2)**1.5
	dZ0x_10 = -((x*y*np.sqrt(1/(1 + xi*(x**2 + y**2)))*(1 + xi*(1 + 12*xi)*(x**2 + y**2)))/((x**2 + y**2)**1.5
		*(1 + xi*(1 + 6*xi)*(x**2 + y**2))*np.sqrt(1 + (6*xi**2*(x**2 + y**2))/(1 + xi*(x**2 + y**2)))))
	dZ0x_11 = (np.sqrt(1/(1 + xi*(x**2 + y**2)))*(y**2 + xi*(x**2 + y**2)**2))/(x**2 + y**2)**1.5
	
	dZ0x = np.array([[dZ0x_00, dZ0x_01], [dZ0x_10, dZ0x_11]])

	dZ0y_00 = -((x*y*np.sqrt(1/(1 + xi*(x**2 + y**2)))*(1 + xi*(1 + 12*xi)*(x**2 + y**2)))/((x**2 + y**2)**1.5
		*(1 + xi*(1 + 6*xi)*(x**2 + y**2))*np.sqrt(1 + (6*xi**2*(x**2 + y**2))/(1 + xi*(x**2 + y**2)))))
	dZ0y_01 = -((np.sqrt(1/(1 + xi*(x**2 + y**2)))*(x**2 + xi*(x**2 + y**2)**2))/(x**2 + y**2)**1.5)
	dZ0y_10 = (np.sqrt(1/(1 + xi*(x**2 + y**2)))*(x**2 + 6*xi**3*(x**2 + y**2)**3 + xi*(x**2 + y**2)
		*(2*x**2 + y**2) + xi**2*(x**2 + y**2)*(x**4 - 6*y**2 + y**4 + 2*x**2*(3 + y**2))))/((x**2 + y**2)**1.5
		*(1 + xi*(1 + 6*xi)*(x**2 + y**2))*np.sqrt(1 + (6*xi**2*(x**2 + y**2))/(1 + xi*(x**2 + y**2))))
	dZ0y_11 = -((x*y*np.sqrt(1/(1 + xi*(x**2 + y**2))))/(x**2 + y**2)**1.5)

	dZ0y = np.array([[dZ0y_00, dZ0y_01], [dZ0y_10, dZ0y_11]])

	return [dZ0x, dZ0y]


def Z(x, y, mllpxilpp, psi=0) :
	"""
	Matrix Z. From Wolfram Mathematica.
	Input:
	* x = value of the field X
	* y = value of the field Y
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi (optional: fixed psi = 0)
	Output: 2x2 matrix Z
	"""

	Z0m = Z0(x=x, y=y, mllpxilpp=mllpxilpp)

	sp = np.sin(psi)
	cp = np.cos(psi)
	Rpsi = [[cp, sp], [-sp, cp]]

	Zm = MatMul(A=Z0m, B=Rpsi) 

	return Zm


def dZ(x, y, mllpxilpp, psi=0) :
	"""
	Derivative of Z matrix with respect to X or Y fields.
	Input:
	* x = value of the field X
	* y = value of the field Y
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* * psi = angle psi (optional: fixed psi = 0)
	Output: [dZx, dZy]
	"""

	m, l, lp, xi, lpp = mllpxilpp

	sp = np.sin(psi)
	cp = np.cos(psi)
	
	dZ0x, dZ0y = dZ0(x=x, y=y, mllpxilpp=mllpxilpp)
	Rpsi = np.array([[cp, sp], [-sp, cp]])

	Z0m = Z0(x=x, y=y, mllpxilpp=mllpxilpp)

	psix, psiy = psix_y(x=x, y=y, mllpxilpp=mllpxilpp)
	dRpsi = np.array([[-sp, cp], [-cp, -sp]])

	dZx = MatMul(A=dZ0x, B=Rpsi) + psix*MatMul(A=Z0m, B=dRpsi)
	dZy = MatMul(A=dZ0y, B=Rpsi) + psiy*MatMul(A=Z0m, B=dRpsi)
	
	return [dZx, dZy]


def UdVd(xpixypiy, mllpxilpp, psi=0) :
	"""
	Time-derivative of the fields U and V.
	Input:
	* xpixypiy = list of 4 elements [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi (optional: fixed psi = 0)
	Output: array([Ud, Vd])
	"""

	x, pix, y, piy = xpixypiy
	
	xdyd1 = xdyd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
	Zm = Z(x=x, y=y, mllpxilpp=mllpxilpp, psi=psi)
	Zi = InvMat(A=Zm)

	return Dot(A=Zi, v=xdyd1)


def Compute_theta(xpixypiy, NN, mllpxilpp, psi=0) :
	"""
	Compute the angle theta (elsewhere is also called phi!) between Vd and Ud and its time-derivative.
	Input:
	* xpixypiy = list of 4 elements [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi (optional: fixed psi = 0)
	Output: [Ud, Vd, theta, costh, sinth, thd, thd_alt]
	"""

	x, pix, y, piy = xpixypiy

	Ud, Vd = UdVd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)
	tanth = Vd/Ud

	tan_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=tanth)
	Ud_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=Ud)
	Vd_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=Vd)

	tanp_tmp = tan_interp.derivative(n=1)
	tanp = tanp_tmp(NN)
	Udp_tmp = Ud_interp.derivative(n=1)
	Udp = Udp_tmp(NN)
	Vdp_tmp = Vd_interp.derivative(n=1)
	Vdp = Vdp_tmp(NN)

	thd = Hubble(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)*tanp/(1+tanth**2)
	thd_alt = Hubble(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)*(Ud*Vdp - Vd*Udp)/(Ud*Ud + Vd*Vd)
	theta = np.arctan(tanth)
	costh = np.cos(theta) # equivalent to cphi below
	sinth = np.sin(theta) # equivalent to sphi below

	return [Ud, Vd, theta, costh, sinth, thd, thd_alt]


def Compute_gamma(xpixypiy, mllpxilpp, psi=0) :
	"""
	Compute gamma, which is the angle between (dB,dC) and (dU,dV). 
	Input:
	* xpixypiy = list of 4 elements [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi (optional: fixed psi = 0)
	Output:
	"""

	x, pix, y, piy = xpixypiy

	Vex, Vey = VEx_y(x=x, y=y, mllpxilpp=mllpxilpp)
	Zm = Z(x=x, y=y, mllpxilpp=mllpxilpp, psi=psi)

	Ud, Vd = UdVd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)

	# Matteo
	VeU = Zm[0][0]*Vex + Zm[1][0]*Vey
	VeV = Zm[0][1]*Vex + Zm[1][1]*Vey

	#idx_Matteo = np.where(VeV == 0.0)[0]
	#tang_Matteo = -VeU/VeV 	# approximate formula because slow-roll approx used
	idx_Matteo = np.where(Vd == 0.0)[0]
	tang_Matteo = - Ud/Vd
	cosg_Matteo = np.cos(np.arctan(tang_Matteo))
	sing_Matteo = np.sin(np.arctan(tang_Matteo))

	for i in range(0, len(idx_Matteo)) :
		tang_Matteo[i] = float('nan')
		cosg_Matteo[i] = 0.0
		sing_Matteo[i] = 1.0 
		
	# Jim
	Num = -Zm[0][1]*Vey + Zm[1][1]*Vex
	Den = -Zm[0][0]*Vey + Zm[1][0]*Vex

	idx_Jim = np.where(Den == 0.0)[0]
	tang_Jim = Num/Den 	# perhaps it is correct
	cosg_Jim = np.cos(np.arctan(tang_Jim))
	sing_Jim = np.sin(np.arctan(tang_Jim))

	for i in range(0, len(idx_Jim)) :
		tang_Jim[i] = float('nan')
		cosg_Jim[i] = 0.0
		sing_Jim[i] = 1.0 

	return [tang_Matteo, cosg_Matteo, sing_Matteo, tang_Jim, cosg_Jim, sing_Jim]


def SR_params_UV(xpixypiy, mllpxilpp, psi=0) :
	"""
	Compute the slow-roll parameters epsilon and eta in terms of U and V fields.
	Input:
	* xpixypiy = list of 4 elements [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi (optional: fixed psi = 0)
	Output: [Ve, VeU, VeV, epsU, epsV, etaUU, etaUV, etaVU, etaVV]
	"""

	x, pix, y, piy = xpixypiy

	#Ve = VE(x=x, y=y, mllpxilpp=mllpxilpp)
	H = Hubble(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
	Ve = 3*H**2 # this makes the transfer function eqs well defined beyond slow roll
	
	Vex, Vey = VEx_y(x=x, y=y, mllpxilpp=mllpxilpp)

	Ve2 = VExx_xy_yy(x=x, y=y, mllpxilpp=mllpxilpp)
	Vexx = Ve2[0][0]
	Veyy = Ve2[1][1]
	Vexy = Ve2[0][1]

	Zm = Z(x=x, y=y, mllpxilpp=mllpxilpp, psi=psi) # approximately Rpsi = 1 since we are free to choose psi = 0 (no its first derivative)
	dZx, dZy = dZ(x=x, y=y, mllpxilpp=mllpxilpp, psi=psi)

	VeU = Zm[0][0]*Vex + Zm[1][0]*Vey
	VeV = Zm[0][1]*Vex + Zm[1][1]*Vey

	epsU = (VeU)**2 /(2*Ve**2)
	epsV = (VeV)**2 /(2*Ve**2)
	etaUU = (Zm[0][0]**2*Vexx + 2*Zm[0][0]*Zm[1][0]*Vexy + Zm[1][0]**2*Veyy + Zm[0][0]*dZx[0][0]*Vex + Zm[0][0]*dZx[1][0]*Vey + Zm[1][0]*dZy[0][0]*Vex + Zm[1][0]*dZy[1][0]*Vey)/Ve
	etaVV = (Zm[0][1]**2*Vexx + 2*Zm[0][1]*Zm[1][1]*Vexy + Zm[1][1]**2*Veyy + Zm[0][1]*dZx[0][1]*Vex + Zm[0][1]*dZx[1][1]*Vey + Zm[1][1]*dZy[0][1]*Vex + Zm[1][1]*dZy[1][1]*Vey)/Ve
	etaUV = (Zm[0][0]*Zm[0][1]*Vexx + (Zm[0][0]*Zm[1][1] + Zm[1][0]*Zm[0][1])*Vexy + Zm[1][0]*Zm[1][1]*Veyy + Zm[0][0]*dZx[0][1]*Vex + Zm[0][0]*dZx[1][1]*Vey + Zm[1][0]*dZy[0][1]*Vex + Zm[1][0]*dZy[1][1]*Vey)/Ve
	etaVU = (Zm[0][1]*Zm[0][0]*Vexx + (Zm[0][1]*Zm[1][0] + Zm[1][1]*Zm[0][0])*Vexy + Zm[1][1]*Zm[1][0]*Veyy + Zm[0][1]*dZx[0][0]*Vex + Zm[0][1]*dZx[1][0]*Vey + Zm[1][1]*dZy[0][0]*Vex + Zm[1][1]*dZy[1][0]*Vey)/Ve

	return [Ve, VeU, VeV, epsU, epsV, etaUU, etaUV, etaVU, etaVV]


def GetNend(xpixypiy, NN, mllpxilpp, psi, Nsteps=1e4) :
	"""
	Find the number of e-folding at which inflation ends, namely when the real scalar field X reaches zero.
	Use interpolation technique for better precision
	Input:
	* xpixypiy = list of 4 elements containing the values of x, pix, y, piy at some NN_i
	* NN = number of e-foldings
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi
	* Nsteps = number of steps in the interpolation procedure to get a better precion (fixed)
	Output: [iend, Nend, xend, pixend, yend, piyend, psiend, rhoend]
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	xx = interpolate.InterpolatedUnivariateSpline(x=NN, y=x)
	pixpix = interpolate.InterpolatedUnivariateSpline(x=NN, y=pix)
	yy = interpolate.InterpolatedUnivariateSpline(x=NN, y=y)
	piypiy = interpolate.InterpolatedUnivariateSpline(x=NN, y=piy)
	psipsi = interpolate.InterpolatedUnivariateSpline(x=NN, y=psi)

	for i in range(0, len(NN)) :
		if x[i] <= 0.0 :
			#print "i = %d, NN_i = %f, x_i = %.8e" % (i, NN[i], x[i])
			iend = i
			if i < len(NN)-1 :
				NNinterp = np.linspace(start=NN[i-1], stop=NN[i+1], num=int(Nsteps))
				xinterp = xx(NNinterp)
				for j in range(0, len(NNinterp)) :
					xinterp_j = xinterp[j]
					if xinterp_j <= 0.0 :
						Nend = NNinterp[j]
						break
						
			else :
				NNinterp = np.linspace(start=NN[i-1], stop=NN[i], num=int(Nsteps))
				xinterp = xx(NNinterp)
				for j in range(0, len(NNinterp)) :
					xinterp_j = xinterp[j]
					if xinterp_j <= 0.0 :
						Nend = NNinterp[j]
						break
			break

	try:
		Nend
	except NameError:
		print "Nend is bigger than Ntot!"
		sys.exit(1)

	xend = xx(Nend)
	pixend = pixpix(Nend)
	yend = yy(Nend)
	piyend = piypiy(Nend)
	psiend = psipsi(Nend)
	rhoend = m*NI(xpixypiy=[xend, pixend, yend, piyend], mllpxilpp=mllpxilpp)

	return [iend, Nend, xend, pixend, yend, piyend, psiend, rhoend]


def nsrOLD(xpixypiy, NN, mllpxilpp, psi, Nend) :
	"""
	Values of ns and r defined as in Bartolo et al. paper, https://journals.aps.org/prd/pdf/10.1103/PhysRevD.64.123504
	Input:
	* xpixypiy = list of 4 elements containing the values of x, pix, y, piy at some NN_i
	* NN = number of e-foldings
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi
	* Nend = number of e-folding when inflation ends
	Output: [ns, r]
	"""

	Ud, Vd = UdVd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)
	Ve, VeU, VeV, epsU, epsV, etaUU, etaUV, etaVU, etaVV = SR_params_UV(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)

	sU = np.sign(Ud)
	sV = np.sign(Vd)
	epstot = epsU + epsV

	ns = []
	r = []
	for i in range(0, len(NN)) :
		Nk = Nend - NN[i]

		# use definitions in https://journals.aps.org/prd/pdf/10.1103/PhysRevD.64.123504
		tmp1 = Nk*(-(epsU[i]*etaUU[i] + epsV[i]*etaVV[i])/epstot[i] - 2*(sU[i]*np.sqrt(epsU[i]))*(sV[i]*np.sqrt(epsV[i]))*etaUV[i]/epstot[i] + 2*epstot[i])
		f = bigfloat.exp(tmp1)
		tmp2 = Nk*(-(epsU[i]*etaVV[i] + epsV[i]*etaUU[i])/epstot[i] + 2*(sU[i]*np.sqrt(epsU[i]))*(sV[i]*np.sqrt(epsV[i]))*etaUV[i]/epstot[i])
		g = bigfloat.exp(tmp2)
		C = (epsU[i] - epsV[i])*(etaVV[i] - etaUU[i])/epstot[i] - 4*(sU[i]*np.sqrt(epsU[i]))*(sV[i]*np.sqrt(epsV[i]))*etaUV[i]/epstot[i] + 2*epstot[i]
		bdoh = (epsU[i] - epsV[i])*etaUV[i]/epstot[i] + (sU[i]*np.sqrt(epsU[i]))*(sV[i]*np.sqrt(epsV[i]))*(etaUU[i] - etaVV[i])/epstot[i]
		expCN = bigfloat.exp(Nk*C)
		Ptilde = 2*bdoh*g*(expCN - 1)/C
		# in the paper, the scalar (or curvature) spectral index ns = nR
		ns_i = 1.0 - 6.0*epstot[i] + 2.0*(epsU[i]*etaUU[i] + epsV[i]*etaVV[i])/epstot[i] + 4.0*(sU[i]*np.sqrt(epsU[i]))*(sV[i]*np.sqrt(epsV[i]))*etaUV[i]/epstot[i] - 8.0\
			*bdoh**2*(1.0 - 1.0/expCN)/(expCN*C*(1 + (Ptilde/f)**2))
		r_i = 16.0*epstot[i]
		ns.append(ns_i)
		r.append(r_i)

	return [np.array(ns), np.array(r)]


def SR_params_AS(xpixypiy, NN, mllpxilpp, psi, altDefVeS=True) :
	"""
	Compute the slow-roll parameters epsilon and eta in terms of A (or sigma) and S fields.
	Input:
	* xpixypiy = list of 4 elements containing the values of x, pix, y, piy at some NN_i
	* NN = number of e-foldings
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi
	* altDefVeS = boolean: True if we want to use alternative definition of Ves (fixed: True)
	Output: [Ad, VeA, VeS_def, VeS, VeS_alt, epsA, epsS, etaAA, etaAS, etaSS]
	"""

	if altDefVeS :
		Ud, Vd, theta, costh, sinth, thd, thd_alt = Compute_theta(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi)
	else :
		Ud, Vd = UdVd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)

	Ve, VeU, VeV, epsU, epsV, etaUU, etaUV, etaVU, etaVV = SR_params_UV(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)

	Ad = np.sqrt(Ud**2 + Vd**2) 	# dot(sigma)

	cphi = Ud/Ad
	sphi = Vd/Ad

	VeA = cphi*VeU + sphi*VeV
	VeS_def = -sphi*VeU + cphi*VeV

	if altDefVeS :
		VeS = -Ad*thd
		VeS_alt = -Ad*thd_alt

	epsA = VeA**2 / (2.0*Ve**2)
	epsS = VeS_def**2 / (2.0*Ve**2)
	etaAA = cphi**2*etaUU + sphi**2*etaVV + cphi*sphi*(etaUV + etaVU)
	etaAS = cphi*sphi*(etaVV - etaUU) + etaUV*(cphi**2 - sphi**2)
	etaSS = sphi**2*etaUU + cphi**2*etaVV - cphi*sphi*(etaUV + etaVU)

	if altDefVeS :
		return [Ad, VeA, VeS_def, VeS, VeS_alt, epsA, epsS, etaAA, etaAS, etaSS]
	else :
		return [Ad, VeA, VeS_def, epsA, epsS, etaAA, etaAS, etaSS]


def nsrNEW(xpixypiy, NN, mllpxilpp, psi, altDefVeS=False) :
	"""
	Values of ns and r defined as in Gordon et al. paper, https://arxiv.org/pdf/astro-ph/0009131.pdf,
	and in Byrnes and Wands' paper, https://arxiv.org/pdf/astro-ph/0605679.pdf.
	Input:
	* xpixypiy = list of 4 elements containing the values of x, pix, y, piy at some NN_i
	* NN = number of e-foldings
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi
	* altDefVeS = boolean: True if we want to use alternative definition of Ves (fixed: False)
	Output: [As_Matteo, ns, r]
	"""

	x, pix, y, piy = xpixypiy

	if altDefVeS :
		Ad, VeA, VeS_def, VeS, VeS_alt, epsA, epsS, etaAA, etaAS, etaSS = SR_params_AS(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi)
	else :
		Ad, VeA, VeS_def, epsA, epsS, etaAA, etaAS, etaSS = SR_params_AS(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi, altDefVeS=False)

	Cem = 2.0 - np.log(2.0) - 0.5772

	cosD = -2*Cem*etaAS 	# at Hubble exit
	sinD = np.sqrt(np.abs(1 - cosD**2))
	ns = 1.0 - (6.0 - 4*cosD**2)*epsA + 2.0*sinD**2*etaAA + 4.0*sinD*cosD*etaAS + 2.0*cosD**2*etaSS
	r = 16.0*epsA 

	H = Hubble(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
	a1 = -2.0*epsA + 2*Cem*(3*epsA - etaAA)

	P0s_Matteo = H**4 / (2*np.pi*Ad)**2
	As_Matteo = P0s_Matteo*(1 + a1)

	#Ve = VE(x=x, y=y, mllpxilpp=mllpxilpp) # same as Matteo's one!
	#P0s_Jim = Ve/(24*np.pi**2*epsA)
	#As_Jim = P0s_Jim*(1 + a1)

	return [As_Matteo, ns, r]


def GetNstar(xpixypiy, NN, mllpxilpp, psi, As, Nsteps=1e4, eps=1.e-10) :
	"""
	Find the number of e-folding at which the mode kstar exit the Hubble horizon. 
	Formula (47) from https://arxiv.org/pdf/1807.06211.pdf
	Input:
	* xpixypiy = list of 4 elements containing the values of x, pix, y, piy
	* NN = number of e-foldings
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi
	* As = amplitude of the scalar power spectrum (independent of pivot scale)
	* Nsteps = number of steps in the interpolation procedure to get a better precion (fixed)
	Output: value of Nstar
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	iend, Nend, xend, pixend, yend, piyend, psiend, rhoend = GetNend(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi, Nsteps=Nsteps)

	kstar_ns = 0.05 # Mpc^-1, Planck pivot scale for As, ns
	kstar_r = 0.002 # Mpc^-1, Planck pivot scale for r
	Asobs = np.exp(3.044)*1.0e-10

	a0H0 = 67.8/3e5 #a0H0 = 67.66/2.99792458e5 # Mpc^-1
	rhoth = m*nrh(mllpxilpp=mllpxilpp)
	gth = 106.75 + 18.0
	cte = 67.0 - np.log(kstar_ns/a0H0) + (1.0/12)*np.log(rhoth/gth)	# wint = 0

	rescale = Asobs/As
	Ve = VE(x=x, y=y, mllpxilpp=mllpxilpp)

	Ve_res = np.abs(rescale*Ve)
	rhoend_res = np.abs(rescale*rhoend)
	fns = cte - (1.0/3)*np.log(rhoend_res) + 0.5*np.log(Ve_res) - NN # main equation
	fr = fns + np.log(kstar_ns/kstar_r)

	fns_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=fns) # cubic interpolation to make the equation a callable function
	Nstar_ns, numIter_ns = secant(fns_interp, x0=60.0, x1=50.0, eps=eps) # solve the main equation by secant method
	ustar_ns = Nend - Nstar_ns
	fr_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=fr)
	Nstar_r, numIter_r = secant(fr_interp, x0=60.0, x1=50.0, eps=eps)
	ustar_r = Nend - Nstar_r

	
	rescale2 = interpolate.InterpolatedUnivariateSpline(x=NN, y=rescale)
	xx = interpolate.InterpolatedUnivariateSpline(x=NN, y=x)
	pixpix = interpolate.InterpolatedUnivariateSpline(x=NN, y=pix)
	yy = interpolate.InterpolatedUnivariateSpline(x=NN, y=y)
	piypiy = interpolate.InterpolatedUnivariateSpline(x=NN, y=piy)
	psipsi = interpolate.InterpolatedUnivariateSpline(x=NN, y=psi)

	# for As, ns quantities:
	xstar_ns = xx(ustar_ns) # Nstar is defined as the number of e-folding before the end of inflation, but we have to consider ustar = Nend - Nstar as the corresponding point in our range NN = 0,..,Nend
	pixstar_ns = pixpix(ustar_ns)
	ystar_ns = yy(ustar_ns)
	piystar_ns = piypiy(ustar_ns)
	psistar_ns = psipsi(ustar_ns)
	rescalestar_ns = rescale2(ustar_ns)
	star_ns = [Nstar_ns, xstar_ns, pixstar_ns, ystar_ns, piystar_ns, psistar_ns, rescalestar_ns]

	# for r quantities:
	xstar_r = xx(ustar_r) # Nstar is defined as the number of e-folding before the end of inflation, but we consider ustar = Nend - Nstar as the corresponding point in our range NN = 0,..,Nend
	pixstar_r = pixpix(ustar_r)
	ystar_r = yy(ustar_r)
	piystar_r = piypiy(ustar_r)
	psistar_r = psipsi(ustar_r)
	rescalestar_r = rescale2(ustar_r)
	Vstar_r = VE(x=xstar_r, y=ystar_r, mllpxilpp=mllpxilpp)
	star_r = [Nstar_r, xstar_r, pixstar_r, ystar_r, piystar_r, psistar_r, rescalestar_r]

	return [star_ns, star_r]


def Adjust_nsr(NN, ns, r, Nshift, Nsteps=1e4) :
	"""
	Shift the values of r(NN) by Nshift. This is needed in order to align ns with r, which are defined at two different pivot scales for Planck.
	Input:
	* NN = number of e-foldings
	* ns = array of ns values
	* r = array of corresponding r values
	* Nshift = Nstar_r - Nstar_ns (computed at two different pivot scales)
	* Nsteps = number of steps in the interpolation procedure to get a better precion (fixed)
	Output: [ns_shift, r_shift]
	"""

	ns_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=ns, k=3)
	r_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=r, k=3)

	NNshift = np.linspace(start=NN[0]-Nshift, stop=NN[-1]-Nshift, num=Nsteps)

	ns_shift = ns_interp(NN)
	r_shift = r_interp(NNshift)

	return [ns_shift, r_shift]


def Up2pVp2pUpVppC1(xpixypiy, NN, mllpxilpp, psi) :
	"""
	Auxiliary variables entering the second-order differential equation for dU and dV, which are the perturbations of U and V fields.
	This equation corresponds to eq.(5) in https://arxiv.org/pdf/astro-ph/0605679.pdf.
	Input:
	* xpixypiy = list of 4 elements containing the values of x, pix, y, piy
	* NN = number of e-foldings
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi 
	Output: [Up2p, Vp2p, UpVpp, C1]
	"""

	x, pix, y, piy = xpixypiy
	m, l, lp, xi, lpp = mllpxilpp

	Ud, Vd = UdVd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)
	Ad = np.sqrt(Ud**2 + Vd**2) 
	H = Hubble(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
	
	Ap = Ad/H
	Up = Ud/H
	Vp = Vd/H

	Up_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=Up, k=3)
	Upp_tmp = Up_interp.derivative(n=1)
	Upp = Upp_tmp(NN) 
	Vp_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=Vp, k=3)
	Vpp_tmp = Vp_interp.derivative(n=1)
	Vpp = Vpp_tmp(NN) 

	Up2p = 2*Up*Upp
	Vp2p = 2*Vp*Vpp
	UpVpp = Vp*Upp + Up*Vpp
	Hp = -H/2 * Ap**2
	C1 = 3.0 + Hp/H

	return [Up, Vp, Up2p, Vp2p, UpVpp, C1]


def Solve_dUdV(xpixypiy, NN, mllpxilpp, psi, ustar_ns, Nsteps=1e4, method='RK45', rtol=1.e-4, atol=1.e-7) :
	"""
	Find the solution of the system of 2nd order differential equations for dU and dV.
	Input:
	* xpixypiy = [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi
	* ustar_ns = Nend - Nstar_ns
	* Nsteps = number of steps in the interpolation procedure to get a better precion (fixed)
	* method = method used to solve the system of differential equations: RK45, RK23, Radau, BDF, LSODA (fixed: 'RK45')
	* rtol = relative accuracy, i.e. number of correct digits (fixed: 1.e-4)
 	* atol = absolute tollerance with local error given by atol+rtol*sol (fixed: 1.e-7)
 	Output: [NN_ustartot, dU, dUp, dV, dVp, T_, T_] x2
	"""

	x, pix, y, piy = xpixypiy
	Ntot = NN[-1]
	
	Up_tmp, Vp_tmp, Up2p_tmp, Vp2p_tmp, UpVpp_tmp, C1_tmp = Up2pVp2pUpVppC1(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi)
	Ve, VeU, VeV, epsU, epsV, etaUU_tmp, etaUV_tmp, etaVU_tmp, etaVV_tmp = SR_params_UV(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)
	
	Up = interpolate.InterpolatedUnivariateSpline(x=NN, y=Up_tmp, k=3)
	Vp = interpolate.InterpolatedUnivariateSpline(x=NN, y=Vp_tmp, k=3)
	Up2p = interpolate.InterpolatedUnivariateSpline(x=NN, y=Up2p_tmp, k=3)
	Vp2p = interpolate.InterpolatedUnivariateSpline(x=NN, y=Vp2p_tmp, k=3)
	UpVpp = interpolate.InterpolatedUnivariateSpline(x=NN, y=UpVpp_tmp, k=3)
	C1 = interpolate.InterpolatedUnivariateSpline(x=NN, y=C1_tmp, k=3)

	etaUU = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaUU_tmp, k=3)
	etaUV = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaUV_tmp, k=3)
	etaVU = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaVU_tmp, k=3)
	etaVV = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaVV_tmp, k=3)  


	def dUppdVpp(NNi, dUdUpdVdVp) :
		"""
		Equations of motion for dU and dV as given by eq.(5) in https://arxiv.org/pdf/astro-ph/0605679.pdf.
		Input:
		* NNi = integration variable: number of e-foldings
		* dUdUpdVdVp = [dU, dUp, dV, dVp] which are the variables we want to compute
		Output: [dUpp, dVpp]
		"""

		dU, dUp, dV, dVp = dUdUpdVdVp

		dC1 = C1(NNi)*(Up(NNi)*dU + Vp(NNi)*dV)

		dUpp = - C1(NNi)*dUp - 3.0*(etaUU(NNi)*dU + etaUV(NNi)*dV) + Up(NNi)*dC1 + Up2p(NNi)*dU + UpVpp(NNi)*dV
		dVpp = - C1(NNi)*dVp - 3.0*(etaVV(NNi)*dV + etaVU(NNi)*dU) + Vp(NNi)*dC1 + Vp2p(NNi)*dV + UpVpp(NNi)*dU

		return [dUp, dUpp, dVp, dVpp]


	NN0_ustartot = np.linspace(start=ustar_ns, stop=Ntot, num=Nsteps)

	
	## find the initial values of dUp and dVp at ustar = Nend - Nstar such that the condition dUpp = dVpp = 0
	#--- dU* = 1.0, dV* = 0.0 ---#
	dU0_10 = 1.0
	dV0_10 = 0.0
	dUdUpdVdVp0_10 = [dU0_10, 0.0, dV0_10, 0.0]
	dUp_10_ustar, dUpp_10_ustar, dVp_10_ustar, dVpp_10_ustar = dUppdVpp(NNi=ustar_ns, dUdUpdVdVp=dUdUpdVdVp0_10)

	if C1(ustar_ns) != 0.0 :
		dUp0_10 = dUpp_10_ustar/C1(ustar_ns)  # this is because if we set dUp = dVp = 0 and we consider dUpp/C1 and dVpp/C1 we get exactly...
		dVp0_10 = dVpp_10_ustar/C1(ustar_ns)  # ...the same situation if we set dUpp = dVpp = 0 with dUp != 0 and dVp != 0 (look at the above equations!)
	else :
		print "C1(ustar_ns) = 0!"
		sys.exit(1)

	print "[dU0_10, dUp0_10, dV0_10, dVp0_10] = [%.8e, %.8e, %.8e, %.8e]" %(Decimal(dU0_10), Decimal(dUp0_10), Decimal(dV0_10), Decimal(dVp0_10))

	sol_10 = integrate.solve_ivp(lambda NNi, dUdUpdVdVp: dUppdVpp(NNi, dUdUpdVdVp), [ustar_ns, Ntot], [dU0_10, dUp0_10, dV0_10, dVp0_10], method=method, dense_output=False, t_eval=NN0_ustartot, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	dU_10 = sol_10.y[0, :]
	dUp_10 = sol_10.y[1, :]
	dV_10 = sol_10.y[2, :]
	dVp_10 = sol_10.y[3, :]
	NN_ustartot_10 = sol_10.t # it should be equal to NN0_ustartot
	
	#--- dU* = 0.0, dV* = 1.0 ---#
	dU0_01 = 0.0
	dV0_01 = 1.0
	dUdUpdVdVp0_01 = [dU0_01, 0.0, dV0_01, 0.0]
	dUp_01_ustar, dUpp_01_ustar, dVp_01_ustar, dVpp_01_ustar = dUppdVpp(NNi=ustar_ns, dUdUpdVdVp=dUdUpdVdVp0_01)
	
	if C1(ustar_ns) != 0.0 :
		dUp0_01 = dUpp_01_ustar/C1(ustar_ns)
		dVp0_01 = dVpp_01_ustar/C1(ustar_ns)
	else :
		print "C1(ustar_ns) = 0!"
		sys.exit(1)

	print "[dU0_01, dUp0_01, dV0_01, dVp0_01] = [%.8e, %.8e, %.8e, %.8e]" %(Decimal(dU0_01), Decimal(dUp0_01), Decimal(dV0_01), Decimal(dVp0_01))

	sol_01 = integrate.solve_ivp(lambda NNi, dUdUpdVdVp: dUppdVpp(NNi, dUdUpdVdVp), [ustar_ns, Ntot], [dU0_01, dUp0_01, dV0_01, dVp0_01], method=method, dense_output=False, t_eval=NN0_ustartot, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	dU_01 = sol_01.y[0, :]
	dUp_01 = sol_01.y[1, :]
	dV_01 = sol_01.y[2, :]
	dVp_01 = sol_01.y[3, :]
	NN_ustartot_01 = sol_01.t # it should be equal to NN0_ustartot

	
	if np.array_equal(NN_ustartot_10, NN_ustartot_01) :
		NN_ustartot = NN_ustartot_10
	
	xx = interpolate.InterpolatedUnivariateSpline(x=NN, y=x, k=3)
	pixpix = interpolate.InterpolatedUnivariateSpline(x=NN, y=pix, k=3)
	yy = interpolate.InterpolatedUnivariateSpline(x=NN, y=y, k=3)
	piypiy = interpolate.InterpolatedUnivariateSpline(x=NN, y=piy, k=3)
	psipsi = interpolate.InterpolatedUnivariateSpline(x=NN, y=psi, k=3)
	x_ustartot = xx(NN_ustartot)
	pix_ustartot = pixpix(NN_ustartot)
	y_ustartot = yy(NN_ustartot)
	piy_ustartot = piypiy(NN_ustartot)
	psi_ustartot = psipsi(NN_ustartot)
	xpixypiy_ustartot = [x_ustartot, pix_ustartot, y_ustartot, piy_ustartot]

	Ud, Vd = UdVd(xpixypiy=xpixypiy_ustartot, mllpxilpp=mllpxilpp, psi=psi_ustartot)
	Ad = np.sqrt(Ud**2 + Vd**2) 	# dot(sigma)
	H = Hubble(xpixypiy=xpixypiy_ustartot, mllpxilpp=mllpxilpp)

	
	#--- dU* = 1.0, dV* = 0.0 ---#
	uu_10 = H/Ad * dU_10
	vv_10 = H/Ad * dV_10
	T_UU = uu_10/uu_10[0]
	T_VU = vv_10/uu_10[0]
 
	dUdV_10 = [NN_ustartot, dU_10, dUp_10, dV_10, dVp_10, T_UU, T_VU] 

	#--- dU* = 0.0, dV* = 1.0 ---#
	uu_01 = H/Ad * dU_01
	vv_01 = H/Ad * dV_01
	T_UV = uu_01/vv_01[0]
	T_VV = vv_01/vv_01[0]
 
	dUdV_01 = [NN_ustartot, dU_01, dUp_01, dV_01, dVp_01, T_UV, T_VV]

	return [dUdV_10, dUdV_01]


def Solve_dAdS(xpixypiy, NN, mllpxilpp, psi, ustar_ns, Nsteps=1e4, method='RK45', rtol=1.e-4, atol=1.e-7, dense_output=False, events=None, ist_eval=True) :
	"""
	Find the solution of the system of 2nd order differential equations for dA and dS.
	Input:
	* xpixypiy = [x, pix, y, piy]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* psi = angle psi
	* ustar_ns = Nend - Nstar_ns
	* Nsteps = number of steps in the interpolation procedure to get a better precion (fixed)
	* method = method used to solve the system of differential equations: RK45, RK23, Radau, BDF, LSODA (fixed: 'RK45')
	* dense_output = False (True if we want a callable function in output)
	* events = None
	* rtol = relative accuracy, i.e. number of correct digits (fixed: 1.e-4)
 	* atol = absolute tollerance with local error given by atol+rtol*sol (fixed: 1.e-7)
 	* ist_eval = True when you want the solution to be computed at the point given in an array
 	Output: [NN_ustartot, dU, dUp, dV, dVp, T_, T_, dA, dS, R, S, T_, T_] x2
	"""

	x, pix, y, piy = xpixypiy
	Ntot = NN[-1]
	
	Up_tmp, Vp_tmp, Up2p_tmp, Vp2p_tmp, UpVpp_tmp, C1_tmp = Up2pVp2pUpVppC1(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi)
	Ve, VeU, VeV, epsU, epsV, etaUU_tmp, etaUV_tmp, etaVU_tmp, etaVV_tmp = SR_params_UV(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)
	
	Up = interpolate.InterpolatedUnivariateSpline(x=NN, y=Up_tmp, k=3)
	Vp = interpolate.InterpolatedUnivariateSpline(x=NN, y=Vp_tmp, k=3)
	Up2p = interpolate.InterpolatedUnivariateSpline(x=NN, y=Up2p_tmp, k=3)
	Vp2p = interpolate.InterpolatedUnivariateSpline(x=NN, y=Vp2p_tmp, k=3)
	UpVpp = interpolate.InterpolatedUnivariateSpline(x=NN, y=UpVpp_tmp, k=3)
	C1 = interpolate.InterpolatedUnivariateSpline(x=NN, y=C1_tmp, k=3)

	etaUU = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaUU_tmp, k=3)
	etaUV = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaUV_tmp, k=3)
	etaVU = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaVU_tmp, k=3)
	etaVV = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaVV_tmp, k=3)  


	def dUppdVpp(NNi, dUdUpdVdVp) :
		"""
		Equations of motion for dU and dV as given by eq.(5) in https://arxiv.org/pdf/astro-ph/0605679.pdf.
		Input:
		* NNi = integration variable: number of e-foldings
		* dUdUpdVdVp = [dU, dUp, dV, dVp] which are the variables we want to compute
		Output: [dUpp, dVpp]
		"""

		dU, dUp, dV, dVp = dUdUpdVdVp

		dC1 = C1(NNi)*(Up(NNi)*dU + Vp(NNi)*dV)

		dUpp = - C1(NNi)*dUp - 3.0*(etaUU(NNi)*dU + etaUV(NNi)*dV) + Up(NNi)*dC1 + Up2p(NNi)*dU + UpVpp(NNi)*dV
		dVpp = - C1(NNi)*dVp - 3.0*(etaVV(NNi)*dV + etaVU(NNi)*dU) + Vp(NNi)*dC1 + Vp2p(NNi)*dV + UpVpp(NNi)*dU

		return [dUp, dUpp, dVp, dVpp]


	NN0_ustartot = np.linspace(start=ustar_ns, stop=Ntot, num=Nsteps)

	Ud_tmp, Vd_tmp = UdVd(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)
	theta_tmp = np.arctan(Vd_tmp/Ud_tmp)
	theta_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=theta_tmp, k=3)
	theta_ustar = theta_interp(ustar_ns)
	costh_ustar = np.cos(theta_ustar) # equivalent to cphi above
	sinth_ustar = np.sin(theta_ustar) # equivalent to sphi above

	
	## find the initial values of dUp and dVp at ustar = Nend - Nstar such that the condition dUpp = dVpp = 0
	#--- dA* = 1.0, dS* = 0.0 ---#
	dA0_10 = 1.0
	dS0_10 = 0.0
	dU0_10 = costh_ustar*dA0_10 - sinth_ustar*dS0_10
	dV0_10 = sinth_ustar*dA0_10 + costh_ustar*dS0_10
	dUdUpdVdVp0_10 = [dU0_10, 0.0, dV0_10, 0.0]
	dUp_10_ustar, dUpp_10_ustar, dVp_10_ustar, dVpp_10_ustar = dUppdVpp(NNi=ustar_ns, dUdUpdVdVp=dUdUpdVdVp0_10)

	if C1(ustar_ns) != 0.0 :
		dUp0_10 = dUpp_10_ustar/C1(ustar_ns)  # this is because if we set dUp = dVp = 0 and we consider dUpp/C1 and dVpp/C1 we get exactly...
		dVp0_10 = dVpp_10_ustar/C1(ustar_ns)  # ...the same situation if we set dUpp = dVpp = 0 with dUp != 0 and dVp != 0 (look at the above equations!)
	else :
		print "C1(ustar_ns) = 0!"
		sys.exit(1)

	print "[dU0_10, dUp0_10, dV0_10, dVp0_10] = [%.8e, %.8e, %.8e, %.8e]" %(Decimal(dU0_10), Decimal(dUp0_10), Decimal(dV0_10), Decimal(dVp0_10))

	if ist_eval :
		sol_10 = integrate.solve_ivp(lambda NNi, dUdUpdVdVp: dUppdVpp(NNi, dUdUpdVdVp), [ustar_ns, Ntot], y0=[dU0_10, dUp0_10, dV0_10, dVp0_10], method=method, t_eval=NN0_ustartot, dense_output=dense_output, events=events, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	else :
		sol_10 = integrate.solve_ivp(lambda NNi, dUdUpdVdVp: dUppdVpp(NNi, dUdUpdVdVp), [ustar_ns, Ntot], y0=[dU0_10, dUp0_10, dV0_10, dVp0_10], method=method, t_eval=None, dense_output=dense_output, events=events, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	dU_10 = sol_10.y[0, :]
	dUp_10 = sol_10.y[1, :]
	dV_10 = sol_10.y[2, :]
	dVp_10 = sol_10.y[3, :]
	NN_ustartot_10 = sol_10.t # it should be equal to NN0_ustartot
	
	#--- dA* = 0.0, dS* = 1.0 ---#
	dA0_01 = 0.0
	dS0_01 = 1.0
	dU0_01 = costh_ustar*dA0_01 - sinth_ustar*dS0_01
	dV0_01 = sinth_ustar*dA0_01 + costh_ustar*dS0_01
	dUdUpdVdVp0_01 = [dU0_01, 0.0, dV0_01, 0.0]
	dUp_01_ustar, dUpp_01_ustar, dVp_01_ustar, dVpp_01_ustar = dUppdVpp(NNi=ustar_ns, dUdUpdVdVp=dUdUpdVdVp0_01)
	
	if C1(ustar_ns) != 0.0 :
		dUp0_01 = dUpp_01_ustar/C1(ustar_ns)
		dVp0_01 = dVpp_01_ustar/C1(ustar_ns)
	else :
		print "C1(ustar_ns) = 0!"
		sys.exit(1)

	print "[dU0_01, dUp0_01, dV0_01, dVp0_01] = [%.8e, %.8e, %.8e, %.8e]" %(Decimal(dU0_01), Decimal(dUp0_01), Decimal(dV0_01), Decimal(dVp0_01))

	if ist_eval :
		sol_01 = integrate.solve_ivp(lambda NNi, dUdUpdVdVp: dUppdVpp(NNi, dUdUpdVdVp), [ustar_ns, Ntot], y0=[dU0_01, dUp0_01, dV0_01, dVp0_01], method=method, t_eval=NN0_ustartot, dense_output=dense_output, events=events, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	else :
		sol_01 = integrate.solve_ivp(lambda NNi, dUdUpdVdVp: dUppdVpp(NNi, dUdUpdVdVp), [ustar_ns, Ntot], y0=[dU0_01, dUp0_01, dV0_01, dVp0_01], method=method, t_eval=None, dense_output=dense_output, events=events, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	dU_01 = sol_01.y[0, :]
	dUp_01 = sol_01.y[1, :]
	dV_01 = sol_01.y[2, :]
	dVp_01 = sol_01.y[3, :]
	NN_ustartot_01 = sol_01.t # it should be equal to NN0_ustartot

	
	if np.array_equal(NN_ustartot_10, NN_ustartot_01) :
		NN_ustartot = NN_ustartot_10
	else :
		print "'NN_ustartot_10' and 'NN_ustartot_01' are unequal!"
		sys.exit(1)
	
	xx = interpolate.InterpolatedUnivariateSpline(x=NN, y=x, k=3)
	pixpix = interpolate.InterpolatedUnivariateSpline(x=NN, y=pix, k=3)
	yy = interpolate.InterpolatedUnivariateSpline(x=NN, y=y, k=3)
	piypiy = interpolate.InterpolatedUnivariateSpline(x=NN, y=piy, k=3)
	psipsi = interpolate.InterpolatedUnivariateSpline(x=NN, y=psi, k=3)
	x_ustartot = xx(NN_ustartot)
	pix_ustartot = pixpix(NN_ustartot)
	y_ustartot = yy(NN_ustartot)
	piy_ustartot = piypiy(NN_ustartot)
	psi_ustartot = psipsi(NN_ustartot)
	xpixypiy_ustartot = [x_ustartot, pix_ustartot, y_ustartot, piy_ustartot]

	Ud, Vd = UdVd(xpixypiy=xpixypiy_ustartot, mllpxilpp=mllpxilpp, psi=psi_ustartot)
	theta = np.arctan(Vd/Ud)
	costh = np.cos(theta) # equivalent to cphi above
	sinth = np.sin(theta) # equivalent to sphi above
	Ad = np.sqrt(Ud**2 + Vd**2) 	# dot(sigma)
	H = Hubble(xpixypiy=xpixypiy_ustartot, mllpxilpp=mllpxilpp)

	
	#--- dA* = 1.0, dS* = 0.0 ---#
	uu_10 = H/Ad * dU_10
	vv_10 = H/Ad * dV_10
	T_UU = uu_10/uu_10[0]
	T_VU = vv_10/uu_10[0]
 
	dA_10 = costh*dU_10 + sinth*dV_10
	dS_10 = -sinth*dU_10 + costh*dV_10

	R_10 = H/Ad * dA_10
	S_10 = H/Ad * dS_10
	T_RR = R_10/R_10[0]
	T_SR = S_10/R_10[0]

	dAdS_10 = [NN_ustartot, dU_10, dUp_10, dV_10, dVp_10, T_UU, T_VU, dA_10, dS_10, R_10, S_10, T_RR, T_SR] 

	#--- dA* = 0.0, dS* = 1.0 ---#
	uu_01 = H/Ad * dU_01
	vv_01 = H/Ad * dV_01
	T_UV = uu_01/vv_01[0]
	T_VV = vv_01/vv_01[0]

	dA_01 = costh*dU_01 + sinth*dV_01
	dS_01 = -sinth*dU_01 + costh*dV_01
 
	R_01 = H/Ad * dA_01
	S_01 = H/Ad * dS_01
	T_RS = R_01/S_01[0]
	T_SS = S_01/S_01[0]

	dAdS_01 = [NN_ustartot, dU_01, dUp_01, dV_01, dVp_01, T_UV, T_VV, dA_01, dS_01, R_01, S_01, T_RS, T_SS]

	return [dAdS_10, dAdS_01]


def betaisoCosD(xpixypiypsi_ustar, mllpxilpp, P0s_ustar, transferAS_stop) :
	"""
	Compute the primordial isocurvature fraction beta_iso and the correlation fraction cos(Delta) at Nstop following the definitions present in https://arxiv.org/pdf/astro-ph/0605679.pdf and https://arxiv.org/pdf/1807.06211.pdf
	Input:
	* xpixypiypsi_ustar = [x, pix, y, piy] at ustar = Nend - Nstar
	* mllpxilpp = (m, l, lp, xi, lpp) parameter list
	* P0s_ustar = amplitude of the scalar power spectrum computed at ustar (https://arxiv.org/pdf/astro-ph/0605679.pdf)
	* transferAS_stop = [T_RR, T_RS, T_SR, T_SS] transfer fuctions computed at Nstop
	Output: [betaiso_stop, cosD_stop, cosD_stop1]
	"""

	x_ustar, pix_ustar, y_ustar, piy_ustar, psi_ustar = xpixypiypsi_ustar

	# compute epsU, epsV, etaUU, etaUV, etaVU, etaVV #
	Ve, VeU, VeV, epsU, epsV, etaUU, etaUV, etaVU, etaVV = SR_params_UV(xpixypiy=[x_ustar, pix_ustar, y_ustar, piy_ustar], mllpxilpp=mllpxilpp, psi=psi_ustar)
	
	# compute epsA, epsS, etaAA, etaAS, etaSS #
	Ud, Vd = UdVd(xpixypiy=[x_ustar, pix_ustar, y_ustar, piy_ustar], mllpxilpp=mllpxilpp, psi=psi_ustar)
	Ad = np.sqrt(Ud**2 + Vd**2) 	# dot(sigma)

	cphi = Ud/Ad
	sphi = Vd/Ad

	VeA = cphi*VeU + sphi*VeV
	VeS = -sphi*VeU + cphi*VeV

	epsA_ustar = VeA**2 / (2.0*Ve**2)
	epsS_ustar = VeS**2 / (2.0*Ve**2)
	etaAA_ustar = cphi**2*etaUU + sphi**2*etaVV + cphi*sphi*(etaUV + etaVU)
	etaAS_ustar = cphi*sphi*(etaVV - etaUU) + etaUV*(cphi**2 - sphi**2)
	etaSS_ustar = sphi**2*etaUU + cphi**2*etaVV - cphi*sphi*(etaUV + etaVU)

	# get the transfer functions AS #
	T_RR_stop, T_RS_stop, T_SR_stop, T_SS_stop = transferAS_stop 

	# compute the power spectra at Nstar #
	Cem = 2.0 - np.log(2.0) - 0.5772
	a1_ustar = -2.0*epsA_ustar + 2*Cem*(3*epsA_ustar - etaAA_ustar)
	a2_ustar = -2.0*Cem*etaAS_ustar
	a3_ustar = -2.0*epsA_ustar + 2*Cem*(epsA_ustar - etaSS_ustar)

	P_RR_ustar = P0s_ustar*(1.0 + a1_ustar)
	P_RS_ustar = P0s_ustar*a2_ustar
	P_SS_ustar = P0s_ustar*(1.0 + a3_ustar)

	# compute modified power spectra at Nstop #
	P_RR_stop = P_RR_ustar + 2.0*P_RS_ustar*T_RS_stop + P_SS_ustar*T_RS_stop**2
	P_RS_stop = P_RS_ustar*T_SS_stop + P_SS_ustar*T_SS_stop*T_RS_stop
	P_SS_stop = P_SS_ustar*T_SS_stop**2

	# compute beta_iso at Nstar since beta_iso is k-dependent and is a measure of the primordial isocurvature perturbations #
	#betaiso_ustar = P_SS_ustar / (P_RR_ustar + P_SS_ustar)
	betaiso_stop = P_SS_stop / (P_RR_stop + P_SS_stop) # as defined in https://arxiv.org/pdf/1807.06211.pdf
	
	# compute cos(Delta) at Nstop, which would remain almost constant till horizon re-entry #
	#cosD_ustar = P_RS_ustar / np.sqrt(P_RR_ustar*P_SS_ustar)
	cosD_stop = P_RS_stop / np.sqrt(P_RR_stop*P_SS_stop) # as defined in https://arxiv.org/pdf/1807.06211.pdf (correct definition!)
	cosD0_stop = T_RS_stop / np.sqrt(1.0 + T_RS_stop) # as defined in https://arxiv.org/pdf/astro-ph/0605679.pdf (zeroth order)

	#print "\nP0s_ustar = ", P0s_ustar
	#print "P_RR_ustar = ", P_RR_ustar
	#print "P_RS_ustar = ", P_RS_ustar
	#print "P_SS_ustar = ", P_SS_ustar
	#print "P_RR_stop = ", P_RR_stop
	#print "P_RS_stop = ", P_RS_stop
	#print "P_SS_stop = ", P_SS_stop
	#print "betaiso_stop = ", betaiso_stop
	#print "cosD_stop = ", cosD_stop
	#print "cosD0_stop = ", cosD0_stop

	return [betaiso_stop, cosD_stop, cosD0_stop]


def inflate(xpixypiypsi0, mllpxilpp, Nstart=0.0, Nsteps=1e3, eps=1e-10, method='RK45', rtol=1.e-4, atol=1.e-7, isrescale=True) :
	"""
	It allows to get all the parameters needed concerning inflation.
	Input:
	* xpixypiypsi0 = initial condition vector: [x0, pix00, y0, piy00, psi0]
	* mllpxilpp = list of constant parameters entering the potential: [m, l, lp, xi, lpp]
	* Nstart = 0.0
	* Nsteps = number of steps in the interpolation procedure to get a better precion (fixed)
	* eps = accuracy in Newton's method
	* method = method used to solve the system of differential equations: RK45, RK23, Radau, BDF, LSODA (fixed: 'RK45')
	* rtol = relative accuracy, i.e. number of correct digits (fixed: 1.e-4)
 	* atol = absolute tollerance with local error given by atol+rtol*sol (fixed: 1.e-7)
 	* rescale = True if we want to rescale Ve and rhoend
 	Output: [Nend, Nstar, As_star, ns_star, r_star, etaBlate, T_RR_stop, T_RS_stop, T_SR_stop, T_SS_stop]
	"""

	m, l, lp, xi, lpp = mllpxilpp
	
	x0, pix00, y0, piy00, psi0 = xpixypiypsi0
	pix0, piy0 = Solve_pi(xpixypiy0=[x0, pix00, y0, piy00], mllpxilpp=mllpxilpp, eps=eps)

	Fsol1, Nend = Solve_EoMs(NN0=[Nstart, 5000.0], xpixypiypsi0=[x0, pix0, y0, piy0, psi0], mllpxilpp=mllpxilpp, method=method, rtol=rtol, atol=atol, dense_output=True, ist_eval=False, stopEnd=True)

	xend, pixend, yend, piyend, psiend = Fsol1(Nend)
	Ntot = Nend + 15.0 # we stop after roughly 15 e-foldings after Nend
	
	Fsol2, Nfin = Solve_EoMs(NN0=[Nend, Ntot], xpixypiypsi0=[xend, pixend, yend, piyend, psiend], mllpxilpp=mllpxilpp, method=method, rtol=rtol, atol=atol, dense_output=True, ist_eval=False, stopEnd=False)
	NN1 = np.linspace(start=Nstart, stop=Nend, num=int(Nsteps*4/5))
	x1, pix1, y1, piy1, psi1 = Fsol1(NN1)
	NN2 = np.linspace(start=Nend, stop=Ntot, num=int(Nsteps/5+1))
	x2, pix2, y2, piy2, psi2 = Fsol2(NN2)
	
	NN = np.concatenate([NN1, np.delete(NN2, 0)])
	x = np.concatenate([x1, np.delete(x2, 0)])
	pix = np.concatenate([pix1, np.delete(pix2, 0)])
	y = np.concatenate([y1, np.delete(y2, 0)])
	piy = np.concatenate([piy1, np.delete(piy2, 0)])
	psi = np.concatenate([psi1, np.delete(psi2, 0)])

	# Compute Nstar #
	As, ns, r = nsrNEW(xpixypiy=[x, pix, y, piy], NN=NN, mllpxilpp=mllpxilpp, psi=psi, altDefVeS=False) # I have to compute it for all the values because As is used in computing Nstar!
	niend = NI(xpixypiy=[xend, pixend, yend, piyend], mllpxilpp=mllpxilpp)
	rhoend = m*niend
	
	kstar_ns = 0.05 # Mpc^-1, Planck pivot scale for As, ns
	kstar_r = 0.002 # Mpc^-1, Planck pivot scale for r
	Asobs = np.exp(3.044)*1.0e-10

	a0H0 = 67.8/3e5 #a0H0 = 67.66/2.99792458e5 # Mpc^-1
	rhoth = m*nrh(mllpxilpp=mllpxilpp)
	gth = 106.75 + 18.0
	cte = 67.0 - np.log(kstar_ns/a0H0) + (1.0/12)*np.log(rhoth/gth)	# wint = 0

	rescale = Asobs/As
	Ve = VE(x=x, y=y, mllpxilpp=mllpxilpp)
	if isrescale :
		Ve_res = np.abs(rescale*Ve)
		rhoend_res = np.abs(rescale*rhoend)
	else :
		Ve_res = np.abs(Ve)
		rhoend_res = np.abs(rhoend)

	fns = cte - (1.0/3)*np.log(rhoend_res) + 0.5*np.log(Ve_res) - NN
	fr = fns + np.log(kstar_ns/kstar_r)

	fns_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=fns)
	Nstar_ns, numIter_ns = secant(fns_interp, x0=60.0, x1=50.0, eps=eps)
	ustar_ns = Nend - Nstar_ns
	fr_interp = interpolate.InterpolatedUnivariateSpline(x=NN, y=fr)
	Nstar_r, numIter_r = secant(fr_interp, x0=60.0, x1=50.0, eps=eps)
	ustar_r = Nend - Nstar_r

	# Compute As, ns, r at Nstar #
	if Nstar_ns <= Nend :
		x_ustar_ns, pix_ustar_ns, y_ustar_ns, piy_ustar_ns, psi_ustar_ns = Fsol1(ustar_ns)
		Vstar_ns = VE(x=x_ustar_ns, y=y_ustar_ns, mllpxilpp=mllpxilpp)  
	else :
		x_ustar_ns, pix_ustar_ns, y_ustar_ns, piy_ustar_ns, psi_ustar_ns = Fsol2(ustar_ns) 
		Vstar_ns = VE(x=x_ustar_ns, y=y_ustar_ns, mllpxilpp=mllpxilpp)
	rescalestar_ns = np.exp(6.0*Nstar_ns - 6.0*cte + 2.0*np.log(rhoend) - 3.0*np.log(Vstar_ns))
	As_star_ns, ns_star_ns, r_star_ns = nsrNEW(xpixypiy=[x_ustar_ns, pix_ustar_ns, y_ustar_ns, piy_ustar_ns], NN=Nstar_ns, mllpxilpp=mllpxilpp, psi=psi_ustar_ns, altDefVeS=False)
	if Nstar_r <= Nend :
		x_ustar_r, pix_ustar_r, y_ustar_r, piy_ustar_r, psi_ustar_r = Fsol1(ustar_r)
	else :
		x_ustar_r, pix_ustar_r, y_ustar_r, piy_ustar_r, psi_ustar_r = Fsol2(ustar_r)
		
	As_star_r, ns_star_r, r_star_r = nsrNEW(xpixypiy=[x_ustar_r, pix_ustar_r, y_ustar_r, piy_ustar_r], NN=Nstar_r, mllpxilpp=mllpxilpp, psi=psi_ustar_r, altDefVeS=False)	

	Nstar = Nstar_ns
	As_star = As_star_ns
	ns_star = ns_star_ns
	r_star = r_star_r

	# Compute etaB #
	Nlate = Nend + 10.0 # compute etaB after 10 e-foldings after Nend
	xlate, pixlate, ylate, piylate, psilate = Fsol2(Nlate)
	etaBlate = etaB(xpixypiye=[xlate, pixlate, ylate, piylate], mllpxilpp=mllpxilpp, corr=True)

	# Compute T_RR, T_RS, T_SR, T_SS #
	Up_tmp, Vp_tmp, Up2p_tmp, Vp2p_tmp, UpVpp_tmp, C1_tmp = Up2pVp2pUpVppC1(xpixypiy=[x, pix, y, piy], NN=NN, mllpxilpp=mllpxilpp, psi=psi)
	Ve, VeU, VeV, epsU, epsV, etaUU_tmp, etaUV_tmp, etaVU_tmp, etaVV_tmp = SR_params_UV(xpixypiy=[x, pix, y, piy], mllpxilpp=mllpxilpp, psi=psi)
	
	Up = interpolate.InterpolatedUnivariateSpline(x=NN, y=Up_tmp, k=3)
	Vp = interpolate.InterpolatedUnivariateSpline(x=NN, y=Vp_tmp, k=3)
	Up2p = interpolate.InterpolatedUnivariateSpline(x=NN, y=Up2p_tmp, k=3)
	Vp2p = interpolate.InterpolatedUnivariateSpline(x=NN, y=Vp2p_tmp, k=3)
	UpVpp = interpolate.InterpolatedUnivariateSpline(x=NN, y=UpVpp_tmp, k=3)
	C1 = interpolate.InterpolatedUnivariateSpline(x=NN, y=C1_tmp, k=3)

	etaUU = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaUU_tmp, k=3)
	etaUV = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaUV_tmp, k=3)
	etaVU = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaVU_tmp, k=3)
	etaVV = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaVV_tmp, k=3)  


	def dUppdVpp(NNi, dUdUpdVdVp) :
		"""
		Equations of motion for dU and dV as given by eq.(5) in https://arxiv.org/pdf/astro-ph/0605679.pdf.
		Input:
		* NNi = integration variable: number of e-foldings
		* dUdUpdVdVp = [dU, dUp, dV, dVp] which are the variables we want to compute
		Output: [dUpp, dVpp]
		"""

		dU, dUp, dV, dVp = dUdUpdVdVp

		dC1 = C1(NNi)*(Up(NNi)*dU + Vp(NNi)*dV)

		dUpp = - C1(NNi)*dUp - 3.0*(etaUU(NNi)*dU + etaUV(NNi)*dV) + Up(NNi)*dC1 + Up2p(NNi)*dU + UpVpp(NNi)*dV
		dVpp = - C1(NNi)*dVp - 3.0*(etaVV(NNi)*dV + etaVU(NNi)*dU) + Vp(NNi)*dC1 + Vp2p(NNi)*dV + UpVpp(NNi)*dU

		return [dUp, dUpp, dVp, dVpp]


	Ud_ustar, Vd_ustar = UdVd(xpixypiy=[x_ustar_ns, pix_ustar_ns, y_ustar_ns, piy_ustar_ns], mllpxilpp=mllpxilpp, psi=psi_ustar_ns)
	theta_ustar = np.arctan(Vd_ustar/Ud_ustar)
	costh_ustar = np.cos(theta_ustar) # equivalent to cphi above
	sinth_ustar = np.sin(theta_ustar) # equivalent to sphi above
	Ad_ustar = np.sqrt(Ud_ustar**2 + Vd_ustar**2) 	# dot(sigma)
	H_ustar = Hubble(xpixypiy=[x_ustar_ns, pix_ustar_ns, y_ustar_ns, piy_ustar_ns], mllpxilpp=mllpxilpp)
	P0s_ustar = H_ustar**4 / (2*np.pi*Ad_ustar)**2

	Nstop = Nend - 5.0 # we consider the values of T_RR, T_RS, T_SR, T_SS at 5 e-folding before Nend
	#NN_ustarend = np.linspace(start=ustar_ns, stop=Nend, num=Nsteps) # we stop to Nend because there is no physical reason to go till Ntot
	
	## find the initial values of dUp and dVp at ustar = Nend - Nstar such that the condition dUpp = dVpp = 0
	#--- dA* = 1.0, dS* = 0.0 ---#
	dA0_10 = 1.0
	dS0_10 = 0.0
	dU0_10 = costh_ustar*dA0_10 - sinth_ustar*dS0_10
	dV0_10 = sinth_ustar*dA0_10 + costh_ustar*dS0_10
	dUp_10_ustar, dUpp_10_ustar, dVp_10_ustar, dVpp_10_ustar = dUppdVpp(NNi=ustar_ns, dUdUpdVdVp=[dU0_10, 0.0, dV0_10, 0.0])

	C1_ustar = C1(ustar_ns)
	if C1_ustar != 0.0 :
		dUp0_10 = dUpp_10_ustar/C1_ustar  # this is because if we set dUp = dVp = 0 and we consider dUpp/C1 and dVpp/C1 we get exactly...
		dVp0_10 = dVpp_10_ustar/C1_ustar  # ...the same situation if we set dUpp = dVpp = 0 with dUp != 0 and dVp != 0 (look at the above equations!)
	else :
		return [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]

	# we stop the intergration at Nend and not at Ntot
	sol_10 = integrate.solve_ivp(lambda NNi, dUdUpdVdVp: dUppdVpp(NNi, dUdUpdVdVp), [ustar_ns, Nend], y0=[dU0_10, dUp0_10, dV0_10, dVp0_10], method=method, t_eval=None, dense_output=True, events=None, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	dU_10_stop, dUp_10_stop, dV_10_stop, dVp_10_stop = sol_10.sol(Nstop)
	dU_10_ustar, dUp_10_ustar, dV_10_ustar, dVp_10_ustar = sol_10.sol(ustar_ns)
	
	#--- dA* = 0.0, dS* = 1.0 ---#
	dA0_01 = 0.0
	dS0_01 = 1.0
	dU0_01 = costh_ustar*dA0_01 - sinth_ustar*dS0_01
	dV0_01 = sinth_ustar*dA0_01 + costh_ustar*dS0_01
	dUp_01_ustar, dUpp_01_ustar, dVp_01_ustar, dVpp_01_ustar = dUppdVpp(NNi=ustar_ns, dUdUpdVdVp=[dU0_01, 0.0, dV0_01, 0.0])
	
	if C1_ustar != 0.0 :
		dUp0_01 = dUpp_01_ustar/C1_ustar
		dVp0_01 = dVpp_01_ustar/C1_ustar
	else :
		return [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]

	# we stop the intergration at Nend and not at Ntot
	sol_01 = integrate.solve_ivp(lambda NNi, dUdUpdVdVp: dUppdVpp(NNi, dUdUpdVdVp), [ustar_ns, Nend], y0=[dU0_01, dUp0_01, dV0_01, dVp0_01], method=method, t_eval=None, dense_output=True, events=None, rtol=rtol, atol=atol)#, rtol=1.e-5, atol=1.e-8)
	dU_01_stop, dUp_01_stop, dV_01_stop, dVp_01_stop = sol_01.sol(Nstop)
	dU_01_ustar, dUp_01_ustar, dV_01_ustar, dVp_01_ustar = sol_01.sol(ustar_ns)

	xstop, pixstop, ystop, piystop, psistop = Fsol1(Nstop)
	
	Ud_stop, Vd_stop = UdVd(xpixypiy=[xstop, pixstop, ystop, piystop], mllpxilpp=mllpxilpp, psi=psistop)
	theta_stop = np.arctan(Vd_stop/Ud_stop)
	costh_stop = np.cos(theta_stop) # equivalent to cphi above
	sinth_stop = np.sin(theta_stop) # equivalent to sphi above
	Ad_stop = np.sqrt(Ud_stop**2 + Vd_stop**2) 	# dot(sigma)
	H_stop = Hubble(xpixypiy=[xstop, pixstop, ystop, piystop], mllpxilpp=mllpxilpp)

	
	#--- dA* = 1.0, dS* = 0.0 ---#
	dA_10_stop = costh_stop*dU_10_stop + sinth_stop*dV_10_stop
	dS_10_stop = -sinth_stop*dU_10_stop + costh_stop*dV_10_stop
	dA_10_ustar = costh_ustar*dU_10_ustar + sinth_ustar*dV_10_ustar
	dS_10_ustar = -sinth_ustar*dU_10_ustar + costh_ustar*dV_10_ustar

	R_10_stop = H_stop/Ad_stop * dA_10_stop
	S_10_stop = H_stop/Ad_stop * dS_10_stop
	R_10_ustar = H_ustar/Ad_ustar * dA_10_ustar
	S_10_ustar = H_ustar/Ad_ustar * dS_10_ustar
	T_RR_stop = R_10_stop/R_10_ustar # if Nstop = Nstar --> T_RR_ustar = 1
	T_SR_stop = S_10_stop/R_10_ustar # if Nstop = Nstar --> T_SR_ustar = 0

	
	#--- dA* = 0.0, dS* = 1.0 ---#
	dA_01_stop = costh_stop*dU_01_stop + sinth_stop*dV_01_stop
	dS_01_stop = -sinth_stop*dU_01_stop + costh_stop*dV_01_stop
	dA_01_ustar = costh_ustar*dU_01_ustar + sinth_ustar*dV_01_ustar
	dS_01_ustar = -sinth_ustar*dU_01_ustar + costh_ustar*dV_01_ustar
 
	R_01_stop = H_stop/Ad_stop * dA_01_stop
	S_01_stop = H_stop/Ad_stop * dS_01_stop
	R_01_ustar = H_ustar/Ad_ustar * dA_01_ustar
	S_01_ustar = H_ustar/Ad_ustar * dS_01_ustar
	T_RS_stop = R_01_stop/S_01_ustar # if Nstop = Nstar --> T_RS_ustar = 0
	T_SS_stop = S_01_stop/S_01_ustar # if Nstop = Nstar --> T_SS_ustar = 1

	# It is useless to adjust the definition of As at Nstar (i.e. ustar_ns) because T_RS_ustar = 0 and so as As_corr #
	betaiso_stop, cosD_stop, cosD0_stop = betaisoCosD(xpixypiypsi_ustar=[x_ustar_ns, pix_ustar_ns, y_ustar_ns, piy_ustar_ns, psi_ustar_ns], mllpxilpp=mllpxilpp, P0s_ustar=P0s_ustar, transferAS_stop=[T_RR_stop, T_RS_stop, T_SR_stop, T_SS_stop])
	
	return [Nend, Nstar, As_star, ns_star, r_star, etaBlate, T_RR_stop, T_RS_stop, T_SR_stop, T_SS_stop, betaiso_stop, cosD_stop]


def chi2(As, ns, r, etaB, betaiso, cosD, isconvfact=True) :
	"""
	Chi2 function for the observed vs found values of As, ns, r and etaB.
	Input:
	* As = found value of As at Nstar
	* ns = found value of ns at Nstar
	* r = found value of r at Nstar
	* etaB = found value of etaB at Nstop
	* betaiso = found value of betaiso at Nstar
	* cosD = found value of cos(Delta) at Nstop
	* isconvfact = True if we want to use the conversion factor between CDI and BDI in Planck data (is it correct?)
	Output: chi2f
	"""

	Asobs = np.exp(3.044)*1.0e-10 # https://arxiv.org/pdf/1807.06211.pdf
	dAsobs = np.exp(0.014)*1.0e-10 # https://arxiv.org/pdf/1807.06211.pdf
	nsobs = 0.9649 # https://arxiv.org/pdf/1807.06211.pdf
	dnsobs = 0.0042 # https://arxiv.org/pdf/1807.06211.pdf
	etaBobs = 8.64e-11 # https://arxiv.org/pdf/1807.06209.pdf, http://background.uchicago.edu/~whu/Courses/Ast321_17/Baryogen_Lecture.pdf 
	detaBobs = 5.79e-13 # https://arxiv.org/pdf/1807.06209.pdf, http://background.uchicago.edu/~whu/Courses/Ast321_17/Baryogen_Lecture.pdf 
	robs = 0.056 # https://arxiv.org/pdf/1807.06209.pdf
	drobs = 0.05 # ???

	# beta_iso CDI bounds are the same as those for BDI
	betaisobs = 0.1375 # average over different data sets with 3 isocurvature parameters for bold CDI Planck TT,TE,EE+lowE+lensing (https://arxiv.org/pdf/1807.06211.pdf, Table 14)
	dbetaisobs = 0.1473 # standard deviation over different data sets with 3 isocurvature parameters for bold CDI Planck TT,TE,EE+lowE+lensing (https://arxiv.org/pdf/1807.06211.pdf, Table 14)
	# cos(Delta) CDI bounds are not the same as those for BDI by I should multiply the bounds by 5.33/28.4
	if isconvfact :
		cosDobs = (5.33/28.4)*0.0225
		dcosDobs = (5.33/28.4)*0.1533
	else :
		cosDobs = 0.0225 # average over different data sets with 3 isocurvature parameters for bold CDI Planck TT,TE,EE+lowE+lensing (https://arxiv.org/pdf/1807.06211.pdf, Table 14)
		dcosDobs = 0.1533 # standard deviation over different data sets with 3 isocurvature parameters for bold CDI Planck TT,TE,EE+lowE+lensing (https://arxiv.org/pdf/1807.06211.pdf, Table 14)

	
	chi2f = (As - Asobs)**2 / dAsobs**2 + (ns - nsobs)**2 / dnsobs**2 + (etaB - etaBobs)**2 / detaBobs**2 + (r - robs)**2 / drobs**2 + (betaiso - betaisobs)**2 / dbetaisobs**2 + (cosD - cosDobs)**2 / dcosDobs**2

	if r >= robs or betaiso >= betaisobs or cosD >= cosDobs :
		chi2f += 15.0 # I disfavour these models

	return chi2f


def genpar(oldparam, parchange=0.05) :
	"""
	Generate new parameters starting from the old ones.
	Input: 
	* oldparam = old parameters
	* parchange = parameter change fraction (optional: 0.05)
	Output: newparam
	"""

	newparam = []
	
	for j in range(0, len(oldparam)) :
		rand = random.uniform(0, 1)
		randp = 2.0*(rand - 0.5)
		newparam_j = oldparam[j]*(1 + randp*parchange)
		newparam.append(newparam_j)	

	return tuple(newparam)


def MCMC(fname, param, matry=6e7, parchange=0.05, Nstart=0.0, Nsteps=1e3, eps=1e-10, method='RK45', rtol=1.e-4, atol=1.e-7, isrescale=False, isconvfact=True) :
	"""
	Markov Chain Monte Carlo used to scan over the parameter space.
	Input:
	* fname = string name of the file fchain
	* param = initial condition parameters: (m, l, lp, xi, lpp, x0, y0)
	* matry = dimension of the initial chain
	* parchange = parameter change fraction
	* Nstart = 0.0
	* Nsteps = number of steps in the interpolation procedure to get a better precion (fixed)
	* eps = accuracy in Newton's method
	* method = method used to solve the system of differential equations: RK45, RK23, Radau, BDF, LSODA (fixed: 'RK45')
	* rtol = relative accuracy, i.e. number of correct digits (fixed: 1.e-4)
 	* atol = absolute tollerance with local error given by atol+rtol*sol (fixed: 1.e-7)
 	* isrescale = True if we want to rescale the potential VE and rho_end parameters such that As = Abobs in order to find Nstar
	* isconvfact = True if we want to use the conversion factor between CDI and BDI in Planck data (is it correct?)
	Output: None
	"""

	str_chain = fname + ".dat"
	fchain = open(str_chain, "w+")
	fchain.write("m\t\t\tl\t\t\tlp\t\t\txi\t\t\tlpp\t\t\tx0\t\t\ty0\t\t\tNend\t\t\tNstar\t\t\tAs\t\t\tns\t\t\tr\t\t\tetaB\t\t\tbeta_iso\t\t\tcosD\t\t\tchi2\n")

	m, l, lp, xi, lpp, x0, y0 = param
	pix00 = 0.0
	piy00 = 0.0
	psi0 = 0.0

	data = inflate(xpixypiypsi0=[x0, pix00, y0, piy00, psi0], mllpxilpp=(m, l, lp, xi, lpp), Nstart=Nstart, Nsteps=Nsteps, eps=eps, method=method, rtol=rtol, atol=atol, isrescale=isrescale)
	Nend, Nstar, As_star, ns_star, r_star, etaBlate, T_RR_stop, T_RS_stop, T_SR_stop, T_SS_stop, betaiso_stop, cosD_stop = data

	if data == [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1] :
		print "C1_ustar == 0.0 in the initial model!"
		sys.exit(1)
	
	chi2f = chi2(As=As_star, ns=ns_star, r=r_star, etaB=etaBlate, betaiso=betaiso_stop, cosD=cosD_stop, isconvfact=isconvfact)
		
	fchain.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(m, l, lp, xi, lpp, x0, y0, Nend, Nstar, As_star, ns_star, r_star, etaBlate, betaiso_stop, cosD_stop, chi2f))
	param_tmp = list(param)

	print "IC"

	ftrial = open('totChain.dat', "w+")
	ftrial.write("m\t\t\tl\t\t\tlp\t\t\txi\t\t\tlpp\t\t\tx0\t\t\ty0\t\t\tNend\t\t\tNstar\t\t\tAs\t\t\tns\t\t\tr\t\t\tetaB\t\t\tbeta_iso\t\t\tcosD\t\t\tchi2\n")
	ftrial.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(m, l, lp, xi, lpp, x0, y0, Nend, Nstar, As_star, ns_star, r_star, etaBlate, betaiso_stop, cosD_stop, chi2f))

	Nchain = 0
	for i in range(0, int(matry)) :
		m_i, l_i, lp_i, xi_i, lpp_i, x0_i, y0_i = genpar(oldparam=param_tmp, parchange=parchange)
		
		data_i = inflate(xpixypiypsi0=[x0_i, pix00, y0_i, piy00, psi0], mllpxilpp=(m_i, l_i, lp_i, xi_i, lpp_i), Nstart=Nstart, Nsteps=Nsteps, eps=eps, method=method, rtol=rtol, atol=atol, isrescale=isrescale)
		Nend_i, Nstar_i, As_star_i, ns_star_i, r_star_i, etaBlate_i, T_RR_stop_i, T_RS_stop_i, T_SR_stop_i, T_SS_stop_i, betaiso_stop_i, cosD_stop_i = data_i
		
		if data_i == [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1] :
			continue
		
		chi2f_i = chi2(As=As_star_i, ns=ns_star_i, r=r_star_i, etaB=etaBlate_i, betaiso=betaiso_stop_i, cosD=cosD_stop_i, isconvfact=isconvfact)
			
		ftrial.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(m_i, l_i, lp_i, xi_i, lpp_i, x0_i, y0_i, Nend_i, Nstar_i, As_star_i, ns_star_i, r_star_i, etaBlate_i, betaiso_stop_i, cosD_stop_i, chi2f_i))

		randn = random.uniform(0, 1)
		if chi2f/chi2f_i > randn :
			if chi2f_i <= 14.067 : # if the chi2 is within 95% C.L. (http://uregina.ca/~gingrich/appchi.pdf)
				fchain.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(m_i, l_i, lp_i, xi_i, lpp_i, x0_i, y0_i, Nend_i, Nstar_i, As_star_i, ns_star_i, r_star_i, etaBlate_i, betaiso_stop_i, cosD_stop_i, chi2f_i))
				param_tmp = [m_i, l_i, lp_i, xi_i, lpp_i, x0_i, y0_i]
				chi2f = chi2f_i
				Nchain += 1

		print i

	return

