import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0, 'font.size': 16})

from matplotlib import rc
rc('text', usetex=True)
rc('font', size=17)
rc('legend', fontsize=15)

from scipy.signal import savgol_filter
from scipy import interpolate

import matplotlib.offsetbox as offsetbox
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import NullFormatter, FixedLocator
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

from mpl_toolkits.mplot3d import Axes3D

import myTools


def InterpAndSmooth(x, y, intkind='linear', Nsteps=1e3) :
	"""
	Data interpolation and smoothing.
	Input:
	* x = list or array corresponding to the x-axis (linear NOT log!)
	* y = list or array corresponding to the y-axis
	* intkind = type of interpolation (fixed: linear)
	* Nsteps = number or interpolation points
	Output: [xx, yy]
	"""

	xx = np.linspace(x.min(), x.max(), Nsteps)
	interp = interpolate.interp1d(x, y, kind=intkind)
	yy = savgol_filter(interp(xx), 101, 1)

	return [xx, yy]


def rnsPlot(outDirName, title, fname1='nsr', fname2=None, xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, nfiles=None) :
	"""
	Plot of r vs n_s. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* fname2 = name of the file containing the data for comparison (optional)
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	colorList = ["brown", "black", "red", "green", "orange", "magenta", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	nsrList = [] # for k* = 0.002/Mpc
	
	ns_ns_starList = [] # for k* = 0.05/Mpc
	r_ns_starList = []
	ns_r_starList = [] # for k* = 0.002/Mpc
	r_r_starList = []

	PlanckData_TTEEBK15BAO_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15BAO_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 1)

	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		nsr_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
		nsrList.append(nsr_i)
		if fname2 != None :
			nsrNew_i = np.loadtxt(fname=outDirName + '/' + fname2 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
			nsrNewList.append(nsrNew_i)

		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "ns_star_ns = " in line:
				value = line[13:-1]
				ns_star_i = float(value)
				ns_ns_starList.append(ns_star_i)
			if "r_star_ns = " in line:
				value = line[12:-1]
				r_star_i = float(value)
				r_ns_starList.append(r_star_i)
			if "ns_star_r = " in line:
				value = line[12:-1]
				ns_star_i = float(value)
				ns_r_starList.append(ns_star_i)
			if "r_star_r = " in line:
				value = line[11:-1]
				r_star_i = float(value)
				r_r_starList.append(r_star_i)

	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1

	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i])
		text1 = r"$\phi_0$ = " + str(x0List[i]) + r" + $i$ " + str(y0List[i])

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$n_s$', fontsize=18)
	ax.set_ylabel(r'$r$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$\xi =$ " + str(xiList[i])
		nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(nsns, r_sg, color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		#ax.plot(ns_ns_starList[i], r_ns_starList[i], color = colorList[i], marker = 's', markersize=12)
		ax.plot(ns_r_starList[i], r_r_starList[i], color = colorList[i], marker = '*', markersize=10)
		if i == 2 :
			ax.plot(9.61792031e-01, 1.33956485e-02, marker = '*', markersize=12, markerfacecolor="white", markeredgecolor="red")
			ax.plot(9.69247809e-01, 1.13397142e-02, marker = '*', markersize=12, markerfacecolor="white", markeredgecolor="red")
			ax.text(x=0.9617, y=0.0138, s=r"$\lambda'' = 1$e$-08$", color='red', fontsize=15)
			ax.text(x=0.965, y=0.01245, s=r"$\lambda'' = 8.67$e$-05$", color='red', fontsize=15)
			ax.text(x=0.9691, y=0.0117, s=r"$\lambda'' = 0.05$", color='red', fontsize=15)

	if fname2 != None :
		for i in range(0, num) :
			ax.plot(nsrNewList[i][:,2], nsrNewList[i][:,3], color = colorList[i], linestyle = '--', marker = '')	

	ax.plot(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color ='blue', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color='blue', alpha=0.1)
	ax.plot(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color ='cyan', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color='cyan', alpha=0.15)

#	ax.plot(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color ='red', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color='red', alpha=0.2)
#	ax.plot(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color ='orange', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color='orange', alpha=0.25)

#	ax.plot(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color ='black', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color='black', alpha=0.2)
#	ax.plot(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color ='gray', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color='gray', alpha=0.25)

	#handles, labels = plt.gca().get_legend_handles_labels() # get existing handles and labels
	#star_ns = Line2D([], [], marker='s', linestyle='', color='black', label=r'$k_{\star} =$ 0.05 Mpc$^{-1}$', markersize=12)
	#star_r = Line2D([], [], marker='*', linestyle='', color='black', label=r'$k_{\star} =$ 0.002 Mpc$^{-1}$', markersize=10)
	#handles.append(star_ns)
	#handles.append(star_r)
	#labels.append(r'$k_{\star} =$ 0.05 Mpc$^{-1}$')
	#labels.append(r'$k_{\star} =$ 0.002 Mpc$^{-1}$')
	#ax.legend(handles, labels, loc='lower left')
	ax.legend(loc='lower left')
	ax.text(x=0.965, y=0.0152, s="Planck", color='blue', fontsize=20)
	ax.text(x=0.956, y=0.0122, s=r"$N_{\star} = 50$", color='orange', fontsize=20)
	ax.text(x=0.9665, y=0.008, s=r"$N_{\star} = 60$", color='orange', fontsize=20)
	ob = offsetbox.AnchoredText(s=text, loc='upper left', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)
	ob1 = offsetbox.AnchoredText(s=text1, loc='lower right', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob1.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob1)
	ax.grid(False)
	fig.show()

	return


def rnsPlot1(outDirName, title, fname1='nsr', fname2=None, xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, nfiles=None) :
	"""
	Plot of r vs n_s. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* fname2 = name of the file containing the data for comparison (optional)
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	nsrList = [] # for k* = 0.002/Mpc
	
	ns_ns_starList = [] # for k* = 0.05/Mpc
	r_ns_starList = []
	ns_r_starList = [] # for k* = 0.002/Mpc
	r_r_starList = []

	PlanckData_TTEEBK15BAO_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15BAO_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 1)

	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		nsr_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
		nsrList.append(nsr_i)
		if fname2 != None :
			nsrNew_i = np.loadtxt(fname=outDirName + '/' + fname2 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
			nsrNewList.append(nsrNew_i)

		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "ns_star_ns = " in line:
				value = line[13:-1]
				ns_star_i = float(value)
				ns_ns_starList.append(ns_star_i)
			if "r_star_ns = " in line:
				value = line[12:-1]
				r_star_i = float(value)
				r_ns_starList.append(r_star_i)
			if "ns_star_r = " in line:
				value = line[12:-1]
				ns_star_i = float(value)
				ns_r_starList.append(ns_star_i)
			if "r_star_r = " in line:
				value = line[11:-1]
				r_star_i = float(value)
				r_r_starList.append(r_star_i)

	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1
	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i])
		text1 = r"$\phi_0$ = " + str(x0List[i]) + r" + $i$ " + str(y0List[i])

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$n_s$', fontsize=18)
	ax.set_ylabel(r'$r$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$\xi =$ " + str(xiList[i])
		nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(nsns, r_sg, color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		#ax.plot(ns_ns_starList[i], r_ns_starList[i], color = colorList[i], marker = 's', markersize=12)
		ax.plot(ns_r_starList[i], r_r_starList[i], color = colorList[i], marker = '*', markersize=10)

	if fname2 != None :
		for i in range(0, num) :
			ax.plot(nsrNewList[i][:,2], nsrNewList[i][:,3], color = colorList[i], linestyle = '--', marker = '')	

	ax.plot(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color ='blue', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color='blue', alpha=0.1)
	ax.plot(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color ='cyan', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color='cyan', alpha=0.15)

#	ax.plot(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color ='red', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color='red', alpha=0.2)
#	ax.plot(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color ='orange', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color='orange', alpha=0.25)

#	ax.plot(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color ='black', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color='black', alpha=0.2)
#	ax.plot(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color ='gray', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color='gray', alpha=0.25)

	#handles, labels = plt.gca().get_legend_handles_labels() # get existing handles and labels
	#star_ns = Line2D([], [], marker='s', linestyle='', color='black', label=r'$k_{\star} =$ 0.05 Mpc$^{-1}$', markersize=12)
	#star_r = Line2D([], [], marker='*', linestyle='', color='black', label=r'$k_{\star} =$ 0.002 Mpc$^{-1}$', markersize=10)
	#handles.append(star_ns)
	#handles.append(star_r)
	#labels.append(r'$k_{\star} =$ 0.05 Mpc$^{-1}$')
	#labels.append(r'$k_{\star} =$ 0.002 Mpc$^{-1}$')
	#ax.legend(handles, labels, loc='upper right')
	ax.legend(loc='upper right')
	ax.text(x=0.967, y=0.04, s="Planck", color='blue', fontsize=20)
	ax.text(x=0.9583, y=0.0428, s=r"$N_{\star} = 50$", color='brown', fontsize=20)
	ax.text(x=0.9688, y=0.030, s=r"$N_{\star} = 60$", color='brown', fontsize=20)
	ob = offsetbox.AnchoredText(s=text, loc='upper center', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)
	ob1 = offsetbox.AnchoredText(s=text1, loc='lower right', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob1.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob1)
	ax.grid(False)
	fig.show()

	return


def rnsPlotX0(outDirName, title, fname1='nsr', fname2=None, xlim=None, ylim=None, PlanckCoor=(0.965,0.0152), Nstar50=None, Nstar60=None, LogXaxis=False, LogYaxis=False, nfiles=None) :
	"""
	Plot of r vs n_s. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* fname2 = name of the file containing the data for comparison (optional)
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	##colorList = ["brown", "black", "red", "green", "yellow", "magenta", "blue"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red", "pink", "yellow", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	nsrList = [] # for k* = 0.002/Mpc
	
	ns_ns_starList = [] # for k* = 0.05/Mpc
	r_ns_starList = []
	ns_r_starList = [] # for k* = 0.002/Mpc
	r_r_starList = []

	PlanckData_TTEEBK15BAO_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15BAO_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 1)

	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		nsr_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
		nsrList.append(nsr_i)
		if fname2 != None :
			nsrNew_i = np.loadtxt(fname=outDirName + '/' + fname2 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
			nsrNewList.append(nsrNew_i)

		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "ns_star_ns = " in line:
				value = line[13:-1]
				ns_star_i = float(value)
				ns_ns_starList.append(ns_star_i)
			if "r_star_ns = " in line:
				value = line[12:-1]
				r_star_i = float(value)
				r_ns_starList.append(r_star_i)
			if "ns_star_r = " in line:
				value = line[12:-1]
				ns_star_i = float(value)
				ns_r_starList.append(ns_star_i)
			if "r_star_r = " in line:
				value = line[11:-1]
				r_star_i = float(value)
				r_r_starList.append(r_star_i)

	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1

	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\xi$ = " + str(xiList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i]) + "\n" + r"$Y_0$ = " + str(y0List[i])
		

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$n_s$', fontsize=18)
	ax.set_ylabel(r'$r$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$X_0 =$ " + str(x0List[i])
		nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(nsns, r_sg, color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		#ax.plot(ns_ns_starList[i], r_ns_starList[i], color = colorList[i], marker = 's', markersize=12)
		ax.plot(ns_r_starList[i], r_r_starList[i], color = colorList[i], marker = '*', markersize=10)

	if fname2 != None :
		for i in range(0, num) :
			ax.plot(nsrNewList[i][:,2], nsrNewList[i][:,3], color = colorList[i], linestyle = '--', marker = '')	

	ax.plot(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color ='blue', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color='blue', alpha=0.1)
	ax.plot(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color ='cyan', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color='cyan', alpha=0.15)

#	ax.plot(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color ='red', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color='red', alpha=0.2)
#	ax.plot(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color ='orange', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color='orange', alpha=0.25)

#	ax.plot(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color ='black', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color='black', alpha=0.2)
#	ax.plot(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color ='gray', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color='gray', alpha=0.25)

	#handles, labels = plt.gca().get_legend_handles_labels() # get existing handles and labels
	#star_ns = Line2D([], [], marker='s', linestyle='', color='black', label=r'$k_{\star} =$ 0.05 Mpc$^{-1}$', markersize=12)
	#star_r = Line2D([], [], marker='*', linestyle='', color='black', label=r'$k_{\star} =$ 0.002 Mpc$^{-1}$', markersize=10)
	#handles.append(star_ns)
	#handles.append(star_r)
	#labels.append(r'$k_{\star} =$ 0.05 Mpc$^{-1}$')
	#labels.append(r'$k_{\star} =$ 0.002 Mpc$^{-1}$')
	#ax.legend(handles, labels, loc='lower left')
	ax.legend(loc='lower left')
	xPlanck, yPlanck = PlanckCoor
	ax.text(x=xPlanck, y=yPlanck, s="Planck", color='blue', fontsize=20)
	if Nstar50 != None :
		xNstar50, yNstar50 = Nstar50
		ax.text(x=xNstar50, y=yNstar50, s=r"$N_{\star} = 50$", color='red', fontsize=20)
	if Nstar60 != None :
		xNstar60, yNstar60 = Nstar60
		ax.text(x=xNstar60, y=yNstar60, s=r"$N_{\star} = 60$", color='red', fontsize=20)
	ob = offsetbox.AnchoredText(s=text, loc='upper left', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)
	ax.grid(False)
	fig.show()

	return


def rnsPlot1X0(outDirName, title, fname1='nsr', fname2=None, xlim=None, ylim=None, PlanckCoor=(0.966,0.009), Nstar50=None, Nstar60=None, LogXaxis=False, textModel=None, LogYaxis=False, nfiles=None) :
	"""
	Plot of r vs n_s. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* fname2 = name of the file containing the data for comparison (optional)
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red", "yellow", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	nsrList = [] # for k* = 0.002/Mpc
	
	ns_ns_starList = [] # for k* = 0.05/Mpc
	r_ns_starList = []
	ns_r_starList = [] # for k* = 0.002/Mpc
	r_r_starList = []

	PlanckData_TTEEBK15BAO_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15BAO_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 1)

	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		nsr_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
		nsrList.append(nsr_i)
		if fname2 != None :
			nsrNew_i = np.loadtxt(fname=outDirName + '/' + fname2 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
			nsrNewList.append(nsrNew_i)

		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "ns_star_ns = " in line:
				value = line[13:-1]
				ns_star_i = float(value)
				ns_ns_starList.append(ns_star_i)
			if "r_star_ns = " in line:
				value = line[12:-1]
				r_star_i = float(value)
				r_ns_starList.append(r_star_i)
			if "ns_star_r = " in line:
				value = line[12:-1]
				ns_star_i = float(value)
				ns_r_starList.append(ns_star_i)
			if "r_star_r = " in line:
				value = line[11:-1]
				r_star_i = float(value)
				r_r_starList.append(r_star_i)

	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1
	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\xi$ = " + str(xiList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i]) + "\n" + r"$Y_0$ = " + str(y0List[i])
		

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$n_s$', fontsize=18)
	ax.set_ylabel(r'$r$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$X_0 =$ " + str(x0List[i])
		nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(nsns, r_sg, color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		#ax.plot(ns_ns_starList[i], r_ns_starList[i], color = colorList[i], marker = 's', markersize=12)
		ax.plot(ns_r_starList[i], r_r_starList[i], color = colorList[i], marker = '*', markersize=10)

	if fname2 != None :
		for i in range(0, num) :
			ax.plot(nsrNewList[i][:,2], nsrNewList[i][:,3], color = colorList[i], linestyle = '--', marker = '')	

	ax.plot(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color ='blue', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color='blue', alpha=0.1)
	ax.plot(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color ='cyan', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color='cyan', alpha=0.15)

#	ax.plot(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color ='red', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color='red', alpha=0.2)
#	ax.plot(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color ='orange', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color='orange', alpha=0.25)

#	ax.plot(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color ='black', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color='black', alpha=0.2)
#	ax.plot(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color ='gray', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color='gray', alpha=0.25)

	#handles, labels = plt.gca().get_legend_handles_labels() # get existing handles and labels
	#star_ns = Line2D([], [], marker='s', linestyle='', color='black', label=r'$k_{\star} =$ 0.05 Mpc$^{-1}$', markersize=12)
	#star_r = Line2D([], [], marker='*', linestyle='', color='black', label=r'$k_{\star} =$ 0.002 Mpc$^{-1}$', markersize=10)
	#handles.append(star_ns)
	#handles.append(star_r)
	#labels.append(r'$k_{\star} =$ 0.05 Mpc$^{-1}$')
	#labels.append(r'$k_{\star} =$ 0.002 Mpc$^{-1}$')
	#ax.legend(handles, labels, loc='upper right')
	ax.legend(loc='upper right')
	#ax.text(x=0.967, y=0.018, s="Planck", color='blue', fontsize=20)
	xPlanck, yPlanck = PlanckCoor
	ax.text(x=xPlanck, y=yPlanck, s="Planck", color='blue', fontsize=20)
	if Nstar50 != None :
		xNstar50, yNstar50 = Nstar50
		ax.text(x=xNstar50, y=yNstar50, s=r"$N_{\star} = 50$", color='red', fontsize=20)
	if Nstar60 != None :
		xNstar60, yNstar60 = Nstar60
		ax.text(x=xNstar60, y=yNstar60, s=r"$N_{\star} = 60$", color='red', fontsize=20)
	ob = offsetbox.AnchoredText(s=text, loc='lower left', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)

	if textModel != None :
		ob1 = offsetbox.AnchoredText(s=textModel, loc='lower right', pad=0.3, borderpad=0.85, frameon=False, prop=dict(color='black', size=15))
		ax.add_artist(ob1)		

	ax.grid(False)
	fig.show()

	return


def rnsPlotY0(outDirName, title, fname1='nsr', fname2=None, xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, PlanckCoor=(0.965,0.0152), Nstar50=None, Nstar60=None, nfiles=None) :
	"""
	Plot of r vs n_s. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* fname2 = name of the file containing the data for comparison (optional)
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	##colorList = ["brown", "black", "red", "green", "yellow", "magenta", "blue", "pink"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red", "pink", "yellow", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	nsrList = [] # for k* = 0.002/Mpc
	
	ns_ns_starList = [] # for k* = 0.05/Mpc
	r_ns_starList = []
	ns_r_starList = [] # for k* = 0.002/Mpc
	r_r_starList = []

	PlanckData_TTEEBK15BAO_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15BAO_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 1)

	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		nsr_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
		nsrList.append(nsr_i)
		if fname2 != None :
			nsrNew_i = np.loadtxt(fname=outDirName + '/' + fname2 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
			nsrNewList.append(nsrNew_i)

		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "ns_star_ns = " in line:
				value = line[13:-1]
				ns_star_i = float(value)
				ns_ns_starList.append(ns_star_i)
			if "r_star_ns = " in line:
				value = line[12:-1]
				r_star_i = float(value)
				r_ns_starList.append(r_star_i)
			if "ns_star_r = " in line:
				value = line[12:-1]
				ns_star_i = float(value)
				ns_r_starList.append(ns_star_i)
			if "r_star_r = " in line:
				value = line[11:-1]
				r_star_i = float(value)
				r_r_starList.append(r_star_i)

	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1

	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\xi$ = " + str(xiList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i]) + "\n" + r"$X_0$ = " + str(x0List[i])
		

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$n_s$', fontsize=18)
	ax.set_ylabel(r'$r$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$Y_0 =$ " + str(y0List[i])
		nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(nsns, r_sg, color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		#ax.plot(ns_ns_starList[i], r_ns_starList[i], color = colorList[i], marker = 's', markersize=12)
		ax.plot(ns_r_starList[i], r_r_starList[i], color = colorList[i], marker = '*', markersize=10)

	if fname2 != None :
		for i in range(0, num) :
			ax.plot(nsrNewList[i][:,2], nsrNewList[i][:,3], color = colorList[i], linestyle = '--', marker = '')	

	ax.plot(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color ='blue', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color='blue', alpha=0.1)
	ax.plot(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color ='cyan', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color='cyan', alpha=0.15)

#	ax.plot(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color ='red', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color='red', alpha=0.2)
#	ax.plot(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color ='orange', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color='orange', alpha=0.25)

#	ax.plot(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color ='black', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color='black', alpha=0.2)
#	ax.plot(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color ='gray', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color='gray', alpha=0.25)

	#handles, labels = plt.gca().get_legend_handles_labels() # get existing handles and labels
	#star_ns = Line2D([], [], marker='s', linestyle='', color='black', label=r'$k_{\star} =$ 0.05 Mpc$^{-1}$', markersize=12)
	#star_r = Line2D([], [], marker='*', linestyle='', color='black', label=r'$k_{\star} =$ 0.002 Mpc$^{-1}$', markersize=10)
	#handles.append(star_ns)
	#handles.append(star_r)
	#labels.append(r'$k_{\star} =$ 0.05 Mpc$^{-1}$')
	#labels.append(r'$k_{\star} =$ 0.002 Mpc$^{-1}$')
	#ax.legend(handles, labels, loc='lower left')
	ax.legend(loc='lower left')
	xPlanck, yPlanck = PlanckCoor
	ax.text(x=xPlanck, y=yPlanck, s="Planck", color='blue', fontsize=20)
	if Nstar50 != None :
		xNstar50, yNstar50 = Nstar50
		ax.text(x=xNstar50, y=yNstar50, s=r"$N_{\star} = 50$", color='red', fontsize=20)
	if Nstar60 != None :
		xNstar60, yNstar60 = Nstar60
		ax.text(x=xNstar60, y=yNstar60, s=r"$N_{\star} = 60$", color='red', fontsize=20)
	ob = offsetbox.AnchoredText(s=text, loc='upper left', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)
	ax.grid(False)
	fig.show()

	return


def rnsPlot1Y0(outDirName, title, fname1='nsr', fname2=None, xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, PlanckCoor=(0.967,0.018), Nstar50=None, Nstar60=None, textModel=None, nfiles=None) :
	"""
	Plot of r vs n_s. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* fname2 = name of the file containing the data for comparison (optional)
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red", "pink", "yellow", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	nsrList = [] # for k* = 0.002/Mpc
	
	ns_ns_starList = [] # for k* = 0.05/Mpc
	r_ns_starList = []
	ns_r_starList = [] # for k* = 0.002/Mpc
	r_r_starList = []

	PlanckData_TTEEBK15BAO_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15BAO_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15BAO_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15BAO_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_1s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEEBK15_ns_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEEBK15_r_2s = np.loadtxt(fname='PlanckData/TTEEBK15_Planck_2sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_1s = np.loadtxt(fname='PlanckData/TTEE_Planck_1sig.dat', delimiter=' ', usecols = 1)
	PlanckData_TTEE_ns_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 0)
	PlanckData_TTEE_r_2s = np.loadtxt(fname='PlanckData/TTEE_Planck_2sig.dat', delimiter=' ', usecols = 1)

	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		nsr_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
		nsrList.append(nsr_i)
		if fname2 != None :
			nsrNew_i = np.loadtxt(fname=outDirName + '/' + fname2 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1)
			nsrNewList.append(nsrNew_i)

		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "ns_star_ns = " in line:
				value = line[13:-1]
				ns_star_i = float(value)
				ns_ns_starList.append(ns_star_i)
			if "r_star_ns = " in line:
				value = line[12:-1]
				r_star_i = float(value)
				r_ns_starList.append(r_star_i)
			if "ns_star_r = " in line:
				value = line[12:-1]
				ns_star_i = float(value)
				ns_r_starList.append(ns_star_i)
			if "r_star_r = " in line:
				value = line[11:-1]
				r_star_i = float(value)
				r_r_starList.append(r_star_i)

	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1
	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\xi$ = " + str(xiList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i]) + "\n" + r"$X_0$ = " + str(x0List[i])
		

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$n_s$', fontsize=18)
	ax.set_ylabel(r'$r$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$Y_0 =$ " + str(y0List[i])
		nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(nsns, r_sg, color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		#ax.plot(ns_ns_starList[i], r_ns_starList[i], color = colorList[i], marker = 's', markersize=12)
		ax.plot(ns_r_starList[i], r_r_starList[i], color = colorList[i], marker = '*', markersize=10)

	if fname2 != None :
		for i in range(0, num) :
			ax.plot(nsrNewList[i][:,2], nsrNewList[i][:,3], color = colorList[i], linestyle = '--', marker = '')	

	ax.plot(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color ='blue', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_1s, PlanckData_TTEEBK15BAO_r_1s, color='blue', alpha=0.1)
	ax.plot(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color ='cyan', linestyle = '-', lw=1)
	ax.fill_between(PlanckData_TTEEBK15BAO_ns_2s, PlanckData_TTEEBK15BAO_r_2s, color='cyan', alpha=0.15)

#	ax.plot(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color ='red', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_1s, PlanckData_TTEEBK15_r_1s, color='red', alpha=0.2)
#	ax.plot(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color ='orange', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEEBK15_ns_2s, PlanckData_TTEEBK15_r_2s, color='orange', alpha=0.25)

#	ax.plot(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color ='black', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_1s, PlanckData_TTEE_r_1s, color='black', alpha=0.2)
#	ax.plot(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color ='gray', linestyle = '-', lw=1)
#	ax.fill_between(PlanckData_TTEE_ns_2s, PlanckData_TTEE_r_2s, color='gray', alpha=0.25)

	#handles, labels = plt.gca().get_legend_handles_labels() # get existing handles and labels
	#star_ns = Line2D([], [], marker='s', linestyle='', color='black', label=r'$k_{\star} =$ 0.05 Mpc$^{-1}$', markersize=12)
	#star_r = Line2D([], [], marker='*', linestyle='', color='black', label=r'$k_{\star} =$ 0.002 Mpc$^{-1}$', markersize=10)
	#handles.append(star_ns)
	#handles.append(star_r)
	#labels.append(r'$k_{\star} =$ 0.05 Mpc$^{-1}$')
	#labels.append(r'$k_{\star} =$ 0.002 Mpc$^{-1}$')
	#ax.legend(handles, labels, loc='upper right')
	ax.legend(loc='lower right')
	xPlanck, yPlanck = PlanckCoor
	ax.text(x=xPlanck, y=yPlanck, s="Planck", color='blue', fontsize=20)
	if Nstar50 != None :
		xNstar50, yNstar50 = Nstar50
		ax.text(x=xNstar50, y=yNstar50, s=r"$N_{\star} = 50$", color='red', fontsize=20)
	if Nstar60 != None :
		xNstar60, yNstar60 = Nstar60
		ax.text(x=xNstar60, y=yNstar60, s=r"$N_{\star} = 60$", color='red', fontsize=20)
	ob = offsetbox.AnchoredText(s=text, loc='lower left', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)

	if textModel != None :
		ob1 = offsetbox.AnchoredText(s=textModel, loc='upper right', pad=0.3, borderpad=0.85, frameon=False, prop=dict(color='black', size=15))
		ax.add_artist(ob1)

	ax.grid(False)
	fig.show()

	return



def VePlotX0(outDirName, title, fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel=None, nfiles=None) :
	"""
	Plot of VE vs NN. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	##colorList = ["brown", "black", "red", "green", "yellow", "magenta", "blue", "pink"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red", "pink", "yellow", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	NList = []
	xList = []
	yList = []

	N_starList = []
	V_starList = []
	x_starList = []
	y_starList = []
	
	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		N_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=0)
		NList.append(N_i)
		x_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=1)
		xList.append(x_i)
		y_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=3)
		yList.append(y_i)

		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "Nstar = " in line:
				value = line[8:-1]
				N_star_i = float(value)
				N_starList.append(N_star_i)
			if "Vstar = " in line:
				value = line[8:-1]
				V_star_i = float(value)
				V_starList.append(V_star_i)
			if "xstar = " in line and "pixstar = " not in line:
				value = line[8:-1]
				x_star_i = float(value)
				x_starList.append(x_star_i)
			if "ystar = " in line and "piystar = " not in line:
				value = line[8:-1]
				y_star_i = float(value)
				y_starList.append(y_star_i)


	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1
	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\xi$ = " + str(xiList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i]) + "\n" + r"$Y_0$ = " + str(y0List[i])
		

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$N$', fontsize=18)
	ax.set_ylabel(r'$V_E$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$X_0 =$ " + str(x0List[i])
		Ve_i = myTools.VE(x=xList[i], y=yList[i], mllpxilpp=(mList[i], lList[i], lpList[i], xiList[i], lppList[i]))
		Vstar_i = myTools.VE(x=x_starList[i], y=y_starList[i], mllpxilpp=(mList[i], lList[i], lpList[i], xiList[i], lppList[i]))
		#nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(NList[i], Ve_i, color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		ax.plot(N_starList[i], V_starList[i], color = colorList[i], marker = '*', markersize=10)
		ax.plot(N_starList[i], Vstar_i, color = colorList[i], marker = '*', markersize=10)
		
	ax.legend(loc='lower right')
	ob = offsetbox.AnchoredText(s=text, loc='upper right', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)
	if textModel != None :
		ob1 = offsetbox.AnchoredText(s=textModel, loc='lower left', pad=0.3, borderpad=0.85, frameon=False, prop=dict(color='black', size=15))
		ax.add_artist(ob1)
	ax.grid(False)
	fig.show()

	return


def VePlotY0(outDirName, title, fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel=None, nfiles=None) :
	"""
	Plot of VE vs NN. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	##colorList = ["brown", "black", "red", "green", "yellow", "magenta", "blue", "pink"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red", "pink", "yellow", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	NList = []
	xList = []
	yList = []

	N_starList = []
	V_starList = []
	x_starList = []
	y_starList = []
	
	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		N_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=0)
		NList.append(N_i)
		x_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=1)
		xList.append(x_i)
		y_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=3)
		yList.append(y_i)

		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "Nstar = " in line:
				value = line[8:-1]
				N_star_i = float(value)
				N_starList.append(N_star_i)
			if "Vstar = " in line:
				value = line[8:-1]
				V_star_i = float(value)
				V_starList.append(V_star_i)
			if "xstar = " in line and "pixstar = " not in line:
				value = line[8:-1]
				x_star_i = float(value)
				x_starList.append(x_star_i)
			if "ystar = " in line and "piystar = " not in line:
				value = line[8:-1]
				y_star_i = float(value)
				y_starList.append(y_star_i)
		
		
	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1

	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\xi$ = " + str(xiList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i]) + "\n" + r"$X_0$ = " + str(x0List[i])
		

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$N$', fontsize=18)
	ax.set_ylabel(r'$V_E$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$Y_0 =$ " + str(y0List[i])
		Ve_i = myTools.VE(x=xList[i], y=yList[i], mllpxilpp=(mList[i], lList[i], lpList[i], xiList[i], lppList[i]))
		Vstar_i = myTools.VE(x=x_starList[i], y=y_starList[i], mllpxilpp=(mList[i], lList[i], lpList[i], xiList[i], lppList[i]))
		#nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(NList[i], Ve_i, color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		ax.plot(N_starList[i], V_starList[i], color = colorList[i], marker = '*', markersize=10)
		ax.plot(N_starList[i], Vstar_i, color = colorList[i], marker = '*', markersize=10)
		
	ax.legend(loc='lower right')
	ob = offsetbox.AnchoredText(s=text, loc='upper right', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)
	if textModel != None :
		ob1 = offsetbox.AnchoredText(s=textModel, loc='lower left', pad=0.3, borderpad=0.85, frameon=False, prop=dict(color='black', size=15))
		ax.add_artist(ob1)
	ax.grid(False)
	fig.show()

	return


def xyPlotX0(outDirName, title, fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel=None, nfiles=None) :
	"""
	Plot of VE vs NN. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	##colorList = ["brown", "black", "red", "green", "yellow", "magenta", "blue", "pink"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red", "pink", "yellow", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	xList = []
	yList = []

	x_starList = []
	y_starList = []
	
	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		x_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=1)
		xList.append(x_i)
		y_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=3)
		yList.append(y_i)


		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "xstar = " in line and "pixstar = " not in line:
				value = line[8:-1]
				x_star_i = float(value)
				x_starList.append(x_star_i)
			if "ystar = " in line and "piystar = " not in line:
				value = line[8:-1]
				y_star_i = float(value)
				y_starList.append(y_star_i)


	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1
	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\xi$ = " + str(xiList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i]) + "\n" + r"$Y_0$ = " + str(y0List[i])
		

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$X$', fontsize=18)
	ax.set_ylabel(r'$Y$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$X_0 =$ " + str(x0List[i])
		ax.plot(xList[i], yList[i], color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		ax.plot(x_starList[i], y_starList[i], color = colorList[i], marker = '*', markersize=10)
		
	ax.legend(loc='upper left')
	ob = offsetbox.AnchoredText(s=text, loc='upper right', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)
	if textModel != None :
		ob1 = offsetbox.AnchoredText(s=textModel, loc='lower right', pad=0.3, borderpad=0.85, frameon=False, prop=dict(color='black', size=15))
		ax.add_artist(ob1)
	ax.grid(False)
	fig.show()

	return


def xyPlotY0(outDirName, title, fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel=None, nfiles=None) :
	"""
	Plot of VE vs NN. 
	Input:
	* outDirName = name of the output directory
	* fname1 = name of the file containing the data
	* xlim = x plot-limits (optional) --> ref: (0.94, 1.0); (0.94, 1.02)
	* ylim = y plot-limits (optional) --> ref: (0.0045, 0.011); (0.004, 0.018)
	* LogXaxis = Boolean used if we want to use log scale for x-axis (optional)
	* LogYaxis = Boolean used if we want to use log scale for y-axis (optional)
	* nfiles = number of files to consider
	Output: plot
	"""

	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	#colorList = ["red", "green", "black", "brown", "yellow", "blue", "magenta"]
	##colorList = ["brown", "black", "red", "green", "yellow", "magenta", "blue", "pink"]
	colorList = ["magenta", "brown", "black", "orange", "green", "purple", "red", "pink", "yellow", "blue"]
	mList = []
	lList = []
	lpList = []
	xiList = []
	lppList = []
	x0List = []
	y0List = []
	xList = []
	yList = []

	x_starList = []
	y_starList = []
	
	for i in range(0, len(param2D)) :
		line_i = param2D[i]
		m_i = line_i[0] # potential parameter
		l_i = line_i[1] # potential parameter
		lp_i = line_i[2] # potential parameter
		xi_i = line_i[3] # potential parameter
		lpp_i = line_i[4]
		x0_i = line_i[5] # initial value of x
		y0_i = line_i[6] # initial value of y

		mList.append(m_i)
		lList.append(l_i)
		lpList.append(lp_i)
		xiList.append(xi_i) # --> ref: [0.2, 0.225, 0.25, 0.275, 0.3, 0.35]; [1.0e-3, 2.5e-3, 5.4e-3, 1.0e-2, 2.0e-2]
		lppList.append(lpp_i)
		x0List.append(x0_i)
		y0List.append(y0_i)

		x_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=1)
		xList.append(x_i)
		y_i = np.loadtxt(fname=outDirName + '/' + fname1 + '_'+ str(i) + '.dat', delimiter='\t', skiprows=1, usecols=3)
		yList.append(y_i)


		with open(outDirName + '/' + 'Nstar' + '_'+ str(i) + '.dat', "r") as f:
			searchlines = f.readlines()
		for line in searchlines:
			if "xstar = " in line and "pixstar = " not in line:
				value = line[8:-1]
				x_star_i = float(value)
				x_starList.append(x_star_i)
			if "ystar = " in line and "piystar = " not in line:
				value = line[8:-1]
				y_star_i = float(value)
				y_starList.append(y_star_i)
		
		
	idNum=0
	for i in range(0, len(mList)) :
		if mList[i] == mList[0] and lList[i] == lList[0] and lpList[i] == lpList[0] and lppList[i] == lppList[0] :
			idNum += 1

	if idNum == len(mList) :
		text = r"$m$ = " + str(mList[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(lList[i]) + "\n" + r"$\lambda'$ = " + str(lpList[i]) + "\n" + r"$\xi$ = " + str(xiList[i]) + "\n" + r"$\lambda''$ = " + str(lppList[i]) + "\n" + r"$X_0$ = " + str(x0List[i])
		

	fig = plt.figure(num=title, figsize=(13, 7.2), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel(r'$X$', fontsize=18)
	ax.set_ylabel(r'$Y$', fontsize=18)
	if xlim != None :
		xMin, xMax = xlim
		ax.set_xlim(xMin, xMax)
	if ylim != None :
		yMin, yMax = ylim
		ax.set_ylim(yMin, yMax)
	if LogXaxis :
		ax.set_xscale('log')
	else :
		ax.set_xscale('linear')
	if LogYaxis :
		ax.set_yscale('log')
	else :
		ax.set_yscale('linear')

	if nfiles != None and nfiles < len(param2D) :
		num = nfiles
	else :
		num = len(param2D)
	for i in range(0, num) :
		label_i = r"$Y_0 =$ " + str(y0List[i])
		#nsns, r_sg = InterpAndSmooth(x=nsrList[i][:,2], y=nsrList[i][:,3], intkind='linear', Nsteps=1e3)
		ax.plot(xList[i], yList[i], color = colorList[i], linestyle = '-', lw=2, marker = '', label = label_i)
		ax.plot(x_starList[i], y_starList[i], color = colorList[i], marker = '*', markersize=10)
		
	ax.legend(loc='upper left')
	ob = offsetbox.AnchoredText(s=text, loc='upper right', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax.add_artist(ob)
	if textModel != None :
		ob1 = offsetbox.AnchoredText(s=textModel, loc='lower right', pad=0.3, borderpad=0.85, frameon=False, prop=dict(color='black', size=15))
		ax.add_artist(ob1)
	ax.grid(False)
	fig.show()

	return





if __name__ == '__main__':

	#outDirName = 'outRK45'
	#outDirName = 'outLSODA'
	#outDirName = 'outLSODA_smallxi'
	#outDirName = 'outLRK45_smallxi'
	#outDirName = 'outLSODA_transfer'
	outDirName = 'outRK45_transfer'
	#outDirName = 'fig1outRK'
	#outDirName = 'fig11outRK'
	
	#outDirName = 'mod1Y0fix'
	#outDirName = 'mod1X0fix'
	#outDirName = 'mod2Y0fix'
	#outDirName = 'mod2X0fix'
	##outDirName = 'out1'
	##textModel="Model 1"
	##outDirName = 'out2'
	##textModel="Model 2"
	#outDirName = 'out'
	i = 0


	param2D = myTools.Create2DListFromTXT(NameTXTinput=outDirName + '/' + 'inf_param.dat')
	m = [line[0] for line in param2D]
	l = [line[1] for line in param2D]
	lp = [line[2] for line in param2D]
	xi = [line[3] for line in param2D]
	lpp = [line[4] for line in param2D]
	x0 = [line[5] for line in param2D]
	y0 = [line[6] for line in param2D]

	text = r"$m$ = " + str(m[i]) + r" $m_{Pl}$" + "\n" + r"$\lambda$ = " + str(l[i]) + "\n" + r"$\lambda'$ = " + str(lp[i]) + "\n" + r"$\xi$ = "+ str(xi[i]) + "\n" + r"$\lambda''$ = " + str(lpp[i]) + "\n" + r"$X_0$ ="  + str(x0[i]) + "\n" + r"$Y_0$ = " + str(y0[i])
	print "m = " + str(m[i]) + "\n" + "l = " + str(l[i]) + "\n" + "lp = " + str(lp[i]) + "\n" + "xi = "+ str(xi[i]) + "\n" + "lpp = " + str(lpp[i]) + "\n" + "x0 = " + str(x0[i]) + "\n" + "y0 = " + str(y0[i])


#	rnsPlot(outDirName, title='nsr_two', xlim=(0.94, 1.0), nfiles=2)
#	rnsPlot(outDirName, title='nsr_six', xlim=(0.95, 0.98), ylim=(0.007,0.0165), nfiles=6)
	#rnsPlot1(outDirName, title='nsr_seven', xlim=(0.958, 0.976), ylim=(0.005,0.07), nfiles=7)

	##rnsPlot1X0(outDirName, title='nsr_seven', xlim=(0.9565, 0.9753), ylim=(0.005,0.02), PlanckCoor=(0.967,0.0187), Nstar50=(0.957,0.0172), Nstar60=(0.9658,0.011), textModel='Model 1', nfiles=8) # Model 1
	##rnsPlot1X0(outDirName, title='nsr_seven', xlim=(0.958, 0.97), ylim=(0.005,0.01), Nstar50=(0.9585,0.0088), Nstar60=(0.9665,0.006), textModel='Model 2', nfiles=8) # Model 2
	#rnsPlot1Y0(outDirName, title='nsr_seven', xlim=(0.956, 0.985), ylim=(0.005,0.022), PlanckCoor=(0.967,0.0187), Nstar50=(0.957,0.0172), Nstar60=(0.9658,0.011), textModel='Model 1', nfiles=8) # Model 1
	#rnsPlot1Y0(outDirName, title='nsr_seven', xlim=(0.9535, 0.9755), ylim=(0.005,0.012), PlanckCoor=(0.9665,0.011), Nstar50=(0.954,0.01055), Nstar60=(0.9665,0.0057), textModel='Model 2', nfiles=8) # Model 2

	##VePlotX0(outDirName, title='Ve_cmp_Y0_fix_mod1', fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel='Model 1', nfiles=8) # Model 1
	#VePlotX0(outDirName, title='Ve_cmp_Y0_fix_mod2', fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel='Model 2', nfiles=8) # Model 2
	##VePlotY0(outDirName, title='Ve_cmp_X0_fix_mod1', fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel='Model 1', nfiles=8) # Model 1
	#VePlotY0(outDirName, title='Ve_cmp_X0_fix_mod2', fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel='Model 2', nfiles=8) # Model 2
	
	##xyPlotX0(outDirName, title='xy_cmp_Y0_fix_mod1', fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel='Model 1', nfiles=8) # Model 1
	#xyPlotX0(outDirName, title='xy_cmp_Y0_fix_mod2', fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel='Model 2', nfiles=8) # Model 2
	##xyPlotY0(outDirName, title='xy_cmp_X0_fix_mod1', fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel='Model 1', nfiles=8) # Model 1
	#xyPlotY0(outDirName, title='xy_cmp_X0_fix_mod2', fname1='inf', xlim=None, ylim=None, LogXaxis=False, LogYaxis=False, textModel='Model 2', nfiles=8) # Model 2


	###
	NN = np.loadtxt(fname=outDirName + '/inf_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 0)
	x = np.loadtxt(fname=outDirName + '/inf_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 1)
	pix = np.loadtxt(fname=outDirName + '/inf_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	y = np.loadtxt(fname=outDirName + '/inf_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)
	piy = np.loadtxt(fname=outDirName + '/inf_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 4)
	psi = np.loadtxt(fname=outDirName + '/inf_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 5)

	fig0 = plt.figure(num='xy', figsize=(13, 7.2), dpi=100)
	ax0 = fig0.add_subplot(111)
	ax0.set_xlabel(r'$N$', fontsize=17)
	ax0.set_ylabel('', fontsize=17)
	ax0.set_xscale('linear')
	ax0.set_yscale('linear')
	ax0.plot(NN, x, color = "red", linestyle = '-', marker = '', label=r'$X$')
	ax0.plot(NN, y, color = "blue", linestyle = '-', marker = '', label=r'$Y$')
	ax0.legend(loc='best', prop={'size': 16})
	ax0.xaxis.set_minor_locator(AutoMinorLocator())
	ax0.yaxis.set_minor_locator(AutoMinorLocator())
	ob1 = offsetbox.AnchoredText(s=textModel, loc='center right', pad=0.3, borderpad=0.85, frameon=False, prop=dict(color='black', size=15))
	ax0.add_artist(ob1)
	ax0.grid(False)

	fig1 = plt.figure(num='pixpiy', figsize=(13, 7.2), dpi=100)
	ax1 = fig1.add_subplot(111)
	ax1.set_xlabel(r'$N$', fontsize=17)
	ax1.set_ylabel('', fontsize=17)
	ax1.set_xscale('linear')
	ax1.set_yscale('linear')
	ax1.plot(NN, pix, color = "red", linestyle = '-', marker = '', label=r'$\pi x$')
	ax1.plot(NN, piy, color = "blue", linestyle = '-', marker = '', label=r'$\pi y$')
	ax1.legend(loc='best', prop={'size': 16})
	ax1.grid(False)

	fig2 = plt.figure(num='psi', figsize=(13, 7.2), dpi=100)
	ax2 = fig2.add_subplot(111)
	ax2.set_xlabel(r'$N$', fontsize=17)
	ax2.set_ylabel(r'$\psi$', fontsize=17)
	ax2.set_xscale('linear')
	ax2.set_yscale('linear')
	ax2.plot(NN, psi, color = "red", linestyle = '-', marker = '')
	#ax2.legend(loc='lower right', prop={'size': 16})
	ax2.grid(False)


	###
	H = np.loadtxt(fname=outDirName + '/nbni_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 1)
	nb = np.loadtxt(fname=outDirName + '/nbni_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	ni = np.loadtxt(fname=outDirName + '/nbni_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)
	log10nbni = np.loadtxt(fname=outDirName + '/nbni_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 4)

	fig3 = plt.figure(num='log10nbni', figsize=(13, 7.2), dpi=100)
	ax3 = fig3.add_subplot(111)
	ax3.set_xlabel(r'$N$', fontsize=17)
	ax3.set_ylabel(r'$\log_{10} (n_B/n_{\phi})$', fontsize=17)
	ax3.set_xscale('linear')
	ax3.set_yscale('linear')
	ax3.plot(NN, log10nbni, color = "red", linestyle = '-', marker = '', label=r"$\lambda' = $" + str(lp[i]))
	ax3.legend(loc='best', prop={'size': 16})
	ax3.xaxis.set_minor_locator(AutoMinorLocator())
	ax3.yaxis.set_minor_locator(AutoMinorLocator())
	ax3.grid(False)

	fig4 = plt.figure(num='nbni_cmp', figsize=(13, 7.2), dpi=100)
	ax4 = fig4.add_subplot(111)
	ax4.set_xlabel(r'$N$', fontsize=17)
	ax4.set_ylabel('', fontsize=17)
	ax4.set_xscale('linear')
	ax4.set_yscale('linear')
	ax4.plot(NN, nb, color = "red", linestyle = '-', marker = '', label=r'$n_B$')
	ax4.plot(NN, ni, color = "blue", linestyle = '-', marker = '', label=r'$n_{\phi}$')
	ax4.plot(NN, nb/ni, color = "green", linestyle = '-', marker = '', label=r'$n_B/n_{\phi}$')
	ax4.plot(NN, np.abs(nb/ni), color = "orange", linestyle = '--', marker = '', label=r'$|n_B/n_{\phi}|$')
	ax4.plot(NN, np.zeros(len(NN)), color = "black", linestyle = '-', marker = '')
	ax4.legend(loc='best', prop={'size': 16})
	ax4.xaxis.set_minor_locator(AutoMinorLocator())
	ax4.yaxis.set_minor_locator(AutoMinorLocator())
	ax4.grid(False)

	fig5 = plt.figure(num='H', figsize=(13, 7.2), dpi=100)
	ax5 = fig5.add_subplot(111)
	ax5.set_xlabel(r'$N$', fontsize=17)
	ax5.set_ylabel(r'$H$', fontsize=17)
	ax5.set_xscale('linear')
	ax5.set_yscale('linear')
	ax5.plot(NN, H, color = "blue", linestyle = '-', marker = '')
	#ax5.legend(loc='best', prop={'size': 16})
	ax5.xaxis.set_minor_locator(AutoMinorLocator())
	ax5.yaxis.set_minor_locator(AutoMinorLocator())
	ax5.grid(False)


	###
	Ud = np.loadtxt(fname=outDirName + '/theta_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 1)
	Vd = np.loadtxt(fname=outDirName + '/theta_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	theta = np.loadtxt(fname=outDirName + '/theta_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)
	costh = np.loadtxt(fname=outDirName + '/theta_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 4)
	sinth = np.loadtxt(fname=outDirName + '/theta_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 5)
	thd = np.loadtxt(fname=outDirName + '/theta_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 6)
	thd_alt = np.loadtxt(fname=outDirName + '/theta_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 7)
		
	fig6 = plt.figure(num='thetad', figsize=(13, 7.2), dpi=100)
	ax6 = fig6.add_subplot(111)
	ax6.set_xlabel(r'$N$', fontsize=17)
	ax6.set_ylabel(r'$\dot{\phi}$', fontsize=17)
	ax6.set_xscale('linear')
	ax6.set_yscale('linear')
	ax6.set_xlim(0, 80.0)
	ax6.plot(NN, thd, color = "red", linestyle = '-', marker = '', label=r'$\frac{H}{1+\tan^2{\phi}} \frac{d\tan{\phi}}{dN}$')
	ax6.plot(NN, thd_alt, color = "blue", linestyle = '-', marker = '', label=r'$\frac{H (\dot{U} \frac{d\dot{V}}{dN} - \dot{V} \frac{d\dot{U}}{dN})}{\dot{U}^2 + \dot{V}^2}$')
	ax6.legend(loc='best', prop={'size': 16})
	ax6.grid(False)

	fig7 = plt.figure(num='sincosth', figsize=(13, 7.2), dpi=100)
	ax7 = fig7.add_subplot(111)
	ax7.set_xlabel(r'$N$', fontsize=17)
	ax7.set_ylabel('', fontsize=17)
	ax7.set_xscale('linear')
	ax7.set_yscale('linear')
	ax7.plot(NN, costh, color = "red", linestyle = '-', marker = '', label=r'$\cos{\phi}$')
	ax7.plot(NN, sinth, color = "blue", linestyle = '-', marker = '', label=r'$\sin{\phi}$')
	ax7.legend(loc='best', prop={'size': 16})
	ax7.grid(False)

	fig8 = plt.figure(num='theta', figsize=(13, 7.2), dpi=100)
	ax8 = fig8.add_subplot(111)
	ax8.set_xlabel(r'$N$', fontsize=17)
	ax8.set_ylabel(r'$\tan{\phi}$', fontsize=17)
	ax8.set_xscale('linear')
	ax8.set_yscale('linear')
	ax8.plot(NN, np.tan(theta), color = "red", linestyle = '-', marker = '')
	#ax6.legend(loc='lower right', prop={'size': 16})
	ax8.grid(False)

	fig9 = plt.figure(num='UdVd', figsize=(13, 7.2), dpi=100)
	ax9 = fig9.add_subplot(111)
	ax9.set_xlabel(r'$N$', fontsize=17)
	ax9.set_ylabel('', fontsize=17)
	ax9.set_xscale('linear')
	ax9.set_yscale('linear')
	ax9.plot(NN, Ud, color = "red", linestyle = '-', marker = '', label=r'$\dot{U}$')
	ax9.plot(NN, Vd, color = "blue", linestyle = '-', marker = '', label=r'$\dot{V}$')
	ax9.legend(loc='best', prop={'size': 16})
	ax9.grid(False)


	###
	epsU = np.loadtxt(fname=outDirName + '/UV_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 4)
	epsV = np.loadtxt(fname=outDirName + '/UV_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 5)
	etaUU = np.loadtxt(fname=outDirName + '/UV_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 6)
	etaUV = np.loadtxt(fname=outDirName + '/UV_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 7)
	etaVV = np.loadtxt(fname=outDirName + '/UV_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 9)

	fig10 = plt.figure(num='epsUepsV', figsize=(13, 7.2), dpi=100)
	ax10 = fig10.add_subplot(111)
	ax10.set_xlabel(r'$N$', fontsize=17)
	ax10.set_ylabel('', fontsize=17)
	ax10.set_xscale('linear')
	ax10.set_yscale('linear')
	ax10.plot(NN, epsU, color = "red", linestyle = '-', marker = '', label=r'$\epsilon_U$')
	ax10.plot(NN, epsV, color = "blue", linestyle = '-', marker = '', label=r'$\epsilon_V$')
	ax10.legend(loc='best', prop={'size': 16})
	ax10.grid(False)

	fig11 = plt.figure(num='etaUV', figsize=(13, 7.2), dpi=100)
	ax11 = fig11.add_subplot(111)
	ax11.set_xlabel(r'$N$', fontsize=17)
	ax11.set_ylabel('', fontsize=17)
	ax11.set_xscale('linear')
	ax11.set_yscale('linear')
	ax11.plot(NN, etaUU, color = "red", linestyle = '-', marker = '', label=r'$\eta_{UU}$')
	ax11.plot(NN, etaUV, color = "green", linestyle = '-', marker = '', label=r'$\eta_{UV}$')
	ax11.plot(NN, etaVV, color = "blue", linestyle = '-', marker = '', label=r'$\eta_{VV}$')
	ax11.legend(loc='best', prop={'size': 16})
	ax11.grid(False)


	###
	cosg_Matteo = np.loadtxt(fname=outDirName + '/gamma_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	sing_Matteo = np.loadtxt(fname=outDirName + '/gamma_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)
	cosg_Jim = np.loadtxt(fname=outDirName + '/gamma_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 5)
	sing_Jim = np.loadtxt(fname=outDirName + '/gamma_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 6)

	fig12 = plt.figure(num='gamma', figsize=(13, 7.2), dpi=100)
	ax12 = fig12.add_subplot(111)
	ax12.set_xlabel(r'$N$', fontsize=17)
	ax12.set_ylabel('', fontsize=17)
	ax12.set_xscale('linear')
	ax12.set_yscale('linear')
	ax12.plot(NN, np.abs(cosg_Matteo), color = "red", linestyle = '-', marker = '', label=r'$\cos{\gamma_{Matteo}}$')
	ax12.plot(NN, np.abs(sing_Matteo), color = "blue", linestyle = '-', marker = '', label=r'$\sin{\gamma_{Matteo}}$')
	ax12.plot(NN, np.abs(cosg_Jim), color = "pink", linestyle = '-', marker = '', label=r'$\cos{\gamma_{Jim}}$')
	ax12.plot(NN, np.abs(sing_Jim), color = "green", linestyle = '-', marker = '', label=r'$\sin{\gamma_{Jim}}$')
	ax12.legend(loc='best', prop={'size': 16})
	ax12.grid(False)

	tang_Matteo = np.loadtxt(fname=outDirName + '/gamma_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 1)
	tang_Jim = np.loadtxt(fname=outDirName + '/gamma_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 4)
	pio2 = np.full(shape=len(NN), fill_value=np.pi/2)
	
	fig121 = plt.figure(num='gamma11', figsize=(13, 7.2), dpi=100)
	ax121 = fig121.add_subplot(111)
	ax121.set_xlabel(r'$N$', fontsize=17)
	ax121.set_ylabel('', fontsize=17)
	ax121.set_xscale('linear')
	ax121.set_yscale('linear')
	ax121.plot(NN, np.arccos(np.abs(cosg_Matteo)), color = "red", linestyle = '-', marker = '', label=r'$\gamma_{Matteo}$')
	ax121.plot(NN, np.arcsin(np.abs(sing_Matteo)), color = "blue", linestyle = '-', marker = '', label=r'$\gamma_{Matteo}$')
	ax121.plot(NN, np.arccos(np.abs(cosg_Jim)), color = "pink", linestyle = '-', marker = '', label=r'$\gamma_{Jim}$')
	ax121.plot(NN, np.arcsin(np.abs(sing_Jim)), color = "green", linestyle = '-', marker = '', label=r'$\gamma_{Jim}$')
	ax121.plot(NN, pio2, color = "black", linestyle = '--', marker = '')
	ax121.legend(loc='best', prop={'size': 16})
	ax121.grid(False)


	###
	Ad = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 1)
	VeA = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	VeS_def = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)
	VeS = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 4)
	VeS_alt = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 5)
	epsA = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 6)
	epsS = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 7)
	etaAA = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 8)
	etaAS = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 9)
	etaSS = np.loadtxt(fname=outDirName + '/AS_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 10)

	Cem = 2.0 - np.log(2.0) - 0.5772
	cosD = -2*Cem*etaAS 

	fig13 = plt.figure(num='VeAVeS', figsize=(13, 7.2), dpi=100)
	ax13 = fig13.add_subplot(111)
	ax13.set_xlabel(r'$N$', fontsize=17)
	ax13.set_ylabel('', fontsize=17)
	ax13.set_xscale('linear')
	ax13.set_yscale('linear')
	ax13.plot(NN, VeA, color = "red", linestyle = '-', marker = '', label=r'$V_{E, A} = \cos{\phi} V_{E, U} + \sin{\phi} V_{E, V}$')
	ax13.plot(NN, VeS_def, color = "blue", linestyle = '-', marker = '', label=r'$V_{E, S} = -\sin{\phi} V_{E, U} + \cos{\phi} V_{E, V}$')
	ax13.legend(loc='best', prop={'size': 16})
	ax13.grid(False)

	fig14 = plt.figure(num='VeS', figsize=(13, 7.2), dpi=100)
	ax14 = fig14.add_subplot(111)
	ax14.set_xlabel(r'$N$', fontsize=17)
	ax14.set_ylabel(r'$V_{E, S}$', fontsize=17)
	ax14.set_xscale('linear')
	ax14.set_yscale('linear')
	ax14.plot(NN, VeS_def, color = "red", linestyle = '-', marker = '', label=r'$-\sin{\phi} V_{E, U} + \cos{\phi} V_{E, V}$')
	ax14.plot(NN, VeS, color = "blue", linestyle = '-', marker = '', label=r'$-\dot{A} \dot{\phi}$')
	ax14.plot(NN, VeS_alt, color = "green", linestyle = '-', marker = '', label=r'$-\dot{A} \dot{\phi}_{alt}$')
	ax14.legend(loc='best', prop={'size': 16})
	ax14.grid(False)

	fig15 = plt.figure(num='epsAepsS', figsize=(13, 7.2), dpi=100)
	ax15 = fig15.add_subplot(111)
	ax15.set_xlabel(r'$N$', fontsize=17)
	ax15.set_ylabel('', fontsize=17)
	ax15.set_xscale('linear')
	ax15.set_yscale('log')
	ax15.set_xlim(0.0,64.5)
	ax15.plot(NN, epsA, color = "red", linestyle = '-', marker = '', label=r'$\epsilon_A$')
	ax15.plot(NN, epsS, color = "blue", linestyle = '-', marker = '', label=r'$\epsilon_S$')
	ax15.legend(loc='best', prop={'size': 16})
	ax15.grid(False)

	fig16 = plt.figure(num='etaAS', figsize=(13, 7.2), dpi=100)
	ax16 = fig16.add_subplot(111)
	ax16.set_xlabel(r'$N$', fontsize=17)
	ax16.set_ylabel('', fontsize=17)
	ax16.set_xscale('linear')
	ax16.set_yscale('linear')
	ax16.set_xlim(0.0,64.5)
	ax16.plot(NN, etaAA, color = "red", linestyle = '-', marker = '', label=r'$\eta_{AA}$')
	ax16.plot(NN, etaAS, color = "green", linestyle = '-', marker = '', label=r'$\eta_{AS}$')
	ax16.plot(NN, etaSS, color = "blue", linestyle = '-', marker = '', label=r'$\eta_{SS}$')	
	ax16.legend(loc='best', prop={'size': 16})
	ax16.grid(False)

	fig17 = plt.figure(num='cosDelta', figsize=(13, 7.2), dpi=100)
	ax17 = fig17.add_subplot(111)
	ax17.set_xlabel(r'$N$', fontsize=17)
	ax17.set_ylabel(r'$cos{\Delta}$', fontsize=17)
	ax17.set_xscale('linear')
	ax17.set_yscale('linear')
	ax17.plot(NN, cosD, color = "red", linestyle = '-', marker = '')
	#ax17.legend(loc='best', prop={'size': 16})
	ax17.grid(False)


	###
	ns_old = np.loadtxt(fname=outDirName + '/nsrOLD_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	r_old = np.loadtxt(fname=outDirName + '/nsrOLD_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)

	fig18 = plt.figure(num='nsrOLD', figsize=(13, 7.2), dpi=100)
	ax18 = fig18.add_subplot(111)
	ax18.set_xlabel(r'$n_s$', fontsize=17)
	ax18.set_ylabel(r'$r$', fontsize=17)
	ax18.set_xscale('linear')
	ax18.set_yscale('linear')
	ax18.plot(ns_old, r_old, color = "red", linestyle = '-', marker = '', label=r"$\xi = $" + str(xi[i]))
	ax18.legend(loc='best', prop={'size': 16})
	ax18.grid(False)


	###
	Nk = np.loadtxt(fname=outDirName + '/nsr_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 1)
	ns = np.loadtxt(fname=outDirName + '/nsr_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	r = np.loadtxt(fname=outDirName + '/nsr_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)	

	fig19 = plt.figure(num='nsr', figsize=(13, 7.2), dpi=100)
	ax19 = fig19.add_subplot(111)
	ax19.set_xlabel(r'$n_s$', fontsize=17)
	ax19.set_ylabel(r'$r$', fontsize=17)
	ax19.set_xscale('linear')
	ax19.set_yscale('linear')
	ax19.plot(ns, r, color = "red", linestyle = '-', marker = '', label=r"$\xi = $" + str(xi[i]))
	ax19.legend(loc='best', prop={'size': 16})
	ax19.grid(False)

	Nk_ns_Jim = np.loadtxt(fname='Jim_1ns.dat', delimiter=' ', usecols = 0)
	log10_1ns_Jim = np.loadtxt(fname='Jim_1ns.dat', delimiter=' ', usecols = 1)
	#
	Nk_r_Jim = np.loadtxt(fname='Jim_r.dat', delimiter=' ', usecols = 0)
	log10_r_Jim = np.loadtxt(fname='Jim_r.dat', delimiter=' ', usecols = 1)

	fig20 = plt.figure(num='nsr_cmp', figsize=(13, 7.2), dpi=100)
	ax20 = fig20.add_subplot(111)
	ax20.set_xlabel(r'$N_{\star}$', fontsize=17)
	ax20.set_ylabel('', fontsize=17)
	ax20.set_xscale('linear')
	ax20.set_yscale('logit')
	ax20.plot(Nk, 1.0-ns, color = "red", linestyle = '-', marker = '', label=r"$log_{10}(1-ns)$ Matteo")
	ax20.plot(Nk_ns_Jim, log10_1ns_Jim, color = "orange", linestyle = '--', marker = '', label=r"$log_{10}(1-ns)$ Jim")
	ax20.plot(Nk, r, color = "blue", linestyle = '-', marker = '', label=r"$log_{10}(r)$ Matteo")
	ax20.plot(Nk_r_Jim, log10_r_Jim, color = "green", linestyle = '--', marker = '', label=r"$log_{10}(r)$ Jim")
	ax20.legend(loc='best', prop={'size': 16})
	ax20.grid(False)
	#ax20.yaxis.set_minor_formatter(NullFormatter())


	###
	NN_ustartot = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 0)
	dU_10 = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 1)
	dUp_10 = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	dV_10 =np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)
	dVp_10 = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 4)
	T_UU = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 5)
	T_VU = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 6)
	dA_10 = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 7)
	dS_10 = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 8)
	R_10 =np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 9)
	S_10 = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 10)
	T_RR = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 11)
	T_SR = np.loadtxt(fname=outDirName + '/dAdS10_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 12)
	
	dU_01 = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 1)
	dUp_01 = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 2)
	dV_01 =np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 3)
	dVp_01 = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 4)
	T_UV = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 5)
	T_VV = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 6)
	dA_01 = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 7)
	dS_01 = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 8)
	R_01 =np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 9)
	S_01 = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 10)
	T_RS = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 11)
	T_SS = np.loadtxt(fname=outDirName + '/dAdS01_' + str(i) + '.dat', delimiter='\t', skiprows=1, usecols = 12)

	fig21 = plt.figure(num='dUdV', figsize=(13, 7.2), dpi=100)
	ax21 = fig21.add_subplot(111)
	ax21.set_xlabel(r'$N_e$', fontsize=17)
	ax21.set_ylabel('', fontsize=17)
	ax21.set_xscale('linear')
	ax21.set_yscale('linear')
	ax21.plot(NN_ustartot, dU_10, color = "red", linestyle = '-', marker = '', label=r"$dU_{10}$")
	ax21.plot(NN_ustartot, dU_01, color = "orange", linestyle = '--', marker = '', label=r"$dU_{01}$")
	ax21.plot(NN_ustartot, dV_10, color = "blue", linestyle = '-', marker = '', label=r"$dV_{10}$")
	ax21.plot(NN_ustartot, dV_01, color = "green", linestyle = '--', marker = '', label=r"$dV_{01}$")
	ax21.legend(loc='best', prop={'size': 16})
	ax21.grid(False)
	#ax21.yaxis.set_minor_formatter(NullFormatter())

	fig22 = plt.figure(num='dUpdVp', figsize=(13, 7.2), dpi=100)
	ax22 = fig22.add_subplot(111)
	ax22.set_xlabel(r'$N_e$', fontsize=17)
	ax22.set_ylabel('', fontsize=17)
	ax22.set_xscale('linear')
	ax22.set_yscale('linear')
	ax22.plot(NN_ustartot, dUp_10, color = "red", linestyle = '-', marker = '', label=r"$dU'_{10}$")
	ax22.plot(NN_ustartot, dUp_01, color = "orange", linestyle = '--', marker = '', label=r"$dU'_{01}$")
	ax22.plot(NN_ustartot, dVp_10, color = "blue", linestyle = '-', marker = '', label=r"$dV'_{10}$")
	ax22.plot(NN_ustartot, dVp_01, color = "green", linestyle = '--', marker = '', label=r"$dV'_{01}$")
	ax22.legend(loc='best', prop={'size': 16})
	ax22.grid(False)
	#ax22.yaxis.set_minor_formatter(NullFormatter())

	fig23 = plt.figure(num='dAdS', figsize=(13, 7.2), dpi=100)
	ax23 = fig23.add_subplot(111)
	ax23.set_xlabel(r'$N_e$', fontsize=17)
	ax23.set_ylabel('', fontsize=17)
	ax23.set_xscale('linear')
	ax23.set_yscale('linear')
	ax23.plot(NN_ustartot, dA_10, color = "red", linestyle = '-', marker = '', label=r"$dA_{10}$")
	ax23.plot(NN_ustartot, dA_01, color = "orange", linestyle = '--', marker = '', label=r"$dA_{01}$")
	ax23.plot(NN_ustartot, dS_10, color = "blue", linestyle = '-', marker = '', label=r"$dS_{10}$")
	ax23.plot(NN_ustartot, dS_01, color = "green", linestyle = '--', marker = '', label=r"$dS_{01}$")
	ax23.legend(loc='best', prop={'size': 16})
	ax23.grid(False)
	#ax23.yaxis.set_minor_formatter(NullFormatter())

	fig24 = plt.figure(num='RS', figsize=(13, 7.2), dpi=100)
	ax24 = fig24.add_subplot(111)
	ax24.set_xlabel(r'$N_e$', fontsize=17)
	ax24.set_ylabel('', fontsize=17)
	ax24.set_xscale('linear')
	ax24.set_yscale('log')
	ax24.plot(NN_ustartot, np.abs(R_10), color = "red", linestyle = '-', marker = '', label=r"$R_{10}$")
	ax24.plot(NN_ustartot, np.abs(R_01), color = "orange", linestyle = '--', marker = '', label=r"$R_{01}$")
	ax24.plot(NN_ustartot, np.abs(S_10), color = "blue", linestyle = '-', marker = '', label=r"$S_{10}$")
	ax24.plot(NN_ustartot, np.abs(S_01), color = "green", linestyle = '--', marker = '', label=r"$S_{01}$")
	ax24.legend(loc='best', prop={'size': 16})
	ax24.grid(False)
	#ax24.yaxis.set_minor_formatter(NullFormatter())

	fig25 = plt.figure(num='T_U_V', figsize=(13, 7.2), dpi=100)
	ax25 = fig25.add_subplot(111)
	ax25.set_xlabel(r'$N_e$', fontsize=17)
	ax25.set_ylabel('', fontsize=17)
	ax25.set_xscale('linear')
	ax25.set_yscale('log')
	ax25.plot(NN_ustartot, np.abs(T_UU), color = "red", linestyle = '-', marker = '', label=r"$T_{UU}$")
	ax25.plot(NN_ustartot, np.abs(T_UV), color = "orange", linestyle = '--', marker = '', label=r"$T_{UV}$")
	ax25.plot(NN_ustartot, np.abs(T_VU), color = "blue", linestyle = '-', marker = '', label=r"$T_{VU}$")
	ax25.plot(NN_ustartot, np.abs(T_VV), color = "green", linestyle = '--', marker = '', label=r"$T_{VV}$")
	ax25.legend(loc='best', prop={'size': 16})
	ax25.grid(False)
	#ax25.yaxis.set_minor_formatter(NullFormatter())

	
	Ne, T_RRe = InterpAndSmooth(x=NN_ustartot, y=np.abs(T_RR), intkind='linear', Nsteps=1e3)
	Ne, T_RSe = InterpAndSmooth(x=NN_ustartot, y=np.abs(T_RS), intkind='linear', Nsteps=1e3)
	Ne, T_SRe = InterpAndSmooth(x=NN_ustartot, y=np.abs(T_SR), intkind='linear', Nsteps=1e3)
	Ne, T_SSe = InterpAndSmooth(x=NN_ustartot, y=np.abs(T_SS), intkind='linear', Nsteps=1e3)

	fig26 = plt.figure(num='T_R_S', figsize=(13, 7.2), dpi=100)
	ax26 = fig26.add_subplot(111)
	ax26.set_xlabel(r'$N_e$', fontsize=17)
	ax26.set_ylabel('', fontsize=17)
	ax26.set_xscale('linear')
	ax26.set_yscale('log')
	#ax26.set_xlim(11.0, 65.0)
#	ax26.plot(NN_ustartot, np.abs(T_RR), color = "red", linestyle = '-', marker = '', label=r"$T_{RR}$")
#	ax26.plot(NN_ustartot, np.abs(T_RS), color = "orange", linestyle = '--', marker = '', label=r"$T_{RS}$")
#	ax26.plot(NN_ustartot, np.abs(T_SR), color = "blue", linestyle = '-', marker = '', label=r"$T_{SR}$")
#	ax26.plot(NN_ustartot, np.abs(T_SS), color= "green", linestyle = '--', marker = '', label=r"$T_{SS}$")
	ax26.plot(Ne, T_RRe, color = "red", linestyle = '-', marker = '', label=r"$T_{RR}$")
	ax26.plot(Ne, T_RSe, color = "orange", linestyle = '--', marker = '', label=r"$T_{RS}$")
	ax26.plot(Ne, T_SRe, color = "blue", linestyle = '-', marker = '', label=r"$T_{SR}$")
	ax26.plot(Ne, T_SSe, color= "green", linestyle = '--', marker = '', label=r"$T_{SS}$")
	ax26.legend(loc='best', prop={'size': 16})
	ax26.grid(False)
	#ax26.yaxis.set_minor_formatter(NullFormatter())

	
	N_TRR_Jim = np.loadtxt(fname='T_RR_Jim.dat', delimiter=' ', usecols = 0)
	TRR_Jim = np.loadtxt(fname='T_RR_Jim.dat', delimiter=' ', usecols = 1)
	#
	N_TRS_Jim = np.loadtxt(fname='T_RS_Jim.dat', delimiter=' ', usecols = 0)
	TRS_Jim = np.loadtxt(fname='T_RS_Jim.dat', delimiter=' ', usecols = 1)
	#
	N_TSR_Jim = np.loadtxt(fname='T_SR_Jim.dat', delimiter=' ', usecols = 0)
	TSR_Jim = np.loadtxt(fname='T_SR_Jim.dat', delimiter=' ', usecols = 1)
	#
	N_TSS_Jim = np.loadtxt(fname='T_SS_Jim.dat', delimiter=' ', usecols = 0)
	TSS_Jim = np.loadtxt(fname='T_SS_Jim.dat', delimiter=' ', usecols = 1)

	fig27 = plt.figure(num='T_R_S_cmp', figsize=(13, 7.2), dpi=100)
	ax27 = fig27.add_subplot(111)
	ax27.set_xlabel(r'$N_e$', fontsize=17)
	ax27.set_ylabel('', fontsize=17)
	ax27.set_xscale('linear')
	ax27.set_yscale('log')
	ax27.set_xlim(20.0, 77.0)
	ax27.plot(NN_ustartot, np.abs(T_RR), color = "red", linestyle = '-', marker = '', label=r"$T_{RR}$ Matteo")
	ax27.plot(NN_ustartot, np.abs(T_RS), color = "orange", linestyle = '-', marker = '', label=r"$T_{RS}$ Matteo")
	ax27.plot(NN_ustartot, np.abs(T_SR), color = "blue", linestyle = '-', marker = '', label=r"$T_{SR}$ Matteo")
	ax27.plot(NN_ustartot, np.abs(T_SS), color= "green", linestyle = '-', marker = '', label=r"$T_{SS}$ Matteo")
	ax27.plot(N_TRR_Jim, np.abs(TRR_Jim), color = "magenta", linestyle = '--', marker = '', label=r"$T_{RR}$ Jim")
	ax27.plot(N_TRS_Jim, np.abs(TRS_Jim), color = "yellow", linestyle = '--', marker = '', label=r"$T_{RS}$ Jim")
	ax27.plot(N_TSR_Jim, np.abs(TSR_Jim), color = "cyan", linestyle = '--', marker = '', label=r"$T_{SR}$ Jim")
	ax27.plot(N_TSS_Jim, np.abs(TSS_Jim), color= "black", linestyle = '--', marker = '', label=r"$T_{SS}$ Jim")
	ax27.legend(loc='best', prop={'size': 16})
	ax27.grid(False)
	#ax27.yaxis.set_minor_formatter(NullFormatter())



	X, Y = np.meshgrid(x, y)
	VE = myTools.VE(x=X, y=Y, mllpxilpp=(m[i], l[i], lp[i], xi[i], lpp[i]))
	
	fig28 = plt.figure(num='VE', figsize=(13, 7.2), dpi=100)
	ax28 = fig28.gca(projection='3d')
	ax28.set_xlabel(r'$X$', fontsize=17)
	ax28.set_ylabel(r'$Y$', fontsize=17)
	ax28.set_zlabel(r'$V_E$', fontsize=17)
	ax28.set_xscale('linear')
	ax28.set_yscale('linear')
	ax28.set_zscale('linear')
	ax28.set_xlim(0, 25.0)
	ax28.set_ylim(0, 8.0)
	ax28.set_zlim(0, 6.0e-10)
	ax28.view_init(elev=30, azim=120)
	ax28.plot_surface(X, Y, VE, cmap='coolwarm', linewidth=0, antialiased=False)
	#ax28.legend(loc='best', prop={'size': 16})
	ax28.grid(False)
	ob = offsetbox.AnchoredText(s=text, loc='upper left', pad=0.3, borderpad=0.85, frameon=True, prop=dict(color='black', size=15))
	ob.patch.set(boxstyle='round', edgecolor='lightgray', alpha=1.0)
	ax28.add_artist(ob)
	ob1 = offsetbox.AnchoredText(s=textModel, loc='lower left', pad=0.3, borderpad=0.85, frameon=False, prop=dict(color='black', size=15))
	ax28.add_artist(ob1)
	#ax28.yaxis.set_minor_formatter(NullFormatter())

	Ve = myTools.VE(x=x, y=y, mllpxilpp=(m[i], l[i], lp[i], xi[i], lpp[i]))

	fig29 = plt.figure(num='Ve', figsize=(13, 7.2), dpi=100)
	ax29 = fig29.add_subplot(111)
	ax29.set_xlabel(r'$N_e$', fontsize=17)
	ax29.set_ylabel(r'$V_E$', fontsize=17)
	ax29.set_xscale('linear')
	ax29.set_yscale('linear')
	ax29.plot(NN, Ve, color = "green", linestyle = '-', marker = '')
	#ax29.legend(loc='best', prop={'size': 16})
	ax29.grid(False)
	#ax29.yaxis.set_minor_formatter(NullFormatter())


	Cem = 2.0 - np.log(2.0) - 0.5772
	As_small = [1.0-2.0*epsA[i]+2.0*Cem*(3.0*epsA[i]-etaAA[i]) for i in range(0, len(NN))]
	etaAS_func = interpolate.InterpolatedUnivariateSpline(x=NN, y=etaAS, k=3)
	etaAS_interp = etaAS_func(NN_ustartot)
	As_corr = [-4.0*Cem*etaAS_interp[i]*T_RS[i] for i in range(0, len(NN_ustartot))]
	Nend = 63.0 # Model 1

	fig30 = plt.figure(num='eq15', figsize=(13, 7.2), dpi=100)
	ax30 = fig30.add_subplot(111)
	ax30.set_xlabel(r'$N_e$', fontsize=17)
	ax30.set_ylabel('', fontsize=17)
	ax30.set_xscale('linear')
	ax30.set_yscale('log')
	ax30.set_xlim(NN_ustartot[0], Nend)
	ax30.plot(NN, As_small, color = "blue", linestyle = '-', marker = '', label=r"$1 - 2 \epsilon_A + 2C (3\epsilon_A - \eta_{AA})$")
	ax30.plot(NN_ustartot, As_corr, color = "red", linestyle = '-', marker = '', label=r"$-4C \eta_{AS} T_{RS}$")
	ax30.legend(loc='best', prop={'size': 16})
	ax30.grid(False)
	#ax30.yaxis.set_minor_formatter(NullFormatter())



#	fig0.show()
#	fig1.show()
#	fig2.show()
#	fig3.show()
#	fig4.show()
#	fig5.show()
#	fig6.show()
#	fig7.show()
#	fig8.show()
#	fig9.show()
#	fig10.show()
#	fig11.show()
#	fig12.show()
#	fig121.show()
#	fig13.show()
#	fig14.show()
#	fig15.show()
#	fig16.show()
#	fig17.show()
#	fig18.show()
#	fig19.show()
#	fig20.show()
#	fig21.show()
#	fig22.show()
#	fig23.show()
#	fig24.show()
#	fig25.show()
	fig26.show()
	fig27.show()
#	fig28.show()
#	fig29.show()
#	fig30.show()


	raw_input()
