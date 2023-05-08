import shutil
from decimal import *
getcontext().prec = 8 # memorize 8 digits after the comma

import numpy as np
from scipy import interpolate
import myTools

if __name__ == '__main__':

	param2D = myTools.Create2DListFromTXT(NameTXTinput='inf_param.dat')
	myTools.CreateNewDir(dirName='out', isLocal=True)
	shutil.copy("../inf_param.dat", "./")

	if len(param2D) > 0 :
		print "Number of lines in parameter file: ", len(param2D)
		ln = 0
		for line in param2D :
			str_a = "inf_" + str(ln) + ".dat"
			file_a = open(str_a, "w+")
			file_a.write("NN\t\t\tx\t\t\tpix\t\t\ty\t\t\tpiy\t\t\tpsi\n")
			str_b = "nbni_" + str(ln) + ".dat"
			file_b = open(str_b, "w+")
			file_b.write("NN\t\t\tH\t\t\tnb\t\t\tni\t\t\tlog10(nb/ni)\n")
			str_c = "theta_" + str(ln) + ".dat"
			file_c = open(str_c, "w+")
			file_c.write("NN\t\t\tUd\t\t\tVd\t\t\ttheta\t\t\tcosth\t\t\tsinth\t\t\tthetadot\t\t\tthetadot_alt\n")
			str_d = "UV_" + str(ln) + ".dat"
			file_d = open(str_d, "w+")
			file_d.write("NN\t\t\tVe\t\t\tVeU\t\t\tVeV\t\t\tepsU\t\t\tepsV\t\t\tetaUU\t\t\tetaUV\t\t\tetaVU\t\t\tetaVV\n")
			str_e = "gamma_" + str(ln) + ".dat"
			file_e = open(str_e, "w+")
			file_e.write("NN\t\t\ttan(gamma_Matteo)\t\t\tcos(gamma_Matteo)\t\t\tsin(gamma_Matteo)\t\t\ttan(gamma_Jim)\t\t\tcos(gamma_Jim)\t\t\tsin(gamma_Jim)\n")
			str_f = "AS_" + str(ln) + ".dat"
			file_f = open(str_f, "w+")
			file_f.write("NN\t\t\tAd\t\t\tVeA\t\t\tVeS_def\t\t\tVeS\t\t\tVeS_alt\t\t\tepsA\t\t\tepsS\t\t\tetaAA\t\t\tetaAS\t\t\tetaSS\n")

			str_g = "nsrOLD_" + str(ln) + ".dat"
			file_g = open(str_g, "w+")
			file_g.write("NN\t\t\tNk\t\t\tns\t\t\tr\n")
			str_h = "nsr_" + str(ln) + ".dat"
			file_h = open(str_h, "w+")
			file_h.write("NN\t\t\tNk\t\t\tns\t\t\tr\n")
			str_l = "nsrTOT_" + str(ln) + ".dat"
			file_l = open(str_l, "w+")
			file_l.write("NN\t\t\tNk\t\t\tns\t\t\tr\n")

			str_m = "Nend_" + str(ln) + ".dat"
			file_m = open(str_m, "w+")
			str_n = "Nstar_" + str(ln) + ".dat"
			file_n = open(str_n, "w+")

			### transfer functions ###
			str_o = "dAdS10_" + str(ln) + ".dat"
			file_o = open(str_o, "w+")
			file_o.write("NN_ustartot\t\t\tdU_10\t\t\tdUp_10\t\t\tdV_10\t\t\tdVp_10\t\t\tT_UU\t\t\tT_VU\t\t\tdA_10\t\t\tdS_10\t\t\tR_10\t\t\tS_10\t\t\tT_RR\t\t\tT_SR\n")
			str_p = "dAdS01_" + str(ln) + ".dat"
			file_p = open(str_p, "w+")
			file_p.write("NN_ustartot\t\t\tdU_01\t\t\tdUp_01\t\t\tdV_01\t\t\tdVp_01\t\t\tT_UV\t\t\tT_VV\t\t\tdA_01\t\t\tdS_01\t\t\tR_01\t\t\tS_01\t\t\tT_RS\t\t\tT_SS\n")

			str_q = "VeLkin_" + str(ln) + ".dat"
			file_q = open(str_q, "w+")
			file_q.write("NN\t\t\tVe\t\t\tLkin\n")

			m = line[0] # potential parameter
			l = line[1] # potential parameter
			lp = line[2] # potential parameter
			xi = line[3] # potential parameter
			lpp = line[4] # baryogenesis parameter
			x0 = line[5] # initial value of x
			pix0 = 0.0 # initial value of pix
			y0 = line[6] # initial value of y
			piy0 = 0.0
			psi0 = 0.0
			Nstart = 0.0 # initial number of e-folding (for odeint)
			Ntot = line[7] # final number of e-folding (for odeint)
			Nsteps = int(1000+1) # number of steps in NN list (its length; +1 is for nice purposes)
			IntMethod = 'RK45'
			#IntMethod = 'LSODA'

			mllpxilpp = (m, l, lp, xi, lpp)

			pix00, piy00 = myTools.Solve_pi(xpixypiy0=[x0, pix0, y0, piy0], mllpxilpp=mllpxilpp, eps=1e-10)
			
			print "\n[pix0, piy0] = [%.8e, %.8e]\n" %(Decimal(pix00), Decimal(piy00))

			xpixypiypsi0 = [x0, pix00, y0, piy00, psi0]
			NN0 = np.linspace(start=Nstart, stop=Ntot, num=Nsteps) # Nsteps = len(NN0)


			NN, x, pix, y, piy, psi, Nend1 = myTools.Solve_EoMs(NN0=NN0, xpixypiypsi0=xpixypiypsi0, mllpxilpp=mllpxilpp, method=IntMethod, rtol=1.e-4, atol=1.e-7, dense_output=False, ist_eval=True)
			# available integrator: RK45, RK23, Radau, BDF, LSODA (used)
			print "Nend1 = ", Nend1
			xpixypiy = [x, pix, y, piy]

			nb = myTools.NB(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
			ni = myTools.NI(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)
			H = myTools.Hubble(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)

			Lkin = myTools.Lkin(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp)

			Ud, Vd, theta, costh, sinth, thd, thd_alt = myTools.Compute_theta(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi)

			iend, Nend, xend, pixend, yend, piyend, psiend, rhoend = myTools.GetNend(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi)
			
			Ve, VeU, VeV, epsU, epsV, etaUU, etaUV, etaVU, etaVV = myTools.SR_params_UV(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)
			
			ns_old, r_old = myTools.nsrOLD(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi, Nend=Nend)

			tang_Matteo, cosg_Matteo, sing_Matteo, tang_Jim, cosg_Jim, sing_Jim = myTools.Compute_gamma(xpixypiy=xpixypiy, mllpxilpp=mllpxilpp, psi=psi)

			Ad, VeA, VeS_def, VeS, VeS_alt, epsA, epsS, etaAA, etaAS, etaSS = myTools.SR_params_AS(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi, altDefVeS=True)
			
			As, ns, r = myTools.nsrNEW(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi, altDefVeS=False)
			
			star_ns, star_r = myTools.GetNstar(xpixypiy=xpixypiy, psi=psi, NN=NN, mllpxilpp=mllpxilpp, As=As, Nsteps=Nsteps, eps=1.e-10)
			Nstar_ns, xstar_ns, pixstar_ns, ystar_ns, piystar_ns, psistar_ns, rescalestar_ns = star_ns
			Nstar_r, xstar_r, pixstar_r, ystar_r, piystar_r, psistar_r, rescalestar_r = star_r

			#Nshift = np.abs(Nstar_r - Nstar_ns)
			#ns_shift, r_shift = myTools.Adjust_nsr(NN=NN, ns=ns, r=r, Nshift=Nshift, Nsteps=Nsteps)


			### transfer functions ###
			# dU and dV #
			ustar_ns = Nend - Nstar_ns
			#dUdV_10, dUdV_01 = myTools.Solve_dUdV(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi, ustar_ns=ustar_ns, Nsteps=Nsteps, method=IntMethod, rtol=1.e-4, atol=1.e-7)
			#NN_ustartot_10, dU_10, dUp_10, dV_10, dVp_10, T_UU, T_VU = dUdV_10
			#NN_ustartot_01, dU_01, dUp_01, dV_01, dVp_01, T_UV, T_VV = dUdV_01

			# dA and dS #
			dAdS_10, dAdS_01 = myTools.Solve_dAdS(xpixypiy=xpixypiy, NN=NN, mllpxilpp=mllpxilpp, psi=psi, ustar_ns=ustar_ns, Nsteps=Nsteps, method=IntMethod, rtol=1.e-4, atol=1.e-7, dense_output=False, events=None, ist_eval=True)
			NN_ustartot_10, dU_10, dUp_10, dV_10, dVp_10, T_UU, T_VU, dA_10, dS_10, R_10, S_10, T_RR, T_SR = dAdS_10
			NN_ustartot_01, dU_01, dUp_01, dV_01, dVp_01, T_UV, T_VV, dA_01, dS_01, R_01, S_01, T_RS, T_SS = dAdS_01

			if np.array_equal(NN_ustartot_10, NN_ustartot_01) :
				NN_ustartot = NN_ustartot_10

			
			for i in range(0, len(NN)) :
				xpixypiy_i = [x[i], pix[i], y[i], piy[i]]

				nbni_i = nb[i]/ni[i] 
				if nbni_i != 0.0 :
					log10nbni_i = np.log10(np.abs(nbni_i))
				else :
					continue

				Nk = Nend - NN[i]

				file_a.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(Decimal(NN[i]), Decimal(x[i]), Decimal(pix[i]), Decimal(y[i]), Decimal(piy[i]), Decimal(psi[i])))
				file_b.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(Decimal(NN[i]), Decimal(H[i]), Decimal(nb[i]), Decimal(ni[i]), Decimal(log10nbni_i)))
				file_c.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(Decimal(NN[i]), Decimal(Ud[i]), Decimal(Vd[i]), Decimal(theta[i]), Decimal(costh[i]), Decimal(sinth[i]), Decimal(thd[i]), Decimal(thd_alt[i])))
				file_d.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(Decimal(NN[i]), Decimal(Ve[i]), Decimal(VeU[i]), Decimal(VeV[i]), Decimal(epsU[i]), Decimal(epsV[i]), Decimal(etaUU[i]), Decimal(etaUV[i]), Decimal(etaVU[i]), Decimal(etaVV[i])))
				file_e.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(Decimal(NN[i]), Decimal(tang_Matteo[i]), Decimal(cosg_Matteo[i]), Decimal(sing_Matteo[i]), Decimal(tang_Jim[i]), Decimal(cosg_Jim[i]), Decimal(sing_Jim[i])))
				file_f.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(Decimal(NN[i]), Decimal(Ad[i]), Decimal(VeA[i]), Decimal(VeS_def[i]), Decimal(VeS[i]), Decimal(VeS_alt[i]), Decimal(epsA[i]), Decimal(epsS[i]), Decimal(etaAA[i]), Decimal(etaAS[i]), Decimal(etaSS[i])))
				if Nk >= 50.0 and Nk <= 60.0 :
					# Planck plots contain ns and r computed at k* = 0.002 Mpc-1 --> ns, r without shift!
					file_g.write("%.8f\t%.8f\t%.8f\t%.8f\n" %(NN[i], Nk, ns_old[i], r_old[i]))
					file_h.write("%.8f\t%.8f\t%.8f\t%.8f\n" %(NN[i], Nk, ns[i], r[i])) # they were shifted
				file_l.write("%.8f\t%.8f\t%.8f\t%.8f\n" %(NN[i], Nk, ns[i], r[i]))

				file_q.write("%.8e\t%.8e\t%.8e\n" %(Decimal(NN[i]), Decimal(Ve[i]), Decimal(Lkin[i])))

			### transfer functions ###
			for i in range(0, len(NN_ustartot)) :	
				file_o.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(Decimal(NN_ustartot[i]), Decimal(dU_10[i]), Decimal(dUp_10[i]), Decimal(dV_10[i]), Decimal(dVp_10[i]), Decimal(T_UU[i]), Decimal(T_VU[i]), Decimal(dA_10[i]), Decimal(dS_10[i]), Decimal(R_10[i]), Decimal(S_10[i]), Decimal(T_RR[i]), Decimal(T_SR[i])))
				file_p.write("%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n" %(Decimal(NN_ustartot[i]), Decimal(dU_01[i]), Decimal(dUp_01[i]), Decimal(dV_01[i]), Decimal(dVp_01[i]), Decimal(T_UV[i]), Decimal(T_VV[i]), Decimal(dA_01[i]), Decimal(dS_01[i]), Decimal(R_01[i]), Decimal(S_01[i]), Decimal(T_RS[i]), Decimal(T_SS[i])))

			
			etaBend = myTools.etaB(xpixypiye=[xend, pixend, yend, piyend], mllpxilpp=mllpxilpp)
			ilate = iend + int(10 * Nsteps / (Ntot - Nstart)) # compute etaB after roughly 10 e-foldings after Nend
			etaBlate = myTools.etaB(xpixypiye=[x[ilate], pix[ilate], y[ilate], piy[ilate]], mllpxilpp=mllpxilpp)

			file_m.write("Nend = %f\nxend = %.8e\npixend = %.8e\nyend = %.8e\npiyend = %.8e\npsiend = %.8e\nrhoend = %.8e\netaBend = %.8e\netaBlate = %.8e\n" %(Nend, xend, pixend, yend, piyend, psiend, rhoend, etaBend, etaBlate))			
			print "\nNsteps = %d\nNend = %f\nxend = %f\npixend = %f\nyend = %f\npiyend = %f\npsiend = %f\nrhoend = %f\netaBend = %f\netaBlate = %f\n" %(Nsteps, Nend, xend, pixend, yend, piyend, psiend, rhoend, etaBend, etaBlate)

			Vstar_ns = myTools.VE(x=xstar_ns, y=ystar_ns, mllpxilpp=mllpxilpp)
			Vstar_r = myTools.VE(x=xstar_r, y=ystar_r, mllpxilpp=mllpxilpp)
			file_n.write("Nstar = %f\nxstar = %.8e\npixstar = %.8e\nystar = %.8e\npiystar = %.8e\npsistar = %.8e\nrescalestar = %.8e\nVstar = %.8e\n" %(Nstar_ns, xstar_ns, pixstar_ns, ystar_ns, piystar_ns, psistar_ns, rescalestar_ns, Vstar_ns))
			
			print "for ns: k* = 0.05 Mpc-1"
			print "Nstar = %f\nxstar = %f\npixstar = %f\nystar = %f\npiystar = %f\npsistar = %f\nrescalestar = %f\nVstar = %f\nVstar_rescale = %f\nrhoend = %f\nrhoend_rescale = %f\n" %(Nstar_ns, xstar_ns, pixstar_ns, ystar_ns, piystar_ns, psistar_ns, rescalestar_ns, Vstar_ns, rescalestar_ns*Vstar_ns, rhoend, rescalestar_ns*rhoend)
			print "for r: k* = 0.002 Mpc-1"
			print "Nstar = %f\nxstar = %f\npixstar = %f\nystar = %f\npiystar = %f\npsistar = %f\nrescalestar = %f\nVstar = %f\nVstar_rescale = %f\nrhoend = %f\nrhoend_rescale = %f\n" %(Nstar_r, xstar_r, pixstar_r, ystar_r, piystar_r, psistar_r, rescalestar_r, Vstar_r, rescalestar_r*Vstar_r, rhoend, rescalestar_r*rhoend)


			# Matteo: ns, r
			As_star_ns, ns_star_ns, r_star_ns = myTools.nsrNEW(xpixypiy=[xstar_ns, pixstar_ns, ystar_ns, piystar_ns], NN=Nstar_ns, mllpxilpp=mllpxilpp, psi=psistar_ns, altDefVeS=False)
			As_star_r, ns_star_r, r_star_r = myTools.nsrNEW(xpixypiy=[xstar_r, pixstar_r, ystar_r, piystar_r], NN=Nstar_r, mllpxilpp=mllpxilpp, psi=psistar_r, altDefVeS=False)

			file_n.write("As_star_ns = %.8e\nns_star_ns = %.8e\nr_star_ns = %.8e\nAs_star_r = %.8e\nns_star_r = %.8e\nr_star_r = %.8e\n" %(As_star_ns, ns_star_ns, r_star_ns, As_star_r, ns_star_r, r_star_r))
			print "As_star = %f\nns_star = %f\nr_star = %f\n" %(As_star_ns, ns_star_ns, r_star_r)


			file_a.close()
			file_b.close()
			file_c.close()
			file_d.close()
			file_e.close()
			file_f.close()
			file_g.close()
			file_h.close()
			file_l.close()
			file_m.close()
			file_n.close()
			file_o.close()
			file_p.close()
			file_q.close()
			print "Done line %d!" % ln
			ln += 1
