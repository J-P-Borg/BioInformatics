#
#	  MRA3_BioInfo_2_Pub.R
#     Copyright (C) 2022, 2023  Jean-Pierre BORG
#
#	  R software used to write the article to publish in BioInformatics. 
#	  The files to download concerning figures 4, 5 and SI, along with their URL, are indicated in the code.
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
#     Jean-Pierre BORG
#     University of Montpellier, France
#     Institut de Recherche en Cancérologie de Montpellier, INSERM, France
#	  Cancer Bioinformatics and Systems Biology
#     jean-pierre.borg@inserm.fr
#
#
#
#	README FIRST : to use this software :
#		- execute the instructions labeled "Initialization" (between two rows of #########################)
#		- IF you want to draw figure 1, 2, 3 or Supplementary Table 1, execute the instructions labeled "Joint statements for figures 1, 2, 3 and "Supplementary Table 1""
#		- Then execute the code corresponding to the figure wanted.
#
#########################################################################################################
#
#	Initialization
#
#	Reset of the R session and environment variables, to start again from scratch.
#
rm(list = ls())		# We have also to restart a new session to unload all the packages.
					# We can also detach every package individually by : detach(package:packagename)
Verbose	<- FALSE
					
#
# 	Declaration of the libraries and of the functions used in this document
#
library ("rootSolve")		# 	To use multiroot
library ("deSolve")			# 	To use ode
library ("Deriv")			#	To use Deriv
library ("ggplot2")			# 	To use ggplot
library ("glmnet")			#	To use glmnet
library ("dplyr")			# 	To use rename
library ("data.table")		# 	To use fread
library	("forcats")			#	To use fct_relevel
library ("pracma")			# 	To use trapz
library ("verification")	# 	To use roc.area
library ("minet")			# 	To use build.mim, clr, aracne, mrnet -- If necessary, install pkg "BiocManager", then BiocManager::install("minet")
library ("parallel")		#	To use parallel
library ("foreach")			#	To use foreach
library ("doParallel")		#	To use doParallel

#
# 	Declaration of the root of the filesystem for this document
#
vRoot		<- c("/home/jpborg/Documents/Publications/BioInformatics/MRA_Regression_ed1/Vers_Rev1/Resultats/")
vRootDC4	<- c("/home/jpborg/Documents/Publications/BioInformatics/MRA_Regression_ed1//Vers_Init/Donnees/DC4/")

#
#	Functions used by this program
#
#	Figures 1, 2 and 3
#
#' F
#'
#' @param 	X 		 A vector of six values representing the concentrations in the MAPK cascade model (see "Details")
#'
#' @details			Rate Expressions Used for the Models of MAPK Cascade (Kinase network , figures 1 to 3)
#' @details			These data come from the article : "Inferring dynamic architecture of cellular networks using time series of gene expression, protein and metabolite data", 
#' @details			(Supplementary Table 2 of their article. Rate expressions, differential equations and parameter values for the MAPK cascade model
#' @details			by Eduardo Sontag et al)
#'
#' @return 	F  		The values of the rate equations described in the article above
#'
F <- function(X)
   c(F1 = (KC1*U*X[1]) 	/ ((K11+X[1]+(MKKK0-X[1]-X[2])*K11/K12)*(1+X[6]/Ki)) - ((V4*(MKKK0-X[1]-X[2]))  / (K31+X[2]+(MKKK0-X[1]-X[2])*K31/K32+X[1]*K31/K33)),
	 F2 = (KC2*U*(MKKK0-X[1]-X[2]))    / ((K11+X[1]+(MKKK0-X[1]-X[2])*K11/K12)*(1+X[6]/Ki)) - ((V3*X[2])  / (K31+X[2]+(MKKK0-X[1]-X[2])*K31/K32+X[1]*K31/K33)),
	 F3 = (KC5*X[3]*X[2]) / (K51+X[3]+(MKK0-X[3]-X[4])*K51/K52) - ((V8*(MKK0-X[3]-X[4])*(1+A*X[6]/Kmp))   / ((K71+X[4]+(MKK0-X[3]-X[4])*K71/K72+X[3]*K71/K73)*(1+X[6]/Kmp))),
	 F4 = (KC6*(MKK0-X[3]-X[4])*X[2])  / (K51+X[3]+(MKK0-X[3]-X[4])*K51/K52) - ((V7*X[4]*(1+A*X[6]/Kmp))  / ((K71+X[4]+(MKK0-X[3]-X[4])*K71/K72+X[3]*K71/K73)*(1+X[6]/Kmp))),
	 F5 = (KC9*X[4]*X[5]) / (K91+X[5]+(MAPK0-X[5]-X[6])*K91/K92) - ((V12*(MAPK0-X[5]-X[6])) / (K111+X[6]+(MAPK0-X[5]-X[6])*K111/K112+X[5]*K111/K113)),
	 F6 = (KC10*X[4]*(MAPK0-X[5]-X[6]))/ (K91+X[5]*(MAPK0-X[5]-X[6])*K91/K92) - ((V11*X[6]) / (K111+X[6]+(MAPK0-X[5]-X[6])*K111/K112+X[5]*K111/K113)))


#' FF
#'
#' @param 	X1..X6 	Six numbers corresponding to the vector X[] above (function F)
#'
#' @details			Same as above
#'
#' @return 	FF 		The values of the rate equations described in the article above
#'
FF <- function(X1,X2,X3,X4,X5,X6)
   c(F1 = (KC1*U*X1) 	/ ((K11+X1+(MKKK0-X1-X2)*K11/K12)*(1+X6/Ki)) - ((V4*(MKKK0-X1-X2))  / (K31+X2+(MKKK0-X1-X2)*K31/K32+X1*K31/K33)),
	 F2 = (KC2*U*(MKKK0-X1-X2))    / ((K11+X1+(MKKK0-X1-X2)*K11/K12)*(1+X6/Ki)) - ((V3*X2)  / (K31+X2+(MKKK0-X1-X2)*K31/K32+X1*K31/K33)),
	 F3 = (KC5*X3*X2) / (K51+X3+(MKK0-X3-X4)*K51/K52) - ((V8*(MKK0-X3-X4)*(1+A*X6/Kmp))   / ((K71+X4+(MKK0-X3-X4)*K71/K72+X3*K71/K73)*(1+X6/Kmp))),
	 F4 = (KC6*(MKK0-X3-X4)*X2)  / (K51+X3+(MKK0-X3-X4)*K51/K52) - ((V7*X4*(1+A*X6/Kmp))  / ((K71+X4+(MKK0-X3-X4)*K71/K72+X3*K71/K73)*(1+X6/Kmp))),
	 F5 = (KC9*X4*X5) / (K91+X5+(MAPK0-X5-X6)*K91/K92) - ((V12*(MAPK0-X5-X6)) / (K111+X6+(MAPK0-X5-X6)*K111/K112+X5*K111/K113)),
	 F6 = (KC10*X4*(MAPK0-X5-X6))/ (K91+X5*(MAPK0-X5-X6)*K91/K92) - ((V11*X6) / (K111+X6+(MAPK0-X5-X6)*K111/K112+X5*K111/K113)))


#' RevSigNet
#'
#' @param 	t		Time (s)
#' @param 	X 		A vector of six values representing the concentrations in the MAPK cascade model (see function F above)
#' @param 	Parms	Parameters, if necessary
#'
#' @return 	RevSigNet  A list with the derivatives dX1/dt ... dX6/dt
#'

RevSigNet <- function(t, X, Parms) {
	with(as.list(c(X, Parms)), {
		y1   <- (KC1*U*X1) 	/ ((K11+X1+(MKKK0-X1-X2)*K11/K12)*(1+X6/Ki))
		y2   <- (KC2*U*(MKKK0-X1-X2)) / ((K11+X1+(MKKK0-X1-X2)*K11/K12)*(1+X6/Ki))
		y3   <- (V3*X2) 	/ (K31+X2+(MKKK0-X1-X2)*K31/K32+X1*K31/K33)
		y4   <- (V4*(MKKK0-X1-X2)) 	  / (K31+X2+(MKKK0-X1-X2)*K31/K32+X1*K31/K33)
		y5   <- (KC5*X3*X2) / (K51+X3+(MKK0-X3-X4)*K51/K52)
		y6   <- (KC6*(MKK0-X3-X4)*X2) / (K51+X3+(MKK0-X3-X4)*K51/K52)
		y7   <- (V7*X4*(1+A*X6/Kmp))  / ((K71+X4+(MKK0-X3-X4)*K71/K72+X3*K71/K73)*(1+X6/Kmp))
		y8   <- (V8*(MKK0-X3-X4)*(1+A*X6/Kmp)) / ((K71+X4+(MKK0-X3-X4)*K71/K72+X3*K71/K73)*(1+X6/Kmp))
		y9   <- (KC9*X4*X5) / (K91+X5+(MAPK0-X5-X6)*K91/K92)
		y10  <- (KC10*X4*(MAPK0-X5-X6))		   / (K91+X5*(MAPK0-X5-X6)*K91/K92)
		y11  <- (V11*X6) 	/ (K111+X6+(MAPK0-X5-X6)*K111/K112+X5*K111/K113)
		y12  <- (V12*(MAPK0-X5-X6))   / (K111+X6+(MAPK0-X5-X6)*K111/K112+X5*K111/K113)		
		
		dX1dt	<- y4-y1
		dX2dt	<- y2-y3
		dX3dt	<- y8-y5
		dX4dt	<- y6-y7
		dX5dt	<- y12-y9
		dX6dt	<- y10-y11
		
		return(list(c(dX1dt, dX2dt, dX3dt, dX4dt, dX5dt, dX6dt)))
	})
}		# RevSigNet


#
#	Figures 4, 5 and SI figures 1, 2 and 3
#
#' Digitalize
#'
#' @param 	MatE 		 Input matrix
#' @param 	Method		 Method used to digitalize ("Threshold", "IThreshold", "Top")
#' @param 	P1			 Parameter used by this method
#'
#' @details	"Threshold"	 Returns 1 iff abs(MatE[i,j] >= P1*max(abs(MatE))
#' @details	"IThreshold" Returns 1 iff MatE[i,j]     <  P1			(MatE must be >=0, if not, an error is returned)
#' @details	"Top"		 Sets a threshold Th such as Prob (abs(MatE[i,j] >= Th) ~ P1
#'
#' @return 	MatS  		 Output matrix filled with 0 and 1. Same size than MatE
#'
#'
Digitalize <- function(MatE, Method, P1) {
	Mat			<- abs(MatE)
	diag(Mat)	<- 0

	if (Method == "Threshold") {
		MatMax	<- max(Mat)
		MatS	<- ifelse(Mat > MatMax*P1, 1, 0)
	}

	if (Method == "IThreshold") {
		if (length(which(MatE < 0)) > 0) {
			print("ERROR : MatE must be >= 0 !")
			return(FALSE)
		}
		
		MatS	<- ifelse(Mat < P1, 1, 0)
	}
	
	if (Method == "Top") {
		Th 		<- quantile(Mat, prob = 1-P1)
		Th		<- max(0.005, Th)						# in case Th = 0
		MatS	<- ifelse(Mat >= Th, 1, 0)	
	}
	
	return(MatS)
}		# Digitalize


#
#' Digitalize2
#'
#' @param 	MatE 		 Input matrix
#' @param 	Method		 Method used to digitalize ("Threshold", "IThreshold", "Top")
#' @param 	P1			 Parameter used by this method
#'
#' @details				 Idem Digitalize, but returns 0 or the value of MatE instead of 1
#'
#' @return 	MatS  		 Output matrix filled with 0 and MatE. Same size than MatE
#'
#'
Digitalize2 <- function(MatE, Method, P1) {
	Mat			<- abs(MatE)
	diag(Mat)	<- 0

	if (Method == "Threshold") {
		MatMax	<- max(Mat)
		MatS	<- ifelse(Mat > MatMax*P1, MatE, 0)
	}

	if (Method == "IThreshold") {
		if (length(which(MatE < 0)) > 0) {
			print("ERROR : MatE must be >= 0 !")
			return(FALSE)
		}
		
		MatS	<- ifelse(Mat < P1, MatE, 0)
	}
	
	if (Method == "Top") {
		Th 		<- quantile(Mat, prob = 1-P1)
		Th		<- max(0.005, Th)						# in case Th = 0
		MatS	<- ifelse(Mat >= Th, MatE, 0)
	}
	
	return(MatS)
}		# Digitalize2


#' Benchmark
#'
#' @param 	Calc 		Computed matrix (0,1), to compare with reference
#' @param 	Ref			Reference matrix (0,1)
#'
#' @details				Calc and Ref must have the same size and contain 0 or 1 only.
#'
#' @return 	Benchmark	A vector : TP, TN, FP, FN, Sensitivity, Specificity
#'
#'
Benchmark <- function(Calc, Ref) {
	if ((dim(Calc)[1] != dim(Ref)[1]) || (dim(Calc)[2] != dim(Ref)[2])) {
		print("ERROR : Calc and Ref have different size !")
		return(FALSE)
	}
	
	diag(Calc)	<- 0														# Diagonal set to 0
	diag(Ref)	<- 0
	
	if (!all(Calc %in% c(0:1)) || !all(Ref %in% c(0:1))) {
		print("ERROR : Calc and Ref must contain 0 or 1 only !")
		return(FALSE)
	}
						
	P 	<- sum(Ref)															# True nbr of 1
	N 	<- dim(Calc)[1] * (dim(Calc)[2]-1) - P								# True nbr of 0
	
	TP	<- round(sum(ifelse(Calc & Ref, 1,0)), digits=0)					# Nbr of true 1 discovered
	TN	<- round(sum(ifelse(!Calc & !Ref, 1,0)) - dim(Calc)[1], digits=0)	# Nbr of true 0 discovered
	FP	<- round(N - TN, digits=0)											# Nbr of false 1 discovered
	FN	<- round(P - TP	, digits=0)											# Nbr of false 0 discovered
	Se 	<- round(TP/P, digits=3)											# Sensibility or sensitivity
	Sp	<- round(TN/N, digits=3)											# Specificity
	#	cat ("P",P," N",N," TP",TP," TN",TN," FP",FP," FN",FN," Se",Se," Sp",Sp, "\n")	
	return (c(TP, TN, FP, FN, Se, Sp))
}		# Benchmark


#' Benchmark2
#'
#' @param 	Calc 		Computed matrix (-1,0,1), to compare with reference
#' @param 	Ref			Reference matrix (-1,0,1)
#'
#' @details				Calc and Ref must have the same size and contain -1, 0 or 1 only.
#'
#' @return 	Benchmark2	A vector : TPS, TMS, T0, FPS, FMS, F0, P, N, Sensitivity, Specificity
#
#' @details : 	TPS (T+ : both Calc & Ref >0)
#' @details : 	TMS (T- : both Calc & Ref <0) 
#' @details : 	T0 (both Calc & Ref =0) , Specificity
#' @details : 	P (nbr of Ref != 0)
#' @details : 	N (nbr of Ref = 0)
#' @details : 	Sensitivity = (T+ + T-) / P
#' @details : 	T0 / N
#' @details : 	The terms in the diagonal are not taken into account
#'
#'
Benchmark2 <- function(Calc, Ref) {
	if ((dim(Calc)[1] != dim(Ref)[1]) || (dim(Calc)[2] != dim(Ref)[2])) {
		print("ERROR : Calc and Ref have different size !")
		return(FALSE)
	}
	
	diag(Calc)	<- 0															# Diagonal set to 0
	diag(Ref)	<- 0
							
	if (!all(Calc %in% c(-1, 0, 1)) || !all(Ref %in% c(-1, 0, 1))) {
		print("ERROR : Calc and Ref must contain -1, 0 or 1 only !")
		return(FALSE)
	}
	
	P 		<- sum(abs(Ref))													# True nbr of elements != 0
	N 		<- dim(Calc)[1] * (dim(Calc)[2]-1) - P								# True nbr of 0
	
	TPS		<- round(sum(ifelse((Calc == 1) & (Ref == 1), 1,0)), digits=0)		# Nbr of true +1 discovered
	TMS		<- round(sum(ifelse((Calc ==-1) & (Ref ==-1), 1,0)), digits=0)		# Nbr of true -1 discovered	
	T0		<- round(sum(ifelse((Calc == 0) & (Ref == 0), 1,0)), digits=0) - NBN	# Nbr of true 0 discovered
	
	FPS		<- round(sum(ifelse((Calc == 1) & (Ref != 1), 1,0)), digits=0)		# Nbr of false +1 discovered
	FMS		<- round(sum(ifelse((Calc ==-1) & (Ref !=-1), 1,0)), digits=0)		# Nbr of false -1 discovered	
	F0		<- round(sum(ifelse((Calc == 0) & (Ref != 0), 1,0)), digits=0)		# Nbr of false 0 discovered
	
	Se 	<- round((TPS+TMS)/P, digits=3)											# Sensibility or sensitivity
	Sp	<- round(T0/N, digits=3)												# Specificity
	#	cat ("P",P," N",N," TPS",TPS," TMS",TMS," T0",T0," Se",Se," Sp",Sp, "\n")	
	return (c(TPS, TMS, T0, FPS, FMS, F0, P, N, Se, Sp))
}		# Benchmark2

#	End of intialization
#
####################################################################################################################################
#
#	Joint statements for figures 1, 2, 3 and "Supplementary Table 1"
#

#	Calculation initialization
#
st <- c(200,  0, 100,  50, 100,  0)			# Start time

# 	Constants in relation with the kinetics of the reactions
KC1	 <- 1			# Catalytic rates in s-1
KC2	 <- 15
KC5	 <- 1
KC6	 <- 15
KC9	 <- 1
KC10 <- 15

K11	 <- 300			# Michaelis constants in nM
K12	 <- 20
K31	 <- 22
K32	 <- 18
K33	 <- 80
K51	 <- 300
K52	 <- 20
K71	 <- 22
K72	 <- 18
K73	 <- 80
K91	 <- 300
K92	 <- 20
K111 <- 22
K112 <- 18
K113 <- 80
Ki	 <- 100
Kmp	 <- 100

A	 <- 5			# Coefficient without dimension

V3	 <- 18.8		# Maximal enzyme rates in nM.s-1
V4	 <- 16.4
V7	 <- 18.8
V8	 <- 16.4
V11	 <- 8.4
V12	 <- 7.3

MKKK0	<- 200		# X1+X7+X2	: total concentration of the protein MKKK
MKK0	<- 180		# X3+X8+X4	: total concentration of the protein MKK
MAPK0	<- 360		# X5+X9+X6	: total concentration of the protein MAPK
U	<- 20			# Injection of smallGTPase Ras from time 0+

# 1/ Research of concentrations at steady state (dXi/dt = 0)

ss 		<- multiroot (f=F, start= st)	# uses rootSolve package
XStat 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])		# Theoritical concentrations at steady state
XMesNPMoy	<- 	mean(XStat)			# Noise standard deviation will be defined from this mean level
									# Measured concentrations in the 6 nodes, non perturbed network, steady state									

# 2/ To check that the steady state is reached at t = 1000 s

Ini	<- c(X1	 = 199.3489,
		 X2	 = 0.0334,
		 X3	 = 179.8226,
		 X4	 = 0.0091,
		 X5	 = 358.2074,
		 X6	 = 0.1519)
		 
#	Long term behavior : looking for steady state
Times	<- seq(0, 1001, by=1)		# time in s.
Parms	<- 0						# No perturbation

out	<-	ode(Ini, Times, RevSigNet, Parms)
out[950:1001,]						# Check that the six values X1..X6 remain constant

# 3/ Calculus of fi's partial derivatives, for t = 1000s 
 
Deriv(FF)
	 
X1 <- XStat[1]
X2 <- XStat[2]
X3 <- XStat[3]
X4 <- XStat[4]
X5 <- XStat[5]
X6 <- XStat[6]
	 
.e1  <- 	MKK0 - (X3 + X4)
.e2  <- 	MKKK0 - (X1 + X2)
.e3  <- 	MAPK0 - (X5 + X6)
.e4  <- 	1 + X6/Kmp
.e8  <- 	1 + X6/Ki
.e10 <- 	K71 * (.e1/K72 + 1 + X3/K73) + X4
.e11 <- 	.e4 * .e10
.e15 <- 	K11 * (.e2/K12 + 1) + X1
.e16 <- 	.e8 * .e15
.e20 <- 	K51 * (.e1/K52 + 1) + X3
.e28 <- 	K111 * (.e3/K112 + 1 + X5/K113) + X6
.e30 <- 	K31 * (.e2/K32 + 1 + X1/K33) + X2
.e31 <- 	.e16^2
.e32 <- 	.e11^2
.e33 <- 	1 + A * X6/Kmp
.e35 <- 	K91 * (1 + X5 * .e3/K92)
.e39 <- 	K91 * (.e3/K92 + 1) + X5
.e40 <- 	1/.e16
.e41 <- 	1/.e11
.e42 <- 	KC2 * U
.e46 <- 	1 - K11/K12
.e47 <- 	1 - K111/K112
.e48 <- 	1 - K31/K32
.e49 <- 	1 - K51/K52
.e50 <- 	1 - K71/K72
.e51 <- 	1/.e35
.e54 <- 	1/K113 - 1/K112
.e57 <- 	1/K33 - 1/K32
.e60 <- 	1/K73 - 1/K72
.e62 <- 	A/.e11 - .e33 * .e10/.e32
.e63 <- 	K12 * .e31
.e64 <- 	K92 * .e35^2
.e65 <- 	KC1 * U
.e66 <- 	KC10 * X4
.e67 <- 	KC6 * X2
.e68 <- 	Ki * .e31
	
dF1dX1	= .e65 * (.e40 - X1 * .e46 * .e8/.e31) + V4 * (1 + K31 * .e57 * .e2/.e30)/.e30
dF2dX1	= K31 * V3 * X2 *  .e57/.e30^2 - .e42 * (.e46 * .e8 * .e2/.e31 + .e40)
dF3dX1	= 0
dF4dX1	= 0
dF5dX1	= 0
dF6dX1	= 0

dF1dX2	=  K11 * KC1 * U * X1 * .e8/.e63 + V4 * (.e48 * .e2/.e30 + 1)/.e30
dF2dX2	=  .e42 * (K11 * .e8 * .e2/.e63 - .e40) - V3 * (1 - X2 * .e48/.e30)/.e30
dF3dX2	=  KC5 * X3/.e20
dF4dX2	=  KC6 * .e1/.e20
dF5dX2	=  0
dF6dX2	=  0

dF1dX3	=  0
dF2dX3	=  0
dF3dX3	= KC5 * X2 * (1 - X3 * .e49/.e20)/.e20 + V8 * .e33 * (.e41 + K71 * .e4 * .e60 * .e1/.e32)
dF4dX3	=  K71 * V7 * X4 * .e33 * .e4 * .e60/.e32 - .e67 * (.e49 * .e1/.e20 + 1) /.e20
dF5dX3	=  0
dF6dX3	=  0

dF1dX4	=  0
dF2dX4	=  0
dF3dX4	=  K51 * KC5 * X2 * X3/(K52 * .e20^2) + V8 * (.e50 * .e4 * .e1/.e32 + .e41) * .e33
dF4dX4	=  .e67 * (K51 * .e1/(K52 * .e20) - 1)/.e20 - V7 * .e33 * (.e41 - X4 * .e50 * .e4/.e32)
dF5dX4	=  KC9 * X5/.e39
dF6dX4	=  KC10 * .e3/.e35

dF1dX5	=  0
dF2dX5	=  0
dF3dX5	=  0
dF4dX5	=  0
dF5dX5	=  KC9 * X4 * (1 - X5 * (1 - K91/K92)/.e39)/.e39 + V12 * (1+ K111* .e54 * .e3/.e28)/.e28
dF6dX5	= K111* V11*X6 * .e54/.e28^2 - .e66 * (.e51 + K91 *(MAPK0 - (2*X5 + X6))* .e3/.e64)

dF1dX6	= -(.e65 * X1 * .e15/.e68)
dF2dX6	= -(.e42 * .e15 * .e2/.e68)
dF3dX6	= -(V8 * .e62 * .e1/Kmp)
dF4dX6	= -(V7 * X4 * .e62/Kmp)
dF5dX6	=  K91 *  KC9 * X4 * X5/(K92 * .e39^2) + V12 * (.e47 * .e3/.e28 +  1)/.e28
dF6dX6	= .e66 * (K91 * X5 * .e3/.e64 - .e51) - V11 * (1 - X6 * .e47/.e28)/.e28

#	Computation of the theoretical "r" matrix

NBN				<-	6								# Number of nodes
MatrTh			<- matrix(, nrow=NBN, ncol=NBN)		# Theoretical "r" matrix: exact value
	
MatrTh[1,]		<- -c(X1*dF1dX1, X2*dF1dX2, X3*dF1dX3, X4*dF1dX4, X5*dF1dX5, X6*dF1dX6) / (X1*dF1dX1)
MatrTh[2,]		<- -c(X1*dF2dX1, X2*dF2dX2, X3*dF2dX3, X4*dF2dX4, X5*dF2dX5, X6*dF2dX6) / (X2*dF2dX2)
MatrTh[3,]		<- -c(X1*dF3dX1, X2*dF3dX2, X3*dF3dX3, X4*dF3dX4, X5*dF3dX5, X6*dF3dX6) / (X3*dF3dX3)
MatrTh[4,]		<- -c(X1*dF4dX1, X2*dF4dX2, X3*dF4dX3, X4*dF4dX4, X5*dF4dX5, X6*dF4dX6) / (X4*dF4dX4)
MatrTh[5,]		<- -c(X1*dF5dX1, X2*dF5dX2, X3*dF5dX3, X4*dF5dX4, X5*dF5dX5, X6*dF5dX6) / (X5*dF5dX5)
MatrTh[6,]		<- -c(X1*dF6dX1, X2*dF6dX2, X3*dF6dX3, X4*dF6dX4, X5*dF6dX5, X6*dF6dX6) / (X6*dF6dX6)

#########################################################################################################################

#
#		Figure 1
#		6 nodes MAP kinase network
#

#	Drawing of the theoretical "r" matrix

NBN				<-	6								# Number of nodes
MatrTh2			<- matrix(, nrow=NBN, ncol=NBN)		# same as abs(MatrTh), but the diagonal is set to 0 (no auto-interaction in this network)
	
MatrTh2			<- abs(MatrTh)
diag(MatrTh2)	<- 0

nbRows			<- sum(ifelse(MatrTh2>0, 1, 0))		# Number of edges of this network
DessMAP			<- matrix(nrow=nbRows, ncol=4)		# Matrix csv used by Cytoscape to draw the network
nodeNames		<- c("MKKK(1)", "MKKK-PP(2)", "MKK(3)", "MKK-PP(4)", "MAPK(5)", "MAPK-PP(6)")	# Name of the nodes = their number (default)
ColNames		<- c("Source", "Target", "Interaction", "Value")		# The first line corresponds to the names
													# Source,target : name of the corresponding nodes
													# Interaction	: "d" (direct), "r" (reverse), "fp" (false positive), "fn" (false negative)
													# Value : value associated with the edge
colnames(DessMAP)	<- ColNames

iSol			<- 0								# Index of the row added in DessMAP
for (iRow in 1:NBN) {
	Col <- 	1:NBN
	Col <-	Col[-iRow]
	for (iCol in Col) {
		if (MatrTh[iRow, iCol] != 0) {
			iSol	<- iSol + 1
			DessMAP[iSol, "Source"]			<- nodeNames[iCol]
			DessMAP[iSol, "Target"]			<- nodeNames[iRow]			
			DessMAP[iSol, "Value"]			<- round(MatrTh[iRow, iCol], digits=3)
		}	
		if (MatrTh[iRow, iCol] > 0) {
			DessMAP[iSol, "Interaction"]	<- "d"
		}
		if (MatrTh[iRow, iCol] < 0) {
			DessMAP[iSol, "Interaction"]	<- "r"
		}		
	}	# iCol
}		# iRow	

write.csv(DessMAP, paste(vRoot, "MRA3.csv", sep=""), quote=FALSE)	
														# File used by Cytoscape v 3.9.1 to draw the figure 1a (import network from file)
														# Use style "MRA_styles.xml"
														# The picture is exported as "Fig1a_MRA3.pdf" and converted to "Fig1a_MRA3.png" (600 dpi) by Inkscape v 1.2.
														
#	Fig1b = MatrTh
#	Fig1c	: see Figure 2, MatrCc[,,3]
#	Fig1d	: see Figure 2, MatrCc[,,11]
#	These values are put in a document Word ("Fig1b.docx" etc...), then exported in "Fig1b.pdf" etc...														


#
#		Figure 2
#		Impact of perturbation and noise level for the 6 nodes MAP kinase network
#

NBN			<-	6										# Number of nodes
Perturbs	<- c(0.2, 0.3, 0.5, 0.7, 0.9, 0.99, 1.01, 1.10, 1.2, 1.3, 1.50, 1.6, 1.7)		# 13 perturbations
NBP			<- length(Perturbs)									
			
XMesNP		<- array(dim=c(NBN))						# Abundance of the tested product (protein, gene ...) WITHOUT perturbation
XMesP		<- array(dim=c(NBN, NBN))					# Abundance of the tested product (protein, gene ...) WITH perturbation
MatR		<- array(dim=c(NBN, NBN))					# Global Response matrix ("R")
pX   		<- array(dim=c(NBN-1, NBN-1))				# row (iRow), column (iCol)	
pY   		<- array(dim=c(NBN-1))						# row
	
MatrCc		<- array(-1, dim=c(NBN, NBN, NBP))			# Computed Connection matrix ("r")  -  we keep the NBP values for tests afterward.
ErrorQ		<- vector(length=NBP)						# Quadratic error

#		Fig. 2a : no noise is added (impact of perturbation level only)

ss <- multiroot (f=F, start= st)						# Non perturbed measure
XMesNP[ ]	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])

for (iPert in c(1:NBP)) {
	cat("iPert ", iPert," st ", st, "\n")
	Perturb	 <- Perturbs[iPert]							# Perturbation level
	
	V4 <- Perturb*V4		# P1
	ss <- multiroot (f=F, start= st)
	XMesP[ ,1] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V4 <- V4/Perturb		# Return to initial value
	
	V3 <- Perturb*V3		# P2
	ss <- multiroot (f=F, start= st)
	XMesP[ ,2] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V3 <- V3/Perturb		# Return to initial value
	
	V8 <- Perturb*V8		# P3
	ss <- multiroot (f=F, start= st)
	XMesP[ ,3] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V8 <- V8/Perturb		# Return to initial value
	
	V7 <- Perturb*V7		# P4
	ss <- multiroot (f=F, start= st)
	XMesP[ ,4] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V7 <- V7/Perturb		# Return to initial value
	
	V12 <- Perturb*V12		# P5
	ss <- multiroot (f=F, start= st)
	XMesP[ ,5] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V12 <- V12/Perturb		# Return to initial value
	
	V11 <- Perturb*V11		# P6
	ss <- multiroot (f=F, start= st)
	XMesP[ ,6] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V11 <- V11/Perturb		# Return to initial value

	for (iCol in 1:NBN) {
		MatR[ ,iCol]  <- 2 * (XMesP[ ,iCol]-XMesNP[ ]) / (XMesP[ ,iCol]+XMesNP[ ])
	}	
	
	for (iMat in 1:NBN) {
		col <- 1:NBN
		col <- col[-iMat]
							
		for (iRow in 1:(NBN-1)) {							# Remove one row and one column from MatR
			pY[iRow] = MatR[iMat, col[iRow]]
								
			for (iCol in 1:NBN-1) {
				pX[iRow, iCol] = MatR[col[iCol], col[iRow]]
			}
		}
							
		Form <- paste("pY[]~pX[,]+0")
		MatrCc[iMat, col, iPert] 	<- (lm(as.formula(Form)))$coefficients
	}		# Loop on iMat
	
	qe2		<- sum((MatrCc[,,iPert] - MatrTh)**2)
	ErrorQ[iPert]	<- sqrt(qe2)
}			# Loop on iPert
	
NomFic 	<- paste(vRoot, "Fig2a_MRA3.pdf", sep="")
pdf(file=NomFic, height=6, width=6) 					# Sizes default to 7

#	par(bg="grey")		# to modify background color
#plot (Perturbs[2:13], ErrorQ[2:13], xlim=c(0.3, 1.8), ylim=c(-0.1,4), main="Error Cc vs. Th", type="l", col="blue", xlab="Perturbation %", ylab="Error", xaxt = "n", axes=TRUE, bty="n")
plot (Perturbs[2:13], ErrorQ[2:13], xlim=c(0.3, 1.8), ylim=c(-0.1,4), main="Error Cc vs. Th", type="b", col="blue", xlab="Perturbation %", ylab="Error", xaxt="n", axes=TRUE, bty="l")
axis(1, at= seq(0.5, 1.5, by=0.5), labels= c("-50", "0", "+50"))		# vs types ="l"
 
segments(x0 = 0, y0 = ErrorQ[3],  x1 = Perturbs[3], y1 = ErrorQ[3], col = "darkgrey", lty="dashed") 
segments(x0 = Perturbs[3], y0 = -0.1, x1 = Perturbs[3], y1 = ErrorQ[3], col = "darkgrey", lty="dashed") 
text(0.3, ErrorQ[3], round(ErrorQ[3], digits=2), font=3)		# font : 1=normal, 2=bold, 3=italic, 4=bold+italic
points(Perturbs[3], ErrorQ[3], pch=16, col="black")
text(Perturbs[3]+0.05, ErrorQ[3]+0.05, "A")

segments(x0 = 0, y0 = 0, x1 = (Perturbs[3]+Perturbs[11])/2, y1 = 0, col = "darkgrey", lty="dashed") 
segments(x0 = (Perturbs[3]+Perturbs[11])/2, y0 = -0.1, x1 = (Perturbs[3]+Perturbs[11])/2, y1 = 0, col = "darkgrey", lty="dashed") 
points((Perturbs[3]+Perturbs[11])/2, 0, pch=16, col="black")
text((Perturbs[3]+Perturbs[11])/2, 0.25, "B")

segments(x0 = 0, y0 = ErrorQ[11],  x1 = Perturbs[11], y1 = ErrorQ[11], col = "darkgrey", lty="dashed") 
segments(x0 = Perturbs[11], y0 = -0.1, x1 = Perturbs[11], y1 = ErrorQ[11], col = "darkgrey", lty="dashed") 
text(0.3, ErrorQ[11], round(ErrorQ[11], digits=2), font=3)		# font : 1=normal, 2=bold, 3=italic, 4=bold+italic
points(Perturbs[11], ErrorQ[11], pch=16, col="black")
text(Perturbs[11], ErrorQ[11]+0.25, "C")

dev.off()


#		Fig. 2b : noise is added (impact of perturbation AND noise level)

NBN 		<- 6								# Number of nodes
Perturbs 	<- c(0.5, 0.7, 0.9, 0.99, 1.01, 1.1, 1.3, 1.5)					# 8 perturbations
PercentPert	<- c("-50", "-30", "-10", "-1", "+1", "+10", "+30", "+50")		# Percentages of perturbation
NBP			<- length(Perturbs)
Noises		<-	c(0.001, 0.005, 0.01)			# Noise  levels
NBT			<- length(Noises)

NBD			<-	100								# Number of random draws

XMesNP		<- array(dim=c(NBN))				# Abundance of the tested product (protein, gene ...) WITHOUT perturbation, WITH noise = abundance measured
XMesP		<- array(dim=c(NBN, NBN))			# Abundance of the tested product (protein, gene ...) WITH perturbation, WITH noise  = abundance measured
MatRN		<- array(dim=c(NBN, NBN))			# Global Response matrix ("R"), with noise
pX   		<- array(dim=c(NBN-1, NBN-1))		# row (iRow), column (iCol)	
pY   		<- array(dim=c(NBN-1))				# row
	
MatrCc		<- array(-1, dim=c(NBN, NBN))		# Computed Connection matrix ("r")
   
ErrorQ		<- array(dim=c(NBD))				# Quadratic error
ErrorQLn	<- array(dim=c(NBD))				# Ln(1+Quad. error)

ErrMoy		<- array(dim=c(NBT, NBP))			# Mean value of ErrorQ (amongst the random draws)
ErrSd		<- array(dim=c(NBT, NBP))			# Standard deviation of ErrorQ (amongst the random draws)
ErrMoyLn	<- array(dim=c(NBT, NBP))			# Mean value of ErrorQLn (amongst the random draws)
ErrSdLn		<- array(dim=c(NBT, NBP))
		 
for (iPert in c(1:NBP)) {
	cat("iPert ", iPert," st ", st, "\n")
	Perturb	 <- Perturbs[iPert]						# Perturbation level
	set.seed(12345)									# So as to generate the same seqences of noise	
	
	for (iTest in c(1:NBT)) {						# Noise level
		sd <- Noises[iTest]*XMesNPMoy				# Noise standard deviation

		for (iDraw in 1:NBD) {
			ss <- multiroot (f=F, start= st)		# Non perturbed measure, with noise
			XMesNP[ ]	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			
			#	Perturbed measures, with noise
			V4 <- Perturb*V4		# P1
			ss <- multiroot (f=F, start= st)
			XMesP[ ,1] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V4 <- V4/Perturb		# Return to initial value	
			
			V3 <- Perturb*V3		# P2
			ss <- multiroot (f=F, start= st)
			XMesP[ ,2] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V3 <- V3/Perturb		# Return to initial value
			
			V8 <- Perturb*V8		# P3
			ss <- multiroot (f=F, start= st)
			XMesP[ ,3]	 <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V8 <- V8/Perturb		# Return to initial value				
			
			V7 <- Perturb*V7		# P4
			ss <- multiroot (f=F, start= st)
			XMesP[ ,4] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V7 <- V7/Perturb		# Return to initial value				
			
			V12 <- Perturb*V12		# P5
			ss <- multiroot (f=F, start= st)
			XMesP[ ,5] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V12 <- V12/Perturb		# Return to initial value				
			
			V11 <- Perturb*V11		# P6
			ss <- multiroot (f=F, start= st)
			XMesP[ ,6] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V11 <- V11/Perturb		# Return to initial valu				
				
			for (iCol in 1:NBN) {
				MatRN[ ,iCol]  <- 2 * (XMesP[ ,iCol]-XMesNP[ ]) / (XMesP[ ,iCol]+XMesNP[ ])
			}

			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
									
				for (iRow in 1:(NBN-1)) {								# Remove one row and one column from MatRN
					pY[iRow] = MatRN[iMat, col[iRow]]
										
					for (iCol in 1:NBN-1) {
						pX[iRow, iCol] = MatRN[col[iCol], col[iRow]]
					}
				}
									
				Form <- paste("pY[]~pX[,]+0")
				MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
			}		# Loop on iMat

			qe2		<- sum((MatrCc - MatrTh)**2)
			ErrorQ  [iDraw]		<- sqrt(qe2)
			ErrorQLn[iDraw]		<- log1p(ErrorQ[iDraw])		# Ln(1+ErrorQ). The function log1p(x= = Ln(1+x) is more accurate than Log(1+x) if |x| << 1
		}			# Loop on iDraw
	
		ErrMoy  [iTest,iPert]		<- mean	(ErrorQ)
		ErrSd	[iTest,iPert]		<- sd  	(ErrorQ)	
		ErrMoyLn[iTest,iPert]		<- mean	(ErrorQLn)
		ErrSdLn	[iTest,iPert]		<- sd  	(ErrorQLn)
	}			# Loop on iTest
}				# Loop on iPert		

matrN2.df	<- data.frame(Pert=1:8, mean=ErrMoyLn[1, 1:8],  lower=ErrMoyLn[1, 1:8]-ErrSdLn[1, 1:8],   upper=ErrMoyLn[1, 1:8]+ErrSdLn[1, 1:8], Coef = "k = 0.002")
matrN5.df	<- data.frame(Pert=1:8, mean=ErrMoyLn[2, 1:8],  lower=ErrMoyLn[2, 1:8]-ErrSdLn[2, 1:8],   upper=ErrMoyLn[2, 1:8]+ErrSdLn[2, 1:8], Coef = "k = 0.005")
matrN10.df	<- data.frame(Pert=1:8, mean=ErrMoyLn[3, 1:8],  lower=ErrMoyLn[3, 1:8]-ErrSdLn[3, 1:8],   upper=ErrMoyLn[3, 1:8]+ErrSdLn[3, 1:8], Coef = "k = 0.01")

NomFic 	<- paste(vRoot, "Fig2b_MRA3.pdf", sep="")
pdf(file=NomFic, height=6, width=6) 					# Sizes default to 7

Titre		<- "Ln (1+Err) as a function of the perturbation"			# Title of the figure
szLBar		<- 7				# Size of the bars  3
spLBar		<- 0.25				# Space between the bars  0.15
szEBar		<- 0.1				# Size of the error bars
ggplot() + 
	geom_linerange(data=matrN2.df,   aes(x=Pert, 			ymin=0, ymax=mean), size=szLBar, color="blue")+ 
	geom_errorbar (data=matrN2.df,   aes(x=Pert, 			ymin=lower, ymax=upper), color = "black", width = szEBar)+
	geom_linerange(data=matrN5.df,   aes(x=Pert+spLBar,		ymin=0, ymax=mean), size=szLBar, color="yellow")+ 
	geom_errorbar (data=matrN5.df,   aes(x=Pert+spLBar,		ymin=lower, ymax=upper), color = "black", width = szEBar)+
	geom_linerange(data=matrN10.df,  aes(x=Pert+2*spLBar, 	ymin=0, ymax=mean), size=szLBar, color="red")+ 
	geom_errorbar (data=matrN10.df,  aes(x=Pert+2*spLBar, 	ymin=lower, ymax=upper), color = "black", width = szEBar)+
	scale_colour_manual(name="Coef", values=c("k = 0.002"="blue", "k = 0.005"="yellow", "k = 0.01"="red"))+
	labs(title=Titre, x="\nPerturbation (%)", y="Ln (1+Err)") + scale_x_continuous(breaks=1:8+spLBar, labels=PercentPert)+
	theme(text = element_text(size = 14)) 
		   
dev.off()
#	To get "Fig2b_MRA3.png", we have to combine this file with "Legende.svg" (Inkscape)


#
#		Figure 3
#		Inference performance for the 6 nodes MAP kinase network
#		Comparison of computation methods, 3 noise levels, perturbation level +50%
#

NBN 		<- 6								# Number of nodes
Perturb 	<- 1.5								# Perturbation level +50%
Noises		<-	c(0.001, 0.005, 0.01)			# Noise  levels
NBT			<- length(Noises)

NBD			<-	100								# Number of random draws
Replics		<-  c(3, 5)							# Number of replicates
NBR			<- length(Replics)
NBRM		<- max(Replics)

Methods		<-  c("MRA", "LSE_CI", "LASSO", "LASSO1", "TLR", "STEP-Fo", "STEP-Ba", "STEP-Bo")	# Methods used (see details inside the functions)
NBM			<- length(Methods)					# Number of méthods
NBRST		<-	3								# Number of STEP methods
ScoreTests	<-	c("T+", "T-", "T0", "F+", "F-", "F0", "P", "N", "Sensib", "Specif")				# Scores to measure
NBSC		<- length(ScoreTests)				# Number of scores

MatrThDig	<- array(dim=c(NBN, NBN))			# Theoretical matix "r", digitalized (-1, 0, 1)
Score		<- array(dim=c(NBSC, NBM, NBT))		# Scores measured
dimnames(Score)	<- list(ScoreTests, Methods, Noises)
DistM 			<- array(0, dim=c(NBM, NBT))	# Distance to the first bisecting line of the point (Se, 1-Sp)
dimnames(DistM) <- list(Methods, Noises)

XMesNP		<- array(dim=c(NBN))				# Abundance of the tested product (protein, gene ...) WITHOUT perturbation, WITH noise = abundance measured
XMesP		<- array(dim=c(NBN, NBN))			# Abundance of the tested product (protein, gene ...) WITH perturbation, WITH noise  = abundance measured
MatRN		<- array(dim=c(NBN, NBN, NBD, NBT))	# Global Response matrix ("R"), with noise
pX   		<- array(dim=c(NBN-1, NBN-1))		# row (iRow), column (iCol)	
pY   		<- array(dim=c(NBN-1))				# row
gX   		<- array(dim=c(NBRM*(NBN-1), NBN-1, NBN, NBT))	# row (iRow), column (iCol), Mat nbr (iMat), Noise (iTest)	 -- 5 replicates
gY   		<- array(dim=c(NBRM*(NBN-1), NBN, NBT))			# row, , Mat nbr (iMat), Noise (iTest)						 -- 5 replicates
gX3   		<- array(dim=c(3*(NBN-1), NBN-1, NBN, NBT))		# row (iRow), column (iCol), , Mat nbr (iMat), Noise (iTest) -- 3 replicates
gY3   		<- array(dim=c(3*(NBN-1), NBN, NBT))			# row, Mat nbr (iMat), Noise (iTest)						 -- 3 replicates

MatrCc		<- array(-1, dim=c(NBN, NBN,NBD))	# Computed Connection matrix ("r")
MatrCcDig	<- array(-1, dim=c(NBN, NBN))		# Computed Connection matrix digitalized ((-1, 0, 1)

ValMeanBP	<- array(dim=c(NBN, NBN, NBT))		# Mean value of rij, bootstrap method
ValMinBP	<- array(dim=c(NBN, NBN, NBT))		# Min value at IC95 of rij, bootstrap method
ValMaxBP	<- array(dim=c(NBN, NBN, NBT))		# Max value at IC95 of rij, bootstrap method

ValMean3	<- array(dim=c(NBN, NBN, NBT))		# Mean value of rij, 3 replicates
ValMin3		<- array(dim=c(NBN, NBN, NBT))		# Min value at IC95 of rij, 3 replicates
ValMax3		<- array(dim=c(NBN, NBN, NBT))		# Max value at IC95 of rij, 3 replicates
ValCc3		<- array(dim=c(NBN, NBN, NBT))		# Coefficient values, 3 replicates

ValMean5	<- array(dim=c(NBN, NBN, NBT))		# Mean value of rij, 5 replicates
ValMin5		<- array(dim=c(NBN, NBN, NBT))		# Min value at IC95 of rij, 5 replicates
ValMax5		<- array(dim=c(NBN, NBN, NBT))		# Max value at IC95 of rij, 5 replicates
ValCc5		<- array(dim=c(NBN, NBN, NBT))		# Coefficient values, 5 replicates

rL 			<- array(dim=c(NBN))				# Result of Lasso method (ie sol. of Yi = Ai * Xi)
rLasso  	<- array(dim=c(NBN, NBN))			# Integration of results Xi, add -1 on the diagonal and transpose the matrix  to get the matrix "r"
rij 		<- matrix(nrow=1, ncol=NBN-1)		# Intermediate computing of rij

set.seed(12345)									# So as to generate the same seqences of noise
Seed	<- .Random.seed

for (iTest in c(1:NBT)) {						# Noise level
	.Random.seed	<- Seed
	sd <- Noises[iTest]*XMesNPMoy				# Noise standard deviation

	for (iDraw in 1:NBD) {
		ss <- multiroot (f=F, start= st)		# Non perturbed measure, with noise
		XMesNP[ ]	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		
		#	Perturbed measures, with noise
		V4 <- Perturb*V4		# P1
		ss <- multiroot (f=F, start= st)
		XMesP[ ,1] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V4 <- V4/Perturb		# Return to initial value	
		
		V3 <- Perturb*V3		# P2
		ss <- multiroot (f=F, start= st)
		XMesP[ ,2] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V3 <- V3/Perturb		# Return to initial value
		
		V8 <- Perturb*V8		# P3
		ss <- multiroot (f=F, start= st)
		XMesP[ ,3]	 <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V8 <- V8/Perturb		# Return to initial value				
		
		V7 <- Perturb*V7		# P4
		ss <- multiroot (f=F, start= st)
		XMesP[ ,4] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V7 <- V7/Perturb		# Return to initial value				
		
		V12 <- Perturb*V12		# P5
		ss <- multiroot (f=F, start= st)
		XMesP[ ,5] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V12 <- V12/Perturb		# Return to initial value				
		
		V11 <- Perturb*V11		# P6
		ss <- multiroot (f=F, start= st)
		XMesP[ ,6] 	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V11 <- V11/Perturb		# Return to initial valu				
			
		for (iCol in 1:NBN) {
			MatRN[ ,iCol,iDraw,iTest]  <- 2 * (XMesP[ ,iCol]-XMesNP[ ]) / (XMesP[ ,iCol]+XMesNP[ ])
		}
	}		# Loop on iDraw
	
	Seed	<- .Random.seed			# Seed is modified by Lasso etc.. in a randomly fashion

#
#		Bootstrap
#
	for (iDraw in 1:NBD) {
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
								
			for (iRow in 1:(NBN-1)) {								# Remove one row and one column from MatRN
				pY[iRow] = MatRN[iMat,col[iRow],iDraw,iTest]
									
				for (iCol in 1:NBN-1) {
					pX[iRow, iCol] = MatRN[col[iCol], col[iRow],iDraw,iTest]
				}
			}
								
			Form <- paste("pY[]~pX[,]+0")
			MatrCc[iMat, col, iDraw] 	<- (lm(as.formula(Form)))$coefficients
		}		# Loop on iMat
	}			# Loop on iDraw
	
	for (iRow in 1:NBN) {
		col <- 1:NBN
		col <- col[-iRow]
		
		for (iCol in col) {	
			ValMeanBP[iRow,iCol,iTest] 		<- mean(MatrCc[iRow,iCol, ])
			qq	<- quantile(MatrCc[iRow,iCol, ], probs=c(0.025,0.975))		# IC 95%
			ValMinBP [iRow,iCol,iTest] 		<- qq[1]
			ValMaxBP [iRow,iCol,iTest] 		<- qq[2]
		}	# Loop on iCol
	}		# Loop on iRow
	
#
#	Use of replicates
#

	for (iDraw in 1:NBRM) {
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
								
			for (iRow in 1:(NBN-1)) {								# Remove one row and one column from MatRN
				pY[iRow] = MatRN[iMat,col[iRow],iDraw,iTest]
									
				for (iCol in 1:NBN-1) {
					pX[iRow, iCol] = MatRN[col[iCol], col[iRow],iDraw,iTest]
				}
			}
		
			gX[((NBN-1)*(iDraw-1)+1) : ((NBN-1)*iDraw), ,iMat,iTest] 	<- pX[ , ]
			gY[((NBN-1)*(iDraw-1)+1) : ((NBN-1)*iDraw),iMat,iTest]  	<- pY[  ]
		}
	}
	
	gX3[ , , ,iTest]	<- gX[1:(3*(NBN-1)), , ,iTest]
	gY3[ , ,iTest]		<- gY[1:(3*(NBN-1)), ,iTest]
	
	for (iRow in 1:NBN) {
		col 	<- 1:NBN
		col 	<- col[-iRow]

		Form3 	<- paste("gY3[ ,", iRow, ",", iTest, "]~gX3[ , ,", iRow, ",", iTest, "]+0")		# 3 replicates
		pQ 		<- lm(as.formula(Form3))
		ValMean3[iRow,col,iTest]  	<- broom::tidy(pQ)$estimate
		ValMin3 [iRow,col,iTest]  	<- broom::tidy(pQ)$estimate - 1.96*broom::tidy(pQ)$std.error	# 1.96 to get an IC 95%
		ValMax3 [iRow,col,iTest]  	<- broom::tidy(pQ)$estimate + 1.96*broom::tidy(pQ)$std.error
		ValCc3  [iRow,col,iTest]  	<- pQ$coefficients		
	
		Form5 	<- paste("gY[ ,", iRow, ",", iTest, "]~gX[ , ,", iRow, ",", iTest, "]+0")			# 5 replicates
		pQ 		<- lm(as.formula(Form5))
		ValMean5[iRow,col,iTest]  	<- broom::tidy(pQ)$estimate
		ValMin5 [iRow,col,iTest]  	<- broom::tidy(pQ)$estimate - 1.96*broom::tidy(pQ)$std.error	# 1.96 to get an IC 95%
		ValMax5 [iRow,col,iTest]  	<- broom::tidy(pQ)$estimate + 1.96*broom::tidy(pQ)$std.error	
		ValCc5  [iRow,col,iTest]  	<- pQ$coefficients				
	}		# Loop on iRow
}			# Loop on iTest


#		Informations for the X axis

OName	<-	"MKKK (1)             MKKK-PP (2)               MKK (3)                 MKK-PP (4)               MAPK (5)                 MAPK-PP (6)"
												# Names of the origins of the edges
Breaks	<-  c(1:(NBN*NBN))

		
#		Fig 3a	: noise level k = 0.1%

VMin	<- min(ValMinBP[ , ,1], ValMin3[ , ,1], ValMin5[ , ,1], na.rm=TRUE)		# Minimum on the Y axis
VMax	<- max(ValMaxBP[ , ,1], ValMax3[ , ,1], ValMax5[ , ,1], na.rm=TRUE)		# Maximum on the Y axis

MatrTh1		<- ifelse (MatrTh == -1, VMin, MatrTh)									# exact values

matrBP.df	<- data.frame(mean=NULL, lower=NULL, upper=NULL)		# To draw figures 3a and 3b
matr3.df	<- data.frame(mean=NULL, lower=NULL, upper=NULL)
matr5.df	<- data.frame(mean=NULL, lower=NULL, upper=NULL)
matrTh1.df	<- data.frame(val=NULL)

for (iMat in 1:NBN) {
	matrBP.tmp.df	<- data.frame(name=1:NBN, mean=as.vector(ValMeanBP[1:NBN,iMat,1]), lower=as.vector(ValMinBP[1:NBN,iMat,1]), upper=as.vector(ValMaxBP[1:NBN,iMat,1]))
	matrBP.tmp.df	<- matrBP.tmp.df[-iMat, ]
	matrBP.df		<- rbind(matrBP.df, matrBP.tmp.df)
	matrBP.tmp.df	<- data.frame(name=" ", mean=VMin, lower=VMin, upper=VMin)
	matrBP.df		<- rbind(matrBP.df, matrBP.tmp.df)

	matr3.tmp.df	<- data.frame(name=1:NBN, mean=as.vector(ValMean3[1:NBN,iMat,1]), lower=as.vector(ValMin3[1:NBN,iMat,1]), upper=as.vector(ValMax3[1:NBN,iMat,1]))
	matr3.tmp.df	<- matr3.tmp.df[-iMat, ]
	matr3.df		<- rbind(matr3.df, matr3.tmp.df)
	matr3.tmp.df	<- data.frame(name=" ", mean=VMin, lower=VMin, upper=VMin)
	matr3.df		<- rbind(matr3.df, matr3.tmp.df)

	matr5.tmp.df	<- data.frame(name=1:NBN, mean=as.vector(ValMean5[1:NBN,iMat,1]), lower=as.vector(ValMin5[1:NBN,iMat,1]), upper=as.vector(ValMax5[1:NBN,iMat,1]))
	matr5.tmp.df	<- matr5.tmp.df[-iMat, ]
	matr5.df		<- rbind(matr5.df, matr5.tmp.df)
	matr5.tmp.df	<- data.frame(name=" ", mean=VMin, lower=VMin, upper=VMin)
	matr5.df		<- rbind(matr5.df, matr5.tmp.df)
	
	matrTh1.tmp.df	<- data.frame(name=1:NBN, val=as.vector(MatrTh1[1:NBN,iMat]))
	matrTh1.tmp.df	<- matrTh1.tmp.df[-iMat, ]
	matrTh1.df		<- rbind(matrTh1.df, matrTh1.tmp.df)
	matrTh1.tmp.df	<- data.frame(name=" ", val=VMin)
	matrTh1.df		<- rbind(matrTh1.df, matrTh1.tmp.df)
}

matrBP.df	<- cbind(ind=1:(NBN*NBN), matrBP.df)
matr3.df	<- cbind(ind=1:(NBN*NBN), matr3.df)
matr5.df	<- cbind(ind=1:(NBN*NBN), matr5.df)
matrTh1.df	<- cbind(ind=1:(NBN*NBN), matrTh1.df)


NomFic 	<- paste(vRoot, "Fig3a_MRA3.pdf", sep="")
pdf(file=NomFic, height=10, width=10) 					# Sizes default to 7

ggplot()+ 
	geom_linerange(data=matr3.df,  aes(x=ind,	   ymin=lower, ymax=upper), size=0.9, color="blue")+
	geom_linerange(data=matr5.df,  aes(x=ind+0.25, ymin=lower, ymax=upper), size=0.9, color="red")+
	geom_linerange(data=matrBP.df, aes(x=ind-0.25, ymin=lower, ymax=upper), size=0.9, color="black")+
	geom_hline(yintercept=0, size=1, color="orange")+
	geom_point(data=matrTh1.df, aes(x=ind, y=val), size=1.8, color="green") +		#vs. 0.9
	ylim(VMin+0.003, VMax+0.3)+														# remove the points corresponding to NA (NA replaced by VMin for computations)								
	labs(title="", x="", y="") + scale_x_continuous(breaks=Breaks)+
	theme(axis.text=element_text(size=12))

dev.off()
#	To get "Fig3a_MRA3.png", we have to combine this file with "Legende.svg" (Inkscape)


#		Fig 3b	: noise level k = 0.5%

VMin	<- min(ValMinBP[ , ,2], ValMin3[ , ,2], ValMin5[ , ,2], na.rm=TRUE)		# Minimum on the Y axis
VMax	<- max(ValMaxBP[ , ,2], ValMax3[ , ,2], ValMax5[ , ,2], na.rm=TRUE)		# Maximum on the Y axis

MatrTh1		<- ifelse (MatrTh == -1, VMin, MatrTh)									# exact values

matrBP.df	<- data.frame(mean=NULL, lower=NULL, upper=NULL)		# To draw figures 3a and 3b
matr3.df	<- data.frame(mean=NULL, lower=NULL, upper=NULL)
matr5.df	<- data.frame(mean=NULL, lower=NULL, upper=NULL)
matrTh1.df	<- data.frame(val=NULL)

for (iMat in 1:NBN) {
	matrBP.tmp.df	<- data.frame(name=1:NBN, mean=as.vector(ValMeanBP[1:NBN,iMat,2]), lower=as.vector(ValMinBP[1:NBN,iMat,2]), upper=as.vector(ValMaxBP[1:NBN,iMat,2]))
	matrBP.tmp.df	<- matrBP.tmp.df[-iMat, ]
	matrBP.df		<- rbind(matrBP.df, matrBP.tmp.df)
	matrBP.tmp.df	<- data.frame(name=" ", mean=VMin, lower=VMin, upper=VMin)
	matrBP.df		<- rbind(matrBP.df, matrBP.tmp.df)

	matr3.tmp.df	<- data.frame(name=1:NBN, mean=as.vector(ValMean3[1:NBN,iMat,2]), lower=as.vector(ValMin3[1:NBN,iMat,2]), upper=as.vector(ValMax3[1:NBN,iMat,2]))
	matr3.tmp.df	<- matr3.tmp.df[-iMat, ]
	matr3.df		<- rbind(matr3.df, matr3.tmp.df)
	matr3.tmp.df	<- data.frame(name=" ", mean=VMin, lower=VMin, upper=VMin)
	matr3.df		<- rbind(matr3.df, matr3.tmp.df)

	matr5.tmp.df	<- data.frame(name=1:NBN, mean=as.vector(ValMean5[1:NBN,iMat,2]), lower=as.vector(ValMin5[1:NBN,iMat,2]), upper=as.vector(ValMax5[1:NBN,iMat,2]))
	matr5.tmp.df	<- matr5.tmp.df[-iMat, ]
	matr5.df		<- rbind(matr5.df, matr5.tmp.df)
	matr5.tmp.df	<- data.frame(name=" ", mean=VMin, lower=VMin, upper=VMin)
	matr5.df		<- rbind(matr5.df, matr5.tmp.df)
	
	matrTh1.tmp.df	<- data.frame(name=1:NBN, val=as.vector(MatrTh1[1:NBN,iMat]))
	matrTh1.tmp.df	<- matrTh1.tmp.df[-iMat, ]
	matrTh1.df		<- rbind(matrTh1.df, matrTh1.tmp.df)
	matrTh1.tmp.df	<- data.frame(name=" ", val=VMin)
	matrTh1.df		<- rbind(matrTh1.df, matrTh1.tmp.df)
}

matrBP.df	<- cbind(ind=1:(NBN*NBN), matrBP.df)
matr3.df	<- cbind(ind=1:(NBN*NBN), matr3.df)
matr5.df	<- cbind(ind=1:(NBN*NBN), matr5.df)
matrTh1.df	<- cbind(ind=1:(NBN*NBN), matrTh1.df)


NomFic 	<- paste(vRoot, "Fig3b_MRA3.pdf", sep="")
pdf(file=NomFic, height=10, width=10) 					# Sizes default to 7

Titre		<- paste("rij distribution : k = ", 100*Noises[2], "%")					# ri,j distribution
ggplot()+ 
	geom_linerange(data=matr3.df,  aes(x=ind,	   ymin=lower, ymax=upper), size=0.9, color="blue")+
	geom_linerange(data=matr5.df,  aes(x=ind+0.25, ymin=lower, ymax=upper), size=0.9, color="red")+
	geom_linerange(data=matrBP.df, aes(x=ind-0.25, ymin=lower, ymax=upper), size=0.9, color="black")+
	geom_hline(yintercept=0, size=1, color="orange")+
	geom_point(data=matrTh1.df, aes(x=ind, y=val), size=1.8, color="green") +		# vs 0.9
	ylim(VMin+0.003, VMax+0.3)+														# remove the points corresponding to NA (NA replaced by VMin for computations)								
#	labs(title=Titre, x="", y="CI 95%") + scale_x_continuous(breaks=Breaks, labels=DName)+	# y = "IC 95%" for ri,j
	labs(title="", x="", y="") + scale_x_continuous(breaks=Breaks, labels=matrBP.df$name)+
	xlab(OName)+
	theme(axis.text=element_text(size=12))
	
dev.off()


	
###		
###		Remark : there was an error on the X axis of Fig 3 a and b, concerning the direction of the edges (initial version of the document). 
###		This error is corrected above and the old version is as follows (comments).
###		
###		ConvertON	<- c(30,7,13,19,25,31,1,6,14,20,26,32,2,8,12,21,27,33,3,9,15,18,28,34,4,10,16,22,24,35,5,11,17,23,29,36)
###		Convert		<- c(7,13,19,25,31,8,2,14,20,26,32,15,3,9,21,27,33,22,4,10,16,28,34,29,5,11,17,23,35,1,6,12,18,24,30,36)	
###		
###		vMeanBP	<- as.vector(ValMeanBP[ , ,1])
###		vMinBP	<- as.vector(ValMinBP [ , ,1])
###		vMaxBP	<- as.vector(ValMaxBP [ , ,1])
###		vMean3	<- as.vector(ValMean3 [ , ,1])
###		vMin3	<- as.vector(ValMin3  [ , ,1])
###		vMax3	<- as.vector(ValMax3  [ , ,1])
###		vMean5	<- as.vector(ValMean5 [ , ,1])
###		vMin5	<- as.vector(ValMin5  [ , ,1])
###		vMax5	<- as.vector(ValMax5  [ , ,1])
###		
###		vMeanBPA <- as.vector(ValMeanBP[ , ,1])
###		vMinBPA	 <- as.vector(ValMinBP [ , ,1])
###		vMaxBPA	 <- as.vector(ValMaxBP [ , ,1])
###		vMean3A	 <- as.vector(ValMean3 [ , ,1])
###		vMin3A	 <- as.vector(ValMin3  [ , ,1])
###		vMax3A	 <- as.vector(ValMax3  [ , ,1])
###		vMean5A	 <- as.vector(ValMean5 [ , ,1])
###		vMin5A	 <- as.vector(ValMin5  [ , ,1])
###		vMax5A	 <- as.vector(ValMax5  [ , ,1])
###		
###		for (i in 1: (NBN*NBN))	{
###			vMeanBP[i]	<- vMeanBPA[Convert[i]]
###			vMinBP [i]	<- vMinBPA [Convert[i]]
###			vMaxBP [i]	<- vMaxBPA [Convert[i]]
###			vMean3 [i]	<- vMean3A [Convert[i]]
###			vMin3  [i]	<- vMin3A  [Convert[i]]
###			vMax3  [i]	<- vMax3A  [Convert[i]]
###			vMean5 [i]	<- vMean5A [Convert[i]]
###			vMin5  [i]	<- vMin5A  [Convert[i]]
###			vMax5  [i]	<- vMax5A  [Convert[i]]
###		}	
###			
###		NmatrBP.df	<- data.frame(ind=1:(NBN*NBN), name=NameX, mean=vMeanBP, lower=vMinBP, upper=vMaxBP)
###		Nmatr3.df	<- data.frame(ind=1:(NBN*NBN), name=NameX, mean=vMean3,  lower=vMin3,  upper=vMax3)
###		Nmatr5.df	<- data.frame(ind=1:(NBN*NBN), name=NameX, mean=vMean5,  lower=vMin5,  upper=vMax5) 
###		NmatrTh1.df	<- data.frame(ind=1:(NBN*NBN), name=NameX, val=as.vector(MatrTh1))
###		
###		NomFic 	<- paste(vRoot, "Fig3aOld_MRA3.pdf", sep="")
###		pdf(file=NomFic, height=10, width=10) 					# Sizes default to 7
###		
###		Titre		<- paste("rij distribution : k = ", Noise[iTest])											# ri,j distribution
###		ggplot()+ 
###			geom_linerange(data=Nmatr3.df,  aes(x=ind,	   ymin=lower, ymax=upper), size=0.9, color="blue")+
###			geom_linerange(data=Nmatr5.df,  aes(x=ind+0.25, ymin=lower, ymax=upper), size=0.9, color="red")+
###			geom_linerange(data=NmatrBP.df, aes(x=ind-0.25, ymin=lower, ymax=upper), size=0.9, color="black")+
###			geom_hline(yintercept=0, size=1, color="orange")+
###			geom_point(data=NmatrTh1.df, aes(x=ind, y=val), size=0.7, color="green") +
###			labs(title=Titre, x="", y="CI 95%") + scale_x_continuous(breaks=Breaks, labels=NameX)	# y = "IC 95%" for ri,j
###			
###		dev.off()
	

#		Fig 3b	: noise level k = 0.5%

VMin	<- min(ValMinBP[ , ,2], ValMin3[ , ,2], ValMin5[ , ,2], na.rm=TRUE)		# Minimum on the Y axis
VMax	<- max(ValMaxBP[ , ,2], ValMax3[ , ,2], ValMax5[ , ,2], na.rm=TRUE)		# Maximum on the Y axis

ValMeanBP[ , ,2]	<- ifelse (is.na(ValMeanBP[ , ,2]),	 VMin, ValMeanBP[ , ,2])	# Replace the values on the diagonal (which were NA) with VMin
ValMinBP [ , ,2]	<- ifelse (is.na(ValMinBP[ , ,2]), 	 VMin, ValMinBP	[ , ,2])
ValMaxBP [ , ,2]	<- ifelse (is.na(ValMaxBP[ , ,2]),	 VMin, ValMaxBP	[ , ,2])

ValMean3[ , ,2]		<- ifelse (is.na(ValMean3[ , ,2]),	 VMin, ValMean3[ , ,2])
ValMin3 [ , ,2]		<- ifelse (is.na(ValMin3[ , ,2]), 	 VMin, ValMin3 [ , ,2])
ValMax3 [ , ,2]		<- ifelse (is.na(ValMax3[ , ,2]),	 VMin, ValMax3 [ , ,2])

ValMean5[ , ,2]		<- ifelse (is.na(ValMean5[ , ,2]),	 VMin, ValMean5[ , ,2])
ValMin5 [ , ,2]		<- ifelse (is.na(ValMin5[ , ,2]), 	 VMin, ValMin5 [ , ,2])
ValMax5 [ , ,2]		<- ifelse (is.na(ValMax5[ , ,2]),	 VMin, ValMax5 [ , ,2])

MatrTh1		<- ifelse (MatrTh == -1, VMin, MatrTh)									# exact values

matrBP.df	<- data.frame(ind=1:(NBN*NBN), name=DName, mean=as.vector(ValMeanBP[ , ,2]), lower=as.vector(ValMinBP[ , ,2]), upper=as.vector(ValMaxBP[ , ,2]))
matr3.df	<- data.frame(ind=1:(NBN*NBN), name=DName, mean=as.vector(ValMean3[ , ,2]),  lower=as.vector(ValMin3[ , ,2]),  upper=as.vector(ValMax3[ , ,2]))
matr5.df	<- data.frame(ind=1:(NBN*NBN), name=DName, mean=as.vector(ValMean5[ , ,2]),  lower=as.vector(ValMin5[ , ,2]),  upper=as.vector(ValMax5[ , ,2]))
matrTh1.df	<- data.frame(ind=1:(NBN*NBN), name=DName, val=as.vector(MatrTh1))

NomFic 	<- paste(vRoot, "Fig3b_MRA3.pdf", sep="")
pdf(file=NomFic, height=10, width=10) 					# Sizes default to 7

Titre		<- paste("rij distribution : k = ", 100*Noises[2], "%")					# ri,j distribution
ggplot()+ 
	geom_linerange(data=matr3.df,  aes(x=ind,	   ymin=lower, ymax=upper), size=0.9, color="blue")+
	geom_linerange(data=matr5.df,  aes(x=ind+0.25, ymin=lower, ymax=upper), size=0.9, color="red")+
	geom_linerange(data=matrBP.df, aes(x=ind-0.25, ymin=lower, ymax=upper), size=0.9, color="black")+
	geom_hline(yintercept=0, size=1, color="orange")+
	geom_point(data=matrTh1.df, aes(x=ind, y=val), size=0.7, color="green") +
	ylim(VMin+0.3, VMax+0.3)+														# remove the points corresponding to NA (NA replaced by VMin for computations)		
	labs(title=Titre, x="", y="CI 95%") + scale_x_continuous(breaks=Breaks, labels=DName)+	# y = "IC 95%" for ri,j
	xlab(OName)
	
dev.off()
#	To get "Fig3b_MRA3.png", we have to combine this file with "Legende.svg" (Inkscape)


#		Fig. 3c : comparison of performance estimation (sensitivity and specificity) according to the computation methods

MatrThDig	<-	ifelse (MatrTh>0, 1, MatrTh)
MatrThDig	<-	ifelse (MatrThDig<0, -1, MatrThDig)

for (iTest in c(1:NBT)) {						# Noise level
	for (iMeth in 1:NBM) {						# Computation methods
		Method	<- Methods[iMeth]
		cat ("Noise ", Noises[iTest], " Method ", Method, "\n")

		if (Method == "MRA") {					# We compare IC95 with 0
			MatrCcDig[]		<- ifelse (ValMinBP[ , ,iTest] > 0,  1, 0)
			MatrCcDig[]		<- ifelse (ValMaxBP[ , ,iTest] < 0, -1, MatrCcDig)	
			Score[ ,iMeth,iTest]	<- Benchmark2(MatrCcDig, MatrThDig)
		}	# Method MRA

		if (Method == "LSE_CI") {				# We compare IC95 with 0  -- 3 replicates
			MatrCcDig[]		<- ifelse (ValMin3[ , ,iTest] > 0,  1, 0)
			MatrCcDig[]		<- ifelse (ValMax3[ , ,iTest] < 0, -1, MatrCcDig)	
			Score[ ,iMeth,iTest]	<- Benchmark2(MatrCcDig, MatrThDig)
		}	# Method LSE-CI

		if (Method == "LASSO") {				# Digitalization = sign of the result	-- 3 replicates  -- automatic choice of lambda
			rLasso[ , ]	<- 0
			for (iMat in 1:NBN) {
				col 	<- 1:NBN
				col 	<- col[-iMat]

				cv_model 	<- cv.glmnet(gX3[ , ,iMat,iTest], gY3[ ,iMat,iTest], alpha = 1)	# Fit lasso regression model using k-fold cross-validation
				best_lambda <- cv_model$lambda.min
				best_model 	<- glmnet(gX3[ , ,iMat,iTest], gY3[ ,iMat,iTest], alpha = 1, lambda = best_lambda)		# View coefficients of best model
				rL			<- coef(best_model)
				rLasso[iMat, col]	<- rL[2:NBN, 1]		# we get the matrix "r"	
			}
			diag(rLasso)	<- -1
			rLasso			<- ifelse (is.na(rLasso), 0, rLasso)
				
			MatrCcDig[]		<- ifelse (rLasso > 0,  1, 0)
			MatrCcDig[]		<- ifelse (rLasso < 0, -1, MatrCcDig)	
			Score[ ,iMeth,iTest]	<- Benchmark2(MatrCcDig, MatrThDig)
		}	# Method LASSO

		if (Method == "LASSO1") {				# Digitalization = sign of the result	-- 3 replicates  -- Optimization of the hyper-parameters
			rLasso[ , ]	<- 0			
			for (iMat in 1:NBN) {
				col 	<- 1:NBN
				col 	<- col[-iMat]

				cv_model 	<- cv.glmnet(gX3[ , ,iMat,iTest], gY3[ ,iMat,iTest], alpha = 1, intercept = FALSE, nfolds=18)	# Fit lasso regression model using k-fold cross-validation
				rL <- coef(cv_model, cv_model$lambda.1se)
				rLasso[iMat, col]	<- rL[2:NBN, 1]		# we get the matrix "r"	
			}
			diag(rLasso)	<- -1
			rLasso		<- ifelse (is.na(rLasso), 0, rLasso)
				
			MatrCcDig[]		<- ifelse (rLasso > 0,  1, 0)
			MatrCcDig[]		<- ifelse (rLasso < 0, -1, MatrCcDig)	
			Score[ ,iMeth,iTest]	<- Benchmark2(MatrCcDig, MatrThDig)
		}	# Method LASSO1
		
		if (Method == "TLR") {
			ValCc3[ , ,iTest]	<- ifelse (is.na(ValCc3[ , ,iTest]), 0, ValCc3[ , ,iTest])
			ThTLR	<- 0.25 * max(ValCc3[ , ,iTest])	
		
			MatrCcDig[]		<- ifelse (ValCc3[ , ,iTest] >=  ThTLR,   1, 0)
			MatrCcDig[]		<- ifelse (ValCc3[ , ,iTest] <= -ThTLR,  -1, MatrCcDig[])
			Score[ ,iMeth,iTest]	<- Benchmark2(MatrCcDig, MatrThDig)	
		}	# Method TLR
		
		if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {	
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
		
				Donnees  <- data.frame(gY3[ ,iMat,iTest], gX3[ , ,iMat,iTest])
				#	cat("iMat ", iMat, " Donnees ", colnames(Donnees), "\n")
				Donnees	 <- rename(Donnees, "Y" ="gY3...iMat..iTest.")		# To simplify coding
				cc 		 <- colnames(Donnees)				# Name of the columns of "Donnees". The first one is "Y"
				cc		 <- cc[-1]							# Delete "Y". The name of the coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1) remains only	
				colnames(rij) <- cc									# step indique le nom des colonnes correspondant aux coefficients conservés
						
				switch(Method, 
					   "STEP-Fo" =
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},	# Forward
					   "STEP-Ba" = 
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward
					   "STEP-Bo" =
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)},		# Both
					)
				
				ll  	<- length(forward$coefficients)				# Number of coefficients kept by STEP
				nn		<- names(forward$coefficients)				# Name of these coefficients
		
				rij[1,]	<- 0
				if (ll >= 2) {
					for (i in 2:ll) {								# nn[1] = "(Intercept)"
						rij[1, nn[i]] <- forward$coefficients[nn[i]]
					}
				}
				if (Verbose) {
					cat (" iMat ", iMat, " ll ", ll, " nn ", nn, " rij ", rij, "\n")
				}
				
				MatrCcDig[iMat, col]	<- rij[1,]
				MatrCcDig[iMat, iMat]  	<- 0
			}	# Loop on iMat

			MatrCcDig	<-	ifelse (MatrCcDig > 0,  1, MatrCcDig)
			MatrCcDig	<-	ifelse (MatrCcDig < 0, -1, MatrCcDig)
			Score[ ,iMeth,iTest]	<- Benchmark2(MatrCcDig, MatrThDig)	
		}	# Methods STEP
		
		DistM[iMeth,iTest]		<- Score["Sensib",iMeth,iTest] - (1-Score["Specif",iMeth,iTest])
		cat(Methods[iMeth], " Noise ", Noises[iTest], " DistM ", DistM[iMeth,iTest], "\n")
	}	# Loop on iMeth
	
	write.csv(Score[ , ,iTest], 	paste(vRoot, "Fig3c_", Noises[iTest], ".csv", sep=""))	
	write.csv(t(Score[ , ,iTest]), 	paste(vRoot, "Fig3c_", Noises[iTest], "_Transp.csv", sep=""))	# File used also by "Supplementary Tables & Figures", Supplementary Table 1		
}		# Loop on iTest


#
#		Figure 4
#		Dream Challenge 4 networks
#
#	To get the data files, go to the site "https://www.synapse.org/#!Synapse:syn3049712/wiki/74630" and download :

##	DREAM4_InSilico_Size10.zip
##	DREAM4_InSilico_Size100.zip
##	Download too :
##	DREAM4_InSilicoNetworks_GoldStandard.zip, thru "https://www.synapse.org/#!Synapse:syn3049736"
##	(Open "DREAM4_Challenge2_GoldStandards", then "Size10", " Size 10 bonus round" and at last " "DREAM4_GoldStandard_InSilico_Size10_1.tsv" 
##	etc… for the other files).
#

NBFIC		<- 5								# There are 5 files for a given size
Sizes		<- c(10, 100)						# Nbr of nodes of the network
NBS			<- length(Sizes) 

Methods		<-  c("MRA-KD", "MRA-KO", "LSE_CI", "LASSO", "LASSO1", "TLR", "STEP-Fo", "STEP-Ba", "STEP-Bo", "CLR", "ARACNE", "MRNET")
NBM			<- length(Methods)					# Number of méthods

ScoreTests	<-	c("TP", "TN", "FP", "FN", "Se", "Sp")		# Scores to measure (Se = "Sensibility", Sp = "Specificity")
NBSC		<- length(ScoreTests)				# Number of scores
ScoreMoyC	<-	c("Se", "Sp")					# Global results for the 5 files of a given size ("Sensibility", "Specificity")
ScoreMoyR	<-	c("mean", "sd")					# Values to display
Score 		<- array(0, dim=c(NBSC, NBM, NBFIC, NBS))		# TP, TN, FP, FN
dimnames(Score)		<- list(ScoreTests, Methods, NULL)
ScoreMoy	<- array(0, dim=c(length(ScoreMoyC), length(ScoreMoyR), NBM, NBS))	# Mean of the scores of a set of files
dimnames(ScoreMoy)	<- list(ScoreMoyC, ScoreMoyR, Methods, Sizes)

vTOP		<- c(0.2, 0.25)						# "top" of the edges to take into account
vPValue		<- 0.05								# pValue expected must be < vPValue

NBK			<- 2								# NBK corresponds to the number of perturbation types (KO or KD)
		
		
#
#		Fig. 4a : drawing of the network topology (insilico_Size10_1), using TLR
#

iSize		<- 1				#insilico_Size10

NBN			<- Sizes[iSize]

Solution	<- array(dim=c(NBN, NBN))			# Solution : the true network structure
Solution2	<- array(dim=c(NBN, NBN))			# Idem "Solution", but row and column numbers correspond to numebers of "levels" ("G1", "G10", "G2", "G3" etc)
Solution3	<- array(dim=c(NBN, NBN))			# Intermediary of calculation

XMesNP		<- array(dim=c(NBN, NBK))			# Measured values, non pertubed system
XMesP		<- array(dim=c(NBN, NBN, NBK))		# Measured values following a perturbation	
	
MatR		<- array(dim=c(NBN, NBN, NBK))		# Matrix R
pX   		<- array(dim=c(NBN-1, NBN-1))				# row (iRow), column (iCol)
pY   		<- array(dim=c(NBN-1))						# row
gX   		<- array(dim=c(NBK*(NBN-1), NBN-1, NBN))	# row (iRow), column (iCol), matrix number (iMat)
gY   		<- array(dim=c(NBK*(NBN-1), NBN))			# row,  matrix number
	
MatrCc		<- array(dim=c(NBN, NBN))			# "r" matrix computed (according to method MRA, TLR, STEP or LASSO)
MatrCcDig	<- array(dim=c(NBN, NBN))			# MatrCc digitalized

ResultNames	<- c("TP", "TN", "FP", "FN", "Se", "Sp")
Results		<- vector(length=length(ResultNames))
# names(Results)	<- ResultNames

iFic		<- 1			# insilico_Size10_1

#	1/ Data reading

val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_wildtype.tsv", sep="")
XMesNPLu 	<- fread(val, data.table=F)		# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_wildtype.tsv"
val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockdowns.tsv", sep="")
XMesPLu1 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_knockdowns.tsv"
val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockouts.tsv", sep="")
XMesPLu2 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1//insilico_size10_1_knockouts.tsv"		
val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/DREAM4_GoldStandard_InSilico_Size", NBN, "_", iFic, ".tsv", sep="")
SolutLu  	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
SolutLu$V1	<- as.factor(SolutLu$V1)
SolutLu$V2	<- as.factor(SolutLu$V2)

for (i in 1:NBN) {
	XMesNP[i,1] = XMesNPLu[1,i]
	XMesNP[i,2] = XMesNPLu[1,i]
	
	for (j in 1:NBN) {
		XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation level -50%	
											# ATTENTION Dream Challenge data are transposed vs. my data 
		XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation level -100%
	}
}

Solution2[,] = 0
for (i in 1:(NBN*(NBN-1))) {
	if(SolutLu[i, ]$V3 == 1) {
		Solution2[as.numeric(SolutLu[i, ]$V1), as.numeric(SolutLu[i, ]$V2)] = 1
	} else {
		break								# We suppose that all the "1" are situated at the beginning of the file
	}
}

numero <- as.numeric(sort(as.character(c(1:NBN))))
for (i in 1:NBN) {
	Solution3[ ,numero[i]] = Solution2[ ,i]
}
for (i in 1:NBN) {
	Solution[ ,numero[i]] = Solution3[i,]	# We find the classical connection matrix (nodes "G1", "G2", ... , "G10")
											# Note that Solution is transposed
}
diag(Solution)	<- 0						# the diagonal is set to 0 (no auto-interaction in this network)

#	2/ Drawing of the solution

nbRows			<- sum(ifelse(Solution>0, 1, 0))	# Number of edges of this network
DessMAP			<- matrix(nrow=nbRows, ncol=4)		# Matrix csv used by Cytoscape to draw the network
nodeNames		<- as.character (seq(1, NBN, by=1))	# Name of the nodes = their number (default)
ColNames		<- c("Source", "Target", "Interaction", "Value")		# The first line corresponds to the names
													# Source,target : name of the corresponding nodes
													# Interaction	: "d" (direct), "r" (reverse), "fp" (false positive), "fn" (false negative)
													# Value : value associated with the edge
colnames(DessMAP)	<- ColNames

iSol			<- 0								# Index of the row added in DessMAP
for (iRow in 1:NBN) {
	Col <- 	1:NBN
	Col <-	Col[-iRow]
	for (iCol in Col) {
		if (Solution[iRow, iCol] != 0) {
			iSol	<- iSol + 1
			DessMAP[iSol, "Source"]			<- nodeNames[iCol]
			DessMAP[iSol, "Target"]			<- nodeNames[iRow]			
			DessMAP[iSol, "Value"]			<- Solution[iRow, iCol]
			DessMAP[iSol, "Interaction"]	<- "d"
		}
	}	# iCol
}		# iRow	

write.csv(DessMAP, paste(vRoot, "silico_10_1_Sol.csv", sep=""), quote=FALSE)	
														# File used by Cytoscape v 3.9.1 to draw the figure 1a (import network from file)
														# Use style "MRA_styles.xml"
														# The picture is exported as "***.pdf" and converted to "***.png" (600 dpi) by Inkscape v 1.2.

#	3/ Discovery of the network

for (iPaq in c(1:NBK)) {					# We take every available data : Knock Down and Knock Out
	#	The "input data", delivered by Dream Challenge, have been loaded before in XMesNP and XMesP
	MatR[ , ,iPaq]  <- 2 * (XMesP[ , ,iPaq]-XMesNP[ ,iPaq]) / (XMesP[ , ,iPaq]+XMesNP[ ,iPaq])
	
	for (iMat in 1:NBN) {
		col <- 1:NBN
		col <- col[-iMat]
		
		for (iRow in 1:(NBN-1)) {
			pY[iRow] = MatR[iMat, col[iRow], iPaq]
			
			for (iCol in 1:NBN-1) {
				pX[iRow, iCol] = MatR[col[iCol], col[iRow], iPaq]
			}
		}
	
		gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,iMat]  	<- pX[ , ]
		gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),iMat]   	<- pY[ ]
	}
}

#		Method : TLR

for (iMat in 1:NBN) {
	col <- 1:NBN
	col <- col[-iMat]

	Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")
	MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
	MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
}
MatrCc[is.na(MatrCc)]  <- 0	
MatrCcDig	<-	Digitalize (MatrCc, "Threshold", 0.25)

#	4/ Drawing of the matrix discovered
											
Results		<-	Benchmark(MatrCcDig, Solution)
names(Results)	<- ResultNames
												
nbRows			<- Results["TP"] + Results["FP"] + Results["FN"]	# Number of edges of this network
DessMAP			<- matrix(nrow=nbRows, ncol=4)		# Matrix csv used by Cytoscape to draw the network
nodeNames		<- as.character (seq(1, NBN, by=1))	# Name of the nodes = their number (default)
ColNames		<- c("Source", "Target", "Interaction", "Value")		# The first line corresponds to the names
													# Source,target : name of the corresponding nodes
													# Interaction	: "d" (direct), "r" (reverse), "fp" (false positive), "fn" (false negative)
													# Value : value associated with the edge
colnames(DessMAP)	<- ColNames

iSol			<- 0								# Index of the row added in DessMAP
for (iRow in 1:NBN) {
	Col <- 	1:NBN
	Col <-	Col[-iRow]
	for (iCol in Col) {
		if (Solution[iRow, iCol]>0 & MatrCcDig[iRow, iCol]>0) {
			iSol	<- iSol + 1
			# cat("TP iRow",iRow," iCol", iCol, " iSol ", iSol, "\n")			
			DessMAP[iSol, "Source"]			<- nodeNames[iCol]
			DessMAP[iSol, "Target"]			<- nodeNames[iRow]			
			DessMAP[iSol, "Value"]			<- Solution[iRow, iCol]
			DessMAP[iSol, "Interaction"]	<- "d"		# TP
		}
		
		if (Solution[iRow, iCol]>0 & MatrCcDig[iRow, iCol]==0) {
			# cat("FN iRow",iRow," iCol", iCol, " iSol ", iSol, "\n")		
			iSol	<- iSol + 1
			DessMAP[iSol, "Source"]			<- nodeNames[iCol]
			DessMAP[iSol, "Target"]			<- nodeNames[iRow]			
			DessMAP[iSol, "Value"]			<- Solution[iRow, iCol]
			DessMAP[iSol, "Interaction"]	<- "fn"		# FN
		}
		
		if (Solution[iRow, iCol]==0 & MatrCcDig[iRow, iCol]>0) {
			iSol	<- iSol + 1
			# cat("FP iRow",iRow," iCol", iCol, " iSol ", iSol, "\n")
			DessMAP[iSol, "Source"]			<- nodeNames[iCol]
			DessMAP[iSol, "Target"]			<- nodeNames[iRow]			
			DessMAP[iSol, "Value"]			<- Solution[iRow, iCol]
			DessMAP[iSol, "Interaction"]	<- "fp"		# FP
		}			
	}	# iCol
}		# iRow	

write.csv(DessMAP, paste(vRoot, "silico_10_1_Res.csv", sep=""), quote=FALSE)
														# File used by Cytoscape v 3.9.1 to draw the figure 1a (import network from file)
														# Use style "MRA_styles.xml"
														# The picture is exported as "***.pdf" and converted to "***.png" (600 dpi) by Inkscape v 1.2.												
		
		
#
#		Fig. 4b	: average performance of the methods
#
		
for (iSize in 1:NBS) {
	NBN			<- Sizes[iSize]
	
	Solution	<- array(dim=c(NBN, NBN))			# Solution : the true network structure
	Solution2	<- array(dim=c(NBN, NBN))			# Idem "Solution", but row and column numbers correspond to numebers of "levels" ("G1", "G10", "G2", "G3" etc)
	Solution3	<- array(dim=c(NBN, NBN))			# Intermediary of calculation
	
	XMesNP		<- array(dim=c(NBN, NBK))			# Measured values, non pertubed system
	XMesP		<- array(dim=c(NBN, NBN, NBK))		# Measured values following a perturbation	
		
	MatR		<- array(dim=c(NBN, NBN, NBK))		# Matrix R
	MatRm1		<- array(dim=c(NBN, NBN, NBK))	# Matrix R	** -1
	dgMatRm1M1	<- array(dim=c(NBN, NBN, NBK))	# Diagonal (R** -1) **-1
	pX   		<- array(dim=c(NBN-1, NBN-1))				# row (iRow), column (iCol)
	pY   		<- array(dim=c(NBN-1))						# row
	gX   		<- array(dim=c(NBK*(NBN-1), NBN-1, NBN))	# row (iRow), column (iCol), matrix number (iMat)
	gY   		<- array(dim=c(NBK*(NBN-1), NBN))			# row,  matrix number
	gXO 		<- array(dim=c(NBN-1, NBN-1, NBN))			# gX matrix truncated so as to keep the values KO only
	gYO 		<- array(dim=c(NBN-1, NBN))					# gY matrix truncated so as to keep the values KO only
	gXD 		<- array(dim=c(NBN-1, NBN-1, NBN))			# gX matrix truncated so as to keep the values KD only
	gYD 		<- array(dim=c(NBN-1, NBN))					# gY matrix truncated so as to keep the values KD only
		
	MatrCc		<- array(dim=c(NBN, NBN))			# "r" matrix computed (according to method MRA, TLR, STEP or LASSO)
	MatrPVal	<- array(dim=c(NBN, NBN))			# p-value
	MatrCcDig	<- array(dim=c(NBN, NBN))			# Computed Connection matrix digitalized (0, 1)	
	
	rij 		<- matrix(nrow=1, ncol=NBN-1)		# Intermediate calculation of rij
	rL 			<- array(dim=c(NBN))				# Result of Lasso method (ie sol. of Yi = Ai * Xi)
	
	for (iFic in 1:NBFIC) {
	
		#	1/ Data reading
		
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_wildtype.tsv", sep="")
		XMesNPLu 	<- fread(val, data.table=F)		# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_wildtype.tsv"
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockdowns.tsv", sep="")
		XMesPLu1 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_knockdowns.tsv"
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockouts.tsv", sep="")
		XMesPLu2 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1//insilico_size10_1_knockouts.tsv"		
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/DREAM4_GoldStandard_InSilico_Size", NBN, "_", iFic, ".tsv", sep="")
		SolutLu  	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
		SolutLu$V1	<- as.factor(SolutLu$V1)
		SolutLu$V2	<- as.factor(SolutLu$V2)

		for (i in 1:NBN) {
			XMesNP[i,1] = XMesNPLu[1,i]
			XMesNP[i,2] = XMesNPLu[1,i]
			
			for (j in 1:NBN) {
				XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation level -50%	
													# ATTENTION Dream Challenge data are transposed vs. my data 
				XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation level -100%
			}
		}
	
		Solution2[,] = 0
		for (i in 1:(NBN*(NBN-1))) {
			if(SolutLu[i, ]$V3 == 1) {
				Solution2[as.numeric(SolutLu[i, ]$V1), as.numeric(SolutLu[i, ]$V2)] = 1
			} else {
				break								# We suppose that all the "1" are situated at the beginning of the file
			}
		}
	
		numero <- as.numeric(sort(as.character(c(1:NBN))))
		for (i in 1:NBN) {
			Solution3[ ,numero[i]] = Solution2[ ,i]
		}
		for (i in 1:NBN) {
			Solution[ ,numero[i]] = Solution3[i,]	# We find the classical connection matrix (nodes "G1", "G2", ... , "G10")
													# Note that Solution is transposed
		}

		for (iPaq in c(1:NBK)) {					# We take every available data : Knock Down and Knock Out
			#	The "input data", delivered by Dream Challenge, have been loaded before in XMesNP and XMesP
			MatR[ , ,iPaq]  <- 2 * (XMesP[ , ,iPaq]-XMesNP[ ,iPaq]) / (XMesP[ , ,iPaq]+XMesNP[ ,iPaq])
			
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
				
				for (iRow in 1:(NBN-1)) {
					pY[iRow] = MatR[iMat, col[iRow], iPaq]
					
					for (iCol in 1:NBN-1) {
						pX[iRow, iCol] = MatR[col[iCol], col[iRow], iPaq]
					}
				}
			
				gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,iMat]  	<- pX[ , ]
				gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),iMat]   	<- pY[ ]
			}		# Loop on iMat
		}			# Loop on iPaq
	
		gXD 	<- gX[1 : (NBN-1), , ] 					# KnockDown
		gYD		<- gY[1 : (NBN-1), ]	
		gXO		<- gX[NBN : ((NBN-1)*2), , ] 			# KnockOut
		gYO		<- gY[NBN : ((NBN-1)*2), ] 
		
		for	(iMeth in 1:NBM) {
			Method  <- Methods[iMeth]
		
			cat("iSize ", iSize, " iFic ", iFic, " Method ", Method, " Heure ", as.character(Sys.time()), "\n")
			
			if (Method %in% c("MRA-KD", "MRA-KO")) {	#	Classical MRA method
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
					
					if (Method == "MRA-KD") {
						Form <- paste("gYD[ ,",iMat,"]~gXD[ , ,",iMat,"]+0")
					} else {
						Form <- paste("gYO[ ,",iMat,"]~gXO[ , ,",iMat,"]+0")
					}	
					
					MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}
				MatrCc[is.na(MatrCc)]  <- 0
			}	# MRA methods

			if (Method == "LSE_CI") {
				MatrCc[ , ] <- 0
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
			
					Form <- paste("gY[ ,",iMat,"]~gX[ , ,",iMat,"]+0")
					pQ 	 <- lm(as.formula(Form))								# MatrCc holds the associated p-values
					MatrPVal[iMat,col] 	<- broom::tidy(pQ)$p.value
					MatrCc [iMat,which(MatrPVal[iMat,] < vPValue)] <- 1
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}
				MatrCc[is.na(MatrCc)]  <- 0	
			}	# LSE_CI method	
	
			if (Method == "TLR") {
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
			
					Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")
					MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}
				MatrCc[is.na(MatrCc)]  <- 0	
			}	# TLR with Linear regression, 2 measures per node		

			if (Method == "LASSO") {				# Digitalization = 1 if the result != 0	-- automatic choice of lambda
				MatrCc[ , ] <- 0
				for (iMat in 1:NBN) {
					col 	<- 1:NBN
					col 	<- col[-iMat]
	
					cv_model 	<- cv.glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1)	# Fit lasso regression model using k-fold cross-validation
					best_lambda <- cv_model$lambda.min
					best_model 	<- glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1, lambda = best_lambda)		# View coefficients of best model
					rL			<-	coef(best_model)
					MatrCc[iMat, col]	<- rL[2:NBN]		# Note that MatrCc is transposed
				}
			}	# Method LASSO
				
			if (Method == "LASSO1") {				# Digitalization = 1 if the result != 0	  -- Optimization of the hyper-parameters
				MatrCc[ , ] <- 0			
				for (iMat in 1:NBN) {
					col 	<- 1:NBN
					col 	<- col[-iMat]
	
					cv_model 	<- cv.glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1, intercept = FALSE, nfolds=3*NBN)	# Fit lasso regression model using k-fold cross-validation
					rL 			<- coef(cv_model, cv_model$lambda.1se)
					MatrCc[iMat, col]	<- rL[2:NBN]		# Note that MatrCc is transposed
				}
			}	# Method LASSO1

			if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
			#	Method "step" (automatic removal of coefficients by optimization of AIC : Akaike's Information Criterion)
			
				MatrCc[ , ] <- 0
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]	
					
					Donnees  <- data.frame(gY[ ,iMat], gX[ , ,iMat])
					#	cat("iMat ", iMat, " Data ", colnames(Donnees), "\n")
					Donnees	 <- rename(Donnees, "Y" ="gY...iMat.")		# The follow-up is clearer like this
					cc 		 <- colnames(Donnees)						# Name of the columns "Data". The first one is "Y"
					cc		 <- cc[-1]									# It remains the name of the coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
					colnames(rij) <- cc									# step gives the column names corresponding to the coefficients that are kept
							
					switch(Method, 
						"STEP-Fo" =
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},	# Forward
						"STEP-Ba" = 
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward
						"STEP-Bo" =
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)},		# Both
						)
					
					ll  	<- length(forward$coefficients)				# Number of coefficients kept by step
					nn		<- names(forward$coefficients)				# Name of these coefficients
			
					rij[1,]	<- 0
					if (ll >= 2) {
						for (i in 2:ll) {								# nn[1] = "(Intercept)"
							rij[1, nn[i]] <- forward$coefficients[nn[i]]
						}
					}

					#	cat (" iMat ", iMat, " ll ", ll, " nn ", nn, " rij ", rij, "\n")

					MatrCc[iMat, col]	<- rij[1, ]
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}		# loop iMat
			}	# STEP methods

			if (Method == "CLR") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))	
				MatrCc  <- clr(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# CLR
					
			if (Method == "ARACNE") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))
				MatrCc  <- aracne(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# ARACNE

			if (Method == "MRNET") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))
				MatrCc  <- mrnet(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# MRNET

			#	Looking for the solution
			#	Digitalization of the result		
			
			if (Method == "TLR") {
				MatrCcDig	<-	Digitalize (MatrCc, "Threshold", 0.25)
			}	else if (Method %in% c("MRA-KD", "MRA-KO")) {	
				Top = vTOP[iSize]
				MatrCcDig	<-	Digitalize (MatrCc, "Top", Top)	
			} else if (Method %in% c("LASSO", "LASSO1")) {
				MatrCcDig	<- Digitalize (MatrCc, "Threshold", 0)
			} else if (Method == "LSE_CI") {
				MatrCcDig	<-	MatrCc
			} else if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
				MatrCcDig	<- Digitalize (MatrCc, "Threshold", 0)
			} else if (Method %in% c("CLR", "ARACNE", "MRNET"))	{
				Top = vTOP[iSize]
				MatrCcDig	<-	Digitalize (MatrCc, "Top", Top)	
			}
			
			Score[ ,iMeth,iFic,iSize]	<- Benchmark(MatrCcDig,Solution)
		}		# Loop on iMeth
	}			# Loop on iFic
	
	for (iMeth in 1:NBM) {
		ScoreMoy["Se","mean",iMeth,iSize]	<-	mean(Score["Se",iMeth, ,iSize])
		ScoreMoy["Se","sd",iMeth,iSize]		<-	(var(Score["Se",iMeth, ,iSize])) ** 0.5
		ScoreMoy["Sp","mean",iMeth,iSize]	<-	mean(Score["Sp",iMeth, ,iSize])
		ScoreMoy["Sp","sd",iMeth,iSize]		<-	(var(Score["Sp",iMeth, ,iSize])) ** 0.5		
	}			# Loop on iMeth
}				# Loop on iSize
cat("END ", as.character(Sys.time()), "\n")


#
#	Figure 8 : Existing knowledge impact on regression methods		(ex Fig. 5)
#	Uses Dream Challenge 4 networks - see fig.4 for details about the differnt files.
#

NBFIC		<- 5								# There are 5 files for a given size
Sizes		<- c(10, 100)						# Nbr of nodes of the network
NBS			<- length(Sizes) 

vZF			<- list(from=c(0,0), to=c(100,60), by=c(10,15), Draws=c(100,10))		# Parameters (from, to, by, draws) of the sequences %ZF  -- Draws : number of draws by sequence of %ZF

Methods		<-  c("LASSO", "TLR", "STEP-Fo")
NBM			<- length(Methods)					# Number of méthods

ScoreTests	<-	c("nbZF", "TP", "TN", "FP", "FN", "Se", "Sp")		# Scores to measure (nbZF = number of forced zeros forced, Se = "Sensibility", Sp = "Specificity")
NBSC		<- length(ScoreTests)				# Number of scores

NBK			<- 2								# NBK corresponds to the number of perturbation types (KO or KD)

Moy				<- c("%ZF", "Mean", "sd")			# Mean values to display
matZFSens.df 	<- data.frame(0, 0, 0)				# Data Frame to use gplot
matZFSpec.df 	<- data.frame(0, 0, 0)				# Data Frame to use gplot


for (iSize in 1:NBS) {
	NBN			<- Sizes[iSize]						# Number of nodes
	NBD			<- vZF$Draws[iSize]					# Number of draws for each %ZF
	
	ZF  		<- seq(vZF$from[iSize], vZF$to[iSize], by=vZF$by[iSize])				# Percentage of forced zeros (%ZF)
	NBZF		<- length(ZF)						# Nbr. of percentages to test
	
	mZFSens			<- array(dim=c(length(Moy), NBZF, NBM))		# %ZF, mean and sd of sensibility
	mZFSpec			<- array(dim=c(length(Moy), NBZF, NBM))		# %ZF, mean and sd of specificity
	dimnames(mZFSens)	<- list(Moy, ZF, Methods)
	dimnames(mZFSpec)	<- list(Moy, ZF, Methods)	
		
	Solution	<- array(dim=c(NBN, NBN))			# Solution : the true network structure
	Solution2	<- array(dim=c(NBN, NBN))			# Idem "Solution", but row and column numbers correspond to numebers of "levels" ("G1", "G10", "G2", "G3" etc)
	Solution3	<- array(dim=c(NBN, NBN))			# Intermediary of calculation
	
	VZero		<- vector(length=NBN*NBN)			# The "true" zeros of the solution : (row number-1)*NBN + (column number-1)
	szVZero		<- 0								# Number of true zeros (excluding those on the diagonal)
	VZero0		<- vector(length=NBN*NBN)			# Copy of VZero at each simulation
	szVZero0	<- 0								# Copy of szVZero at each simulation	
	
	XMesNP		<- array(dim=c(NBN, NBK))			# Measured values, non pertubed system
	XMesP		<- array(dim=c(NBN, NBN, NBK))		# Measured values following a perturbation	
		
	MatR		<- array(dim=c(NBN, NBN, NBK))		# Matrix R
	pX   		<- array(dim=c(NBN-1, NBN-1))				# row (iRow), column (iCol)
	pY   		<- array(dim=c(NBN-1))						# row
	gX   		<- array(dim=c(NBK*(NBN-1), NBN-1, NBN))	# row (iRow), column (iCol), matrix number (iMat)
	gY   		<- array(dim=c(NBK*(NBN-1), NBN))			# row,  matrix number
	gXP	 		<- array(dim=c(NBK*(NBN-1), NBN-1, NBN, NBD, NBZF))	# gX modified thru the knowledge of some zeros
		
	Zero <- array(0, dim=c(NBN, NBN, NBD, NBZF))	# A 1 in position i,j says that the corresponding ri,j is null	
	
	MatrCc		<- array(dim=c(NBN, NBN))			# "r" matrix computed (according to method MRA, TLR, STEP or LASSO)
	MatrCcDig	<- array(dim=c(NBN, NBN))			# Computed Connection matrix digitalized (0, 1)	
	
	rij 		<- matrix(nrow=1, ncol=NBN-1)		# Intermediate calculation of rij
	rL 			<- array(dim=c(NBN))				# Result of Lasso method (ie sol. of Yi = Ai * Xi)
	
	Score 		<- array(0, dim=c(NBSC, NBD, NBZF, NBM, NBFIC, NBS))	# Nbr of zeros forced, TP, TN, FP, FN, Se, Sp
	dimnames(Score)		<- list(ScoreTests, NULL, ZF, Methods, NULL, Sizes)
	
	NomFic <- paste(vRoot, "ZeroForce.txt", sep="")
	
	for (iFic in 1:NBFIC) {
		valw <- paste("\n\n\n", "insilico_", NBN, "_", iFic, "\n", sep="")
		write(valw, NomFic, append=TRUE)
	
		#	1/ Data reading
		
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_wildtype.tsv", sep="")
		XMesNPLu 	<- fread(val, data.table=F)		# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_wildtype.tsv"
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockdowns.tsv", sep="")
		XMesPLu1 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_knockdowns.tsv"
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockouts.tsv", sep="")
		XMesPLu2 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1//insilico_size10_1_knockouts.tsv"		
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/DREAM4_GoldStandard_InSilico_Size", NBN, "_", iFic, ".tsv", sep="")
		SolutLu  	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
		SolutLu$V1	<- as.factor(SolutLu$V1)
		SolutLu$V2	<- as.factor(SolutLu$V2)

		for (i in 1:NBN) {
			XMesNP[i,1] = XMesNPLu[1,i]
			XMesNP[i,2] = XMesNPLu[1,i]
			
			for (j in 1:NBN) {
				XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation level -50%	
													# ATTENTION Dream Challenge data are transposed vs. my data 
				XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation level -100%
			}
		}
	
		Solution2[ , ] = 0
		for (i in 1:(NBN*(NBN-1))) {
			if(SolutLu[i, ]$V3 == 1) {
				Solution2[as.numeric(SolutLu[i, ]$V1), as.numeric(SolutLu[i, ]$V2)] = 1
			} else {
				break								# We suppose that all the "1" are situated at the beginning of the file
			}
		}
	
		numero <- as.numeric(sort(as.character(c(1:NBN))))
		for (i in 1:NBN) {
			Solution3[ ,numero[i]] = Solution2[ ,i]
		}
		for (i in 1:NBN) {
			Solution[ ,numero[i]] = Solution3[i,]	# We find the classical connection matrix (nodes "G1", "G2", ... , "G10")
													# Note that Solution is transposed
		}

		for (iPaq in c(1:NBK)) {					# We take every available data : Knock Down and Knock Out
			#	The "input data", delivered by Dream Challenge, have been loaded before in XMesNP and XMesP
			MatR[ , ,iPaq]  <- 2 * (XMesP[ , ,iPaq]-XMesNP[ ,iPaq]) / (XMesP[ , ,iPaq]+XMesNP[ ,iPaq])
			
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
				
				for (iRow in 1:(NBN-1)) {
					pY[iRow] = MatR[iMat, col[iRow], iPaq]
					
					for (iCol in 1:NBN-1) {
						pX[iRow, iCol] = MatR[col[iCol], col[iRow], iPaq]
					}
				}
			
				gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,iMat]  	<- pX[ , ]
				gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),iMat]   	<- pY[ ]
			}
		}
		
		#	2/ Looking for the position of "true" zeros
		
		VZero[ ] <- 0
		ind 	 <- 1		# Storage index in VZero
		for (i in 1:NBN) {
			for (j in 1:NBN) {
				if ((Solution[i,j] == 0) && (i != j)) {
					VZero[ind] = (i-1)*NBN + (j-1)
					ind  = ind +1
				}
			}
		}
		szVZero = ind -1		
		
		#	3/ Loop on the % of zeros and the simulations
	
		set.seed(12345)									# So as to generate the same sequences
		for (iZF in 1:NBZF) {							# Loop on the zeros to force				
			nbZF = round(ZF[iZF] * szVZero / 100)		# Nbr of zeros to force
			Name <- as.character(ZF[iZF])	
			
			valw <- paste("\n", "% forced zeros : ", ZF[iZF], " nbr. forced zeros : ", nbZF, sep="")
			write(valw, NomFic, append=TRUE)
		
			for (iDraw in 1:NBD) {							# NBD simulations
				nbZF0			<- nbZF			
				Zero[ , ,iDraw,iZF]	<- 0
				VZero0 			<- VZero
				szVZero0		<- szVZero
				
				#	3.1/ Choice of the zeros
				
				while (nbZF0 > 0) {
					indZF = round(runif(1, min=1, max=szVZero0))	# Position of THE real zero to take in the array
					val = VZero0[indZF]
					#	cat("nbZF0", nbZF0, " indZF ", indZF, " szVZero0 ", szVZero0, " val ", val, " Row ", trunc(val/NBN)+1, " Col ", val%%NBN+1, "\n")
					Zero[trunc(val/NBN)+1, val%%NBN+1, iDraw, iZF] <- 1
					
					for (i in indZF:(szVZero0-1)) {		# We remove the value used in VZero (draw without replacement)
						VZero0[i] = VZero0[i+1]
					}	# draw without replacement
					szVZero0 = szVZero0 - 1
					nbZF0 = nbZF0 - 1
				}		# while (nbZF0 > 0)
				
				#	3.2/ Taking into account the selected zeros
				
				gXP[ , , , iDraw, iZF]  <- gX[ , , ]
				for (Row in 1:NBN) {
					for (Col in 1:NBN) {
						if (Zero[Row, Col, iDraw, iZF] == 1) {			# We know that r[Row,Col] = 0
							if (Col > Row) {
								cc = Col -1
							} else {
								cc = Col
							}
							gXP[ ,cc,Row,iDraw,iZF] 	<- 1E6
							#	cat("gXP ", cc, Row, iDraw, iZF, "\n")
						}
					}	# for Col
				}		# for Row
		
				#	3.3/ Looking for MatrCc (depends on the chosen method)
				
				for	(iMeth in 1:NBM) {
					Method  <- Methods[iMeth]		
					# cat("File :", iFic, " %ZF ", ZF[iZF], " Test ", iDraw, " Method ", Method, as.character(Sys.time()), "\n")			
		
					if (Method == "TLR") {
						for (iMat in 1:NBN) {
							col <- 1:NBN
							col <- col[-iMat]

							Form <- paste("gY[ ,",iMat,"]~gXP[ , ,",iMat,",",iDraw,",",iZF,"]+0")					
							
							MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
							MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
						}
						MatrCc[is.na(MatrCc)]  <- 0
					}	# TLR method	
			
					if (Method == "LASSO") { 					#	LASSO's method, automatic choice of Lambda
						MatrCc[ , ] <- 0		# Initialization
						for (iMat in 1:NBN) {
							col <- 1:NBN
							col <- col[-iMat]
							
							exc  <- vector(length=0)						# column to exclude
							for (i in 1:NBN) {
								if (Zero[iMat, i, iDraw, iZF] == 1) {		# We know that r[iMat,i] = 0
									if (i > iMat) {
										cc = i -1
									} else {
										cc = i
									}
									exc <- cbind(exc, cc)
								}
							}	# for i					
					
							if (length(exc) < (NBN-1)) {					
								cv_model 	<- cv.glmnet(gXP[ , ,iMat,iDraw,iZF], gY[ ,iMat], alpha = 1, exclude = exc)			# Fit lasso regression model using k-fold cross-validation
								best_lambda <- cv_model$lambda.min
								best_model 	<- glmnet(gXP[ , ,iMat,iDraw,iZF], gY[ ,iMat], alpha = 1, lambda = best_lambda, exclude = exc)	# View coefficients of best model
								rL			<-	coef(best_model)								
				
								MatrCc[iMat, col]	<- rL[2:NBN]
							} else {
								MatrCc[iMat, col]	<- 0		# All the elements of row iMat are null
							}			
						}		
					}	# LASSO method			
			
					if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
					#	Method "step" (automatic removal of coefficients by optimization of AIC : Akaike's Information Criterion)
					
						MatrCc[ , ] <- 0
						for (iMat in 1:NBN) {
							col <- 1:NBN
							col <- col[-iMat]
					
							Donnees  <- data.frame(gY[ ,iMat], gXP[ , ,iMat,iDraw,iZF])
							#	cat("iMat ", iMat, " Data ", colnames(Donnees), "\n")
							Donnees	 <- rename(Donnees, "Y" ="gY...iMat.")		# The follow-up is clearer like this				
							cc 		 <- colnames(Donnees)						# Name of the columns "Data". The first one is "Y"				
							cc		 <- cc[-1]									# It remains the name of the coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)	
							colnames(rij) <- cc									# step gives the column names corresponding to the coefficients that are kept	
									
							switch(Method, 
								"STEP-Fo" =
								{intercept_only <- lm(Y ~ 1, data=Donnees)
									all 	<- lm(Y ~ ., data=Donnees)
									forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},		# Forward
								"STEP-Ba" = 
								{intercept_only <- lm(Y ~ 1, data=Donnees)
									all 	<- lm(Y ~ ., data=Donnees)
									forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward
								"STEP-Bo" =
								{intercept_only <- lm(Y ~ 1, data=Donnees)
									all 	<- lm(Y ~ ., data=Donnees)
									forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)},		# Both
								)
			
							ll  	<- length(forward$coefficients)				# Nbr of coefficients kept by step
							nn		<- names(forward$coefficients)				# Name of these coefficients
					
							rij[1, ]	<- 0
							if (ll >= 2) {
								for (i in 2:ll) {								# nn[1] = "(Intercept)"
									rij[1,nn[i]] <- forward$coefficients[nn[i]]
								}
							}
							#	cat (" iZF ", iZF, " iDraw ", iDraw, " iMat ", iMat, " ll ", ll, " nn ", nn, " rij ", rij, "\n")

							MatrCc[iMat, col]	<- rij[1, ]
							MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
						}		# loop iMat
					}	# STEP methods		

					#	3.4/ 	Looking for the solution
					#	Digitalization of the result
					
					if (Method == "TLR") {
						MatrCcDig	<-	Digitalize (MatrCc, "Threshold", 0.25)
					} else if (Method %in% c("LASSO", "LASSO1")) {
						MatrCcDig	<- Digitalize (MatrCc, "Threshold", 0)
					} else if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
						MatrCcDig	<- Digitalize (MatrCc, "Threshold", 0)
					}			
					
					#	3.5/	Test of the solution
									
					Score[1,iDraw,iZF,iMeth,iFic,iSize]			<- 	nbZF			
					Score[2:NBSC,iDraw,iZF,iMeth,iFic,iSize]	<- Benchmark(MatrCcDig,Solution)				
			
					matZFSens.df  	<- rbind(matZFSens.df, c(Score["Se",iDraw,iZF,iMeth,iFic,iSize], Method, Name))
					matZFSpec.df 	<- rbind(matZFSpec.df, c(Score["Sp",iDraw,iZF,iMeth,iFic,iSize], Method, Name))
	
					# cat("insilico_Size", NBN, "_", iFic, " Method ", Method, " Score ", Score[ ,iDraw,iZF,iMeth,iFic,iSize], "\n")
				
					valw <- paste(" Method : ", Method, " iDraw ", iDraw, " Score", paste(Score[ ,iDraw,iZF,iMeth,iFic,iSize], collapse=";"), sep=";")
					write(valw, NomFic, append=TRUE)
				}	# Loop on iMeth	(methods)
				
				write("\n", NomFic, append=TRUE)
			}		# Loop on iDraw	(100 if nb nodes = 10 or 10 trials if nb nodes = 100)
			
			cat("insilico_Size", NBN, "_", iFic, " Method ", Method, " %ZF ", nbZF, "\n")
		}			# Loop on iZF	(% of known zeros)
	}				# Loop on iFic	(5 files per size)	
}					# Loop on iSize (10 or 100 nodes)

cat("END ", as.character(Sys.time()), "\n")


#
#	Fig. 5a (iSize = 1)
#

matZFSens.df 	 		<- matZFSens.df[-1, , ]			# The first line, created during declaration , is (0,0,0,0)
names(matZFSens.df)		<- c("Sensib", "Meth", "ZF")
if (ZF[NBZF] == 100) {
	matZFSens.df$ZF		<- fct_relevel(matZFSens.df$ZF, "100", after=Inf)	#  100% is put at the end (R follows by defaut the lexicographic order)
}
matZFSens.df$Sensib		<- as.numeric(matZFSens.df$Sensib)
matZFSens.df$Meth		<- as.factor (matZFSens.df$Meth)
matZFSens.df$ZF			<- as.factor (matZFSens.df$ZF)

NomFic 	<- paste(vRoot, "Fig5a_Sens10.pdf", sep="")
#NomFic 	<- paste(vRoot, "Fig5c_Sens100.pdf", sep="")
pdf(file=NomFic, height=4, width=4) 					# Sizes default to 7

ggplot(matZFSens.df, aes(x=ZF, y=Sensib, fill=Meth), ylimit=c(0.35,1)) + 
	geom_boxplot(outlier.colour="black", outlier.size=0.3) + 
	labs(title="Sensibility vs. % values forced to 0", x="% values forced to 0", y="Sensibility") +
	theme(legend.position = "none") +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black"))	
#    scale_colour_manual(name = "Meth",
#		values = c("TLR" = "green", "LASSO" = "red", "STEP-Fo" = "orange"),
#		labels = c("TLR" = "TLR", "LASSO" = "LASSO", "STEP-Fo" = "STEP-Fo"))
		
dev.off()


#
#	Fig. 5b (iSize = 1)
#

matZFSpec.df 	 		<- matZFSpec.df[-1, , ]			# The first line, created during declaration , is (0,0, 0)
names(matZFSpec.df)		<- c("Specif", "Meth", "ZF")
if (ZF[NBZF] == 100) {
	matZFSpec.df$ZF		<- fct_relevel(matZFSpec.df$ZF, "100", after=Inf)	#  100% is put at the end (R follows by defaut the lexicographic order)
}
matZFSpec.df$Specif		<- as.numeric(matZFSpec.df$Specif)
matZFSpec.df$Meth		<- as.factor (matZFSpec.df$Meth)
matZFSpec.df$ZF			<- as.factor (matZFSpec.df$ZF)

NomFic 	<- paste(vRoot, "Fig5b_Spec10.pdf", sep="")
#NomFic 	<- paste(vRoot, "Fig5d_Spec100.pdf", sep="")
pdf(file=NomFic, height=4, width=4) 					# Sizes default to 7

ggplot(matZFSpec.df, aes(x=ZF, y=Specif, fill=Meth)) + 
	geom_boxplot(outlier.colour="black", outlier.size=0.3) + 
	labs(title="Specificity vs. % values forced to 0", x="% values forced to 0", y="Specificity") +
	theme(legend.position = "none")	+
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black"))
	
#    scale_colour_manual(name = "Meth",
#		values = c("TLR" = "green", "LASSO" = "red", "STEP-Fo" = "orange"),
#		labels = c("TLR" = "TLR", "LASSO" = "LASSO", "STEP-Fo" = "STEP-Fo"))
		
dev.off()


#
#	Fig. 5c and 5d
#	Idem fig 5a and 5b. Corresponds to iSize =2.
#	Replace the name of the files with  "Fig5c_Sens100.pdf" and "Fig5d_Spec100.pdf"
#


#
#	Supp. Info
#	1/ A 6 node model of MAP kinases
#
#	Rate Expressions Used for the Models of MAPK Cascade (Kinase network , figures 1 to 3) 
#	These data come from the article : 
#	"Untangling the wires: A strategy to trace functional interactions in signaling and gene networks". Kholodenko et al, see ref.
#
#
#	2/ Supplementary Tables & Figures
#
#	Supplementary table 1  : see fig. 3c
#

#
#	Supplementary Figure 1
#
#	ROC curves for DREAM 4 Challenge, corresponding to challenges "insilico_size10_1", "**_size100_1" and many methods. 
#	Update the variables "iSize" and "iFic" to draw the corresponding pictures.
#	See fig. 4 for details about the Dream Challenge files.
#

NBFIC		<- 5									# There are 5 files for a given size
Sizes		<- c(10, 100)							# Nbr of nodes of the network
NBS			<- length(Sizes) 

Methods		<-  c("MRA-KD", "MRA-KO", "LSE_CI", "LASSO", "TLR", "STEP-Fo", "STEP-Ba", "STEP-Bo", "CLR", "ARACNE", "MRNET")
NBM			<- length(Methods)						# Number of méthods

vTOP		<- c(0.2, 0.25)							# "top" of the edges to take into account
vPValue		<- 0.05									# pValue expected must be < vPValue

NBK			<- 2									# NBK corresponds to the number of perturbation types (KO or KD)

Seuils		<- seq(0, 1, by=0.05)					# To plot the curves ROC (KO, KD, TLR)
													# For LASSO's method, we will multiply every value by KLASSO
NBSeuils	<- length(Seuils)
KLASSO		<- 0.05									# arbitrary choice, to get many dota, may be modified if needed
KSTEP		<- 3									# arbitrary choice, to get many dota, may be modified if needed

ScoreTests	<-	c("TP", "TN", "FP", "FN", "Se", "Sp")		# Scores to measure (Se = "Sensibility", Sp = "Specificity")
NBSC		<- length(ScoreTests)					# Number of scores
Score 		<- array(0, dim=c(NBSC, NBSeuils, NBM, NBFIC, NBS))		# TP, TN, FP, FN, Se, Sp
dimnames(Score)		<- list(ScoreTests, Seuils, Methods, NULL, Sizes)

Pts			<- array(0, dim=c(2, NBM, NBFIC, NBS))			# Spec and Sensib corresponding to the points to draw on the ROC curves (not plotted methods)
dimnames(Pts)		<- list(c("Se", "Sp"), Methods, NULL, Sizes)
AuROC		<- array(0, dim=c(NBM, NBFIC, NBS))				# AUROC corresponding to the plotted methods, computed with roc.area
dimnames(AuROC)		<- list(Methods, NULL, Sizes)
AuROCZ		<- array(0, dim=c(NBM, NBFIC, NBS))				# AUROC corresponding to the plotted methods, computed with trapz
dimnames(AuROCZ)	<- list(Methods, NULL, Sizes)
PValue		<- array(0, dim=c(NBM, NBFIC, NBS))				# PValue corresponding to the plotted methods (vs. random selection)
dimnames(PValue)	<- list(Methods, NULL, Sizes)

		
for (iSize in 1:NBS) {
	NBN			<- Sizes[iSize]
	
	Solution	<- array(dim=c(NBN, NBN))			# Solution : the true network structure
	Solution2	<- array(dim=c(NBN, NBN))			# Idem "Solution", but row and column numbers correspond to numebers of "levels" ("G1", "G10", "G2", "G3" etc)
	Solution3	<- array(dim=c(NBN, NBN))			# Intermediary of calculation
	
	XMesNP		<- array(dim=c(NBN, NBK))			# Measured values, non pertubed system
	XMesP		<- array(dim=c(NBN, NBN, NBK))		# Measured values following a perturbation	
		
	MatR		<- array(dim=c(NBN, NBN, NBK))		# Matrix R
	MatRm1		<- array(dim=c(NBN, NBN, NBK))		# Matrix R	** -1
	dgMatRm1M1	<- array(dim=c(NBN, NBN, NBK))		# Diagonal (R** -1) **-1
	pX   		<- array(dim=c(NBN-1, NBN-1))				# row (iRow), column (iCol)
	pY   		<- array(dim=c(NBN-1))						# row
	gX   		<- array(dim=c(NBK*(NBN-1), NBN-1, NBN))	# row (iRow), column (iCol), matrix number (iMat)
	gY   		<- array(dim=c(NBK*(NBN-1), NBN))			# row,  matrix number
	gXO 		<- array(dim=c(NBN-1, NBN-1, NBN))			# gX matrix truncated so as to keep the values KO only
	gYO 		<- array(dim=c(NBN-1, NBN))					# gY matrix truncated so as to keep the values KO only
	gXD 		<- array(dim=c(NBN-1, NBN-1, NBN))			# gX matrix truncated so as to keep the values KD only
	gYD 		<- array(dim=c(NBN-1, NBN))					# gY matrix truncated so as to keep the values KD only
		
	MatrCc		<- array(dim=c(NBN, NBN))			# "r" matrix computed (according to method MRA, TLR, STEP or LASSO)
	MatrPVal	<- array(dim=c(NBN, NBN))			# p-value
	MatrCcDig	<- array(dim=c(NBN, NBN))			# Computed Connection matrix digitalized (0, 1)
	
	vMatrCc  	<- vector(length = NBN*NBN)			# MatrCc as a vector (to use roc.area)
	vSolution	<- vector(length = NBN*NBN)			# Solution as a vector (to use roc.area)	
	
	best_lambda	<- array(dim=c(NBN))				# Lasso hyper parameter
	rij 		<- matrix(nrow=1, ncol=NBN-1)		# Intermediate calculation of rij
	rL 			<- array(dim=c(NBN))				# Result of Lasso method (ie sol. of Yi = Ai * Xi)
	
	for (iFic in 1:NBFIC) {
		#	1/ Data reading
		
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_wildtype.tsv", sep="")
		XMesNPLu 	<- fread(val, data.table=F)		# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_wildtype.tsv"
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockdowns.tsv", sep="")
		XMesPLu1 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_knockdowns.tsv"
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockouts.tsv", sep="")
		XMesPLu2 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1//insilico_size10_1_knockouts.tsv"		
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/DREAM4_GoldStandard_InSilico_Size", NBN, "_", iFic, ".tsv", sep="")
		SolutLu  	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
		SolutLu$V1	<- as.factor(SolutLu$V1)
		SolutLu$V2	<- as.factor(SolutLu$V2)

		for (i in 1:NBN) {
			XMesNP[i,1] = XMesNPLu[1,i]
			XMesNP[i,2] = XMesNPLu[1,i]
			
			for (j in 1:NBN) {
				XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation level -50%	
													# ATTENTION Dream Challenge data are transposed vs. my data 
				XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation level -100%
			}
		}
	
		Solution2[,] = 0
		for (i in 1:(NBN*(NBN-1))) {
			if(SolutLu[i, ]$V3 == 1) {
				Solution2[as.numeric(SolutLu[i, ]$V1), as.numeric(SolutLu[i, ]$V2)] = 1
			} else {
				break								# We suppose that all the "1" are situated at the beginning of the file
			}
		}
	
		numero <- as.numeric(sort(as.character(c(1:NBN))))
		for (i in 1:NBN) {
			Solution3[ ,numero[i]] = Solution2[ ,i]
		}
		for (i in 1:NBN) {
			Solution[ ,numero[i]] = Solution3[i,]	# We find the classical connection matrix (nodes "G1", "G2", ... , "G10")
													# Note that Solution is transposed
		}

		for (iPaq in c(1:NBK)) {					# We take every available data : Knock Down and Knock Out
			#	The "input data", delivered by Dream Challenge, have been loaded before in XMesNP and XMesP
			MatR[ , ,iPaq]  <- 2 * (XMesP[ , ,iPaq]-XMesNP[ ,iPaq]) / (XMesP[ , ,iPaq]+XMesNP[ ,iPaq])
			
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
				
				for (iRow in 1:(NBN-1)) {
					pY[iRow] = MatR[iMat, col[iRow], iPaq]
					
					for (iCol in 1:NBN-1) {
						pX[iRow, iCol] = MatR[col[iCol], col[iRow], iPaq]
					}
				}
			
				gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,iMat]  	<- pX[ , ]
				gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),iMat]   	<- pY[ ]
			}
		}
	
		gXD 	<- gX[1 : (NBN-1), , ] 					# KnockDown
		gYD		<- gY[1 : (NBN-1), ]	
		gXO		<- gX[NBN : ((NBN-1)*2), , ] 			# KnockOut
		gYO		<- gY[NBN : ((NBN-1)*2), ] 
		
		for	(iMeth in 1:NBM) {
			Method  <- Methods[iMeth]
		
			cat("iSize ", iSize, " iFic ", iFic, " Method ", Method, " Heure ", as.character(Sys.time()), "\n")
			
			if (Method %in% c("MRA-KD", "MRA-KO")) {	#	Classical MRA method
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
					
					if (Method == "MRA-KD") {
						Form <- paste("gYD[ ,",iMat,"]~gXD[ , ,",iMat,"]+0")
					} else {
						Form <- paste("gYO[ ,",iMat,"]~gXO[ , ,",iMat,"]+0")
					}	
					
					MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}
				MatrCc[is.na(MatrCc)]  <- 0
			}	# MRA methods

			if (Method == "LSE_CI") {
				MatrCc[ , ] <- 0
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
			
					Form <- paste("gY[ ,",iMat,"]~gX[ , ,",iMat,"]+0")
					pQ 	 <- lm(as.formula(Form))								# MatrCc contient les p-values associées
					MatrPVal[iMat,col] 	<- broom::tidy(pQ)$p.value
				}
				MatrCc	<- MatrPVal
				diag(MatrCc)	<- 0
				MatrCc[is.na(MatrCc)]  <- 0	
			}	# LSE_CI method	
	
			if (Method == "TLR") {
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
			
					Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")
					MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}
				MatrCc[is.na(MatrCc)]  <- 0	
			}	# TLR with Linear regression, 2 measures per node		
			
#			if (Method == "LASSO") {				# Digitalization = sign of the result	-- 3 replicates  -- automatic choice of lambda
#				for (iMat in 1:NBN) {
#					col 	<- 1:NBN
#					col 	<- col[-iMat]
#	
#					cv_model 	<- cv.glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1)	# Fit lasso regression model using k-fold cross-validation
#					best_lambda [iMat] <- cv_model$lambda.min	
#				}
#			}
					
###				
###			if (Method == "LASSO1") {					# Digitalization = 1 if the result != 0	  -- Optimization of the hyper-parameters
###				MatrCc[ , ] <- 0			
###				for (iMat in 1:NBN) {
###					col 	<- 1:NBN
###					col 	<- col[-iMat]
###	
###					cv_model 	<- cv.glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1, intercept = FALSE, nfolds=3*NBN)	# Fit lasso regression model using k-fold cross-validation
###					rL 			<- coef(cv_model, cv_model$lambda.1se)
###					MatrCc[iMat, col]	<- rL[2:NBN]	# Note that MatrCc is transposed
###				}
###			}	# Method LASSO1

			if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
			#	Method "step" (automatic removal of coefficients by optimization of AIC : Akaike's Information Criterion)
			
				MatrCc[ , ] <- 0
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]	
					
					Donnees  <- data.frame(gY[ ,iMat], gX[ , ,iMat])
					#	cat("iMat ", iMat, " Data ", colnames(Donnees), "\n")
					Donnees	 <- rename(Donnees, "Y" ="gY...iMat.")		# The follow-up is clearer like this
					cc 		 <- colnames(Donnees)						# Name of the columns "Data". The first one is "Y"
					cc		 <- cc[-1]									# It remains the name of the coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
					colnames(rij) <- cc									# step gives the column names corresponding to the coefficients that are kept
							
					switch(Method, 
						"STEP-Fo" =
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},	# Forward
						"STEP-Ba" = 
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward
						"STEP-Bo" =
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)},		# Both
						)
					
					ll  	<- length(forward$coefficients)				# Number of coefficients kept by step
					nn		<- names(forward$coefficients)				# Name of these coefficients
			
					rij[1, ]	<- 0
					if (ll >= 2) {
						for (i in 2:ll) {								# nn[1] = "(Intercept)"
							rij[1, nn[i]] <- forward$coefficients[nn[i]]
						}
					}

					#	cat (" iMat ", iMat, " ll ", ll, " nn ", nn, " rij ", rij, "\n")

					MatrCc[iMat, col]	<- rij[1, ]
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}		# loop iMat
			}	# STEP methods

			if (Method == "CLR") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))	
				MatrCc  <- clr(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# CLR
					
			if (Method == "ARACNE") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))
				MatrCc  <- aracne(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# ARACNE

			if (Method == "MRNET") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))
				MatrCc  <- mrnet(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# MRNET

			#	Digitalization of the result and scoring	
			
			for (iSeuil in (1:NBSeuils)) {
				Seuil	<- Seuils[iSeuil]
				
				if (Method %in% c("MRA-KD", "MRA-KO")) {
					MatrCcDig	<-	Digitalize (MatrCc, "Top", Seuil)
				} else if (Method == "TLR") {
					MatrCcDig	<-	Digitalize (MatrCc, "Threshold", Seuil)		
				} else if (Method %in% c("LASSO", "LASSO1")) {
					BLambda = Seuil*KLASSO
					MatrCc[ , ] <- 0
					for (iMat in 1:NBN) {
						col 	<- 1:NBN
						col 	<- col[-iMat]

						best_model 	<- glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1, lambda = BLambda)		# View coefficients of best model
						rL			<- coef(best_model)
						MatrCc[iMat, iMat] 	<- -1						
						MatrCc[iMat, col]	<- rL[2:NBN]		# Note that MatrCc is transposed
					}
					MatrCcDig	<- Digitalize (MatrCc, "Threshold", 0)					
				} else if (Method == "LSE_CI") {
					MatrCcDig	<-	Digitalize (MatrCc, "IThreshold", Seuil)
				} else if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
					MatrCcDig	<- Digitalize (MatrCc, "Threshold", Seuil*KSTEP)
				} else if (Method %in% c("CLR", "ARACNE", "MRNET")) {
					MatrCcDig	<- Digitalize (MatrCc, "Threshold", Seuil)
				}			
				
				Score[ ,iSeuil,iMeth,iFic,iSize]	<- Benchmark(MatrCcDig,Solution)
			}	# Loop on iSeuil	
				
			vMatrCc[]  		<- as.numeric(abs(MatrCc[]))
			vSolution[]		<- as.numeric(Solution[])
			rb 	<- roc.area(vSolution, vMatrCc)
			AuROC [iMeth, iFic, iSize]	 = rb$A	
			PValue[iMeth, iFic, iSize]	 = rb$p.value
														# Use this method whenever possible. It's more accurate
			AuROCZ[iMeth,iFic,iSize]	 =  abs(trapz(1-Score["Sp", ,iMeth,iFic,iSize], Score["Se", ,iMeth,iFic,iSize]))
														# When roc.area doesn't work : bad p-value ...
			
			AuROC [iMeth, iFic, iSize]	 = max(AuROC [iMeth, iFic, iSize], AuROCZ [iMeth, iFic, iSize])
														# AUC ROC for the plotted methods				
		}		# Loop on iMeth
		
###			for (iMat in 1:NBN) {
###				col 	<- 1:NBN
###				col 	<- col[-iMat]
###			
###				cv_model 	<- cv.glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1)	# Fit lasso regression model using k-fold cross-validation
###				best_lambda <- cv_model$lambda.min			
###				best_model 	<- glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1, lambda = best_lambda)	# View coefficients of best model
###				rL			<-	coef(best_model)
###				MatrCc[iMat, col]	<- rL[2:NBN]			# Note that MatrCc is transposed
###			}			
###			MatrCcDig	<- Digitalize (MatrCc, "Threshold", 0)
###			A = Benchmark(MatrCcDig,Solution)
###			Pts["Sp","LASSO",iFic,iSize] = A[6]				# Best Lambda for LASSO's method  -- Sp
###			Pts["Se","LASSO",iFic,iSize] = A[5]				# Best Lambda for LASSO's method  -- Se
	}			# Loop on iFic
}				# Loop on iSize
cat("END ", as.character(Sys.time()), "\n")	


#
#	Answers to Reviewers
#	§ 1.1.1 2nd table : comparison of the methods (AUROC, pValue)
#	New Figure 4c
#

for (iMeth in 1:NBM) {
	cat("Method ", Methods[iMeth], "Size 10 :  AUROC ", mean(AuROC [iMeth, ,1]), " pValue ", mean(PValue[iMeth, ,1]), " AUROCZ ", mean(AuROCZ [iMeth, ,1]), "\n")
	cat("Method ", Methods[iMeth], "Size 100 : AUROC ", mean(AuROC [iMeth, ,2]), " pValue ", mean(PValue[iMeth, ,2]), " AUROCZ ", mean(AuROCZ [iMeth, ,2]), "\n")
}		# Loop on iMeth				


ScoreM <- array(0, dim=c(2,NBSeuils+2,NBM,1,NBS))
dimnames(ScoreM)	<- list(c("Se","Sp"), c(0,Seuils,1), Methods, NULL, Sizes)
ScoreM[ ,2:(NBSeuils+1), ,1, ]  <- Score[c("Se","Sp"), , ,1, ]

for (iMeth in 1:NBM) {
	Method	<- Methods[iMeth]
	if (Method %in% c("MRA-KD", "MRA-KO", "LSE_CI")) {
		ScoreM["Se",1,iMeth,1, ]  			<- 0
		ScoreM["Se",NBSeuils+2,iMeth,1, ]	<- 1
		ScoreM["Sp",1,iMeth,1, ]  			<- 1
		ScoreM["Sp",NBSeuils+2,iMeth,1, ]	<- 0
	} else {
		ScoreM["Se",1,iMeth,1, ]  			<- 1
		ScoreM["Se",NBSeuils+2,iMeth,1, ]	<- 0
		ScoreM["Sp",1,iMeth,1, ]  			<- 0
		ScoreM["Sp",NBSeuils+2,iMeth,1, ]	<- 1
	} 
}


#
#	Supp. Info. Fig 1 a
#
NomFic 	<- paste(vRoot, "Fig_SI1_ROC10.pdf", sep="")
pdf(file=NomFic, height=9, width=9) 					# Sizes default to 7
plot (1-ScoreM["Sp", ,"MRA-KD",1,1], ScoreM["Se", ,"MRA-KD",1,1], 
	main= "ROC : insilico_size10_1", type="l", col="blue", 
	xlab=paste("1 - Specificity\nAU ROC  KD = ", round(AuROC["MRA-KD",1,1], digits=3), " KO = ", round(AuROC["MRA-KO",1,1], digits=3),
			   " TLR = ", round(AuROC["TLR",1,1], digits=3), " LASSO = ", round(AuROC["LASSO",1,1], digits=3), "\n",
			   "STEP-Fo = ", round(AuROC["STEP-Fo",1,1], digits=3), "STEP-Ba = ", round(AuROC["STEP-Ba",1,1], digits=3), "STEP-Bo = ", round(AuROC["STEP-Bo",1,1], digits=3)),
	ylab="Sensibility", xlim=c(0,1), ylim=c(0,1))
#			   " TLR = ", round(AuROC["TLR",1,1], digits=3), " LASSO = ", round(AuROC["LASSO",1,1], digits=3), "LSE_CI =", round(AuROC["LSE_CI",1,1], digits=3), "\n",	
lines(1-ScoreM["Sp", ,"MRA-KO",1,1],	ScoreM["Se", ,"MRA-KO",1,1],	col="black", type="l")
lines(1-ScoreM["Sp", ,"TLR",1,1],   ScoreM["Se", ,"TLR",1,1],   col="red", 	 type="l")
lines(1-ScoreM["Sp", ,"LASSO",1,1], ScoreM["Se", ,"LASSO",1,1], col="green",	 type="l")
#lines(1-ScoreM["Sp", ,"LSE_CI",1,1], ScoreM["Se", ,"LSE_CI",1,1], col="#CCCCCC",	 type="l")
lines(1-ScoreM["Sp", ,"STEP-Fo",1,1], ScoreM["Se", ,"STEP-Fo",1,1], col="#CC3FFF",	 type="l")
lines(1-ScoreM["Sp", ,"STEP-Ba",1,1], ScoreM["Se", ,"STEP-Ba",1,1], col="#FF9933",	 type="l")
lines(1-ScoreM["Sp", ,"STEP-Bo",1,1], ScoreM["Se", ,"STEP-Bo",1,1], col="#FFFF00",	 type="l")
abline(c(0,0), c(1,1), col="grey", lty="dashed")
legend("bottomright", c("MRA-KD","MRA-KO","TLR","LASSO","STEP-Fo","STEP-Ba","STEP-Bo"), fill=c("blue","black","red","green","orange","lightpink3","#FFFF00"))
#legend("bottomright", c("MRA-KD","MRA-KO","TLR","LASSO","LSE_CI","STEP-Fo","STEP-Ba","STEP-Bo"), fill=c("blue","black","red","green","#CCCCCC","#CC3FFF","#FF9933","#FFFF00"))
dev.off()

# points(1-Pts["Sp","LSE_CI",1,1],	Pts["Se","LSE_CI",1,1],		pch=19, col="maroon4", 	bg="maroon4")
# points(1-Pts["Sp","STEP-Fo",1,1], 	Pts["Se","STEP-Fo",1,1], 	pch=22, col="maroon4", 	bg="maroon4")
# points(1-Pts["Sp","STEP-Ba",1,1], 	Pts["Se","STEP-Ba",1,1], 	pch=24, col="maroon4", 	bg="maroon4")
# points(1-Pts["Sp","STEP-Bo",1,1], 	Pts["Se","STEP-Bo",1,1], 	pch=25, col="maroon4", 	bg="maroon4")
# points(1-Pts["Sp","LASSO",1,1],   	Pts["Se","LASSO",1,1],   	pch=23, col="orange", 	bg="orange")
# legend("bottom", c("LSE_CI","STEP-Fo","STEP-Ba","STEP-Bo","LASSO"), pch=c(19,22,24,25,23),
#	   col=c("maroon4","maroon4","maroon4","maroon4","maroon4"))
   
#
#	Supp. Info. Fig 1 b
#
NomFic 	<- paste(vRoot, "Fig_SI1_ROC100.pdf", sep="")
pdf(file=NomFic, height=9, width=9) 					# Sizes default to 7
plot (1-ScoreM["Sp", ,"MRA-KD",1,2], ScoreM["Se", ,"MRA-KD",1,2], 
	main= "ROC : insilico_size100_1", type="l", col="blue", 
	xlab=paste("1 - Specificity\nAU ROC  KD = ", round(AuROC["MRA-KD",1,2], digits=3), " KO = ", round(AuROC["MRA-KO",1,2], digits=3),
			   " TLR = ", round(AuROC["TLR",1,2], digits=3), " LASSO = ", round(AuROC["LASSO",1,2], digits=3), "\n",
			   "STEP-Fo = ", round(AuROC["STEP-Fo",1,2], digits=3), "STEP-Ba = ", round(AuROC["STEP-Ba",1,2], digits=3), "STEP-Bo = ", round(AuROC["STEP-Bo",1,2], digits=3)), 
	ylab="Sensibility", xlim=c(0,1), ylim=c(0,1))
#			   " TLR = ", round(AuROC["TLR",1,2], digits=3), " LASSO = ", round(AuROC["LASSO",1,2], digits=3), "LSE_CI =", round(AuROC["LSE_CI",1,2], digits=3), "\n",
lines(1-ScoreM["Sp", ,"MRA-KO",1,2],	ScoreM["Se", ,"MRA-KO",1,2],	col="black", type="l")
lines(1-ScoreM["Sp", ,"TLR",1,2],   ScoreM["Se", ,"TLR",1,2],   col="red", 	 type="l")
lines(1-ScoreM["Sp", ,"LASSO",1,2], ScoreM["Se", ,"LASSO",1,2], col="green",	 type="l")
#lines(1-ScoreM["Sp", ,"LSE_CI",1,2], ScoreM["Se", ,"LSE_CI",1,2], col="#CCCCCC",	 type="l")
lines(1-ScoreM["Sp", ,"STEP-Fo",1,2], ScoreM["Se", ,"STEP-Fo",1,2], col="#CC3FFF",	 type="l")
lines(1-ScoreM["Sp", ,"STEP-Ba",1,2], ScoreM["Se", ,"STEP-Ba",1,2], col="#FF9933",	 type="l")
lines(1-ScoreM["Sp", ,"STEP-Bo",1,2], ScoreM["Se", ,"STEP-Bo",1,2], col="#FFFF00",	 type="l")
abline(c(0,0), c(1,1), col="grey", lty="dashed")
legend("bottomright", c("MRA-KD","MRA-KO","TLR","LASSO","STEP-Fo","STEP-Ba","STEP-Bo"), fill=c("blue","black","red","green","orange","lightpink3","#FFFF00"))
#legend("bottomright", c("MRA-KD","MRA-KO","TLR","LASSO","LSE_CI","STEP-Fo","STEP-Ba","STEP-Bo"), fill=c("blue","black","red","green","#CCCCCC","#CC3FFF","#FF9933","#FFFF00"))
dev.off()
	   
# points(1-Pts["Sp","LSE_CI",1,2],	Pts["Se","LSE_CI",1,2],		pch=19, col="maroon4", 	bg="maroon4")
# points(1-Pts["Sp","STEP-Fo",1,2], 	Pts["Se","STEP-Fo",1,2], 	pch=22, col="maroon4", 	bg="maroon4")
# points(1-Pts["Sp","STEP-Ba",1,2], 	Pts["Se","STEP-Ba",1,2], 	pch=24, col="maroon4", 	bg="maroon4")
# points(1-Pts["Sp","STEP-Bo",1,2], 	Pts["Se","STEP-Bo",1,2], 	pch=25, col="maroon4", 	bg="maroon4")
# points(1-Pts["Sp","LASSO",1,2],   	Pts["Se","LASSO",1,2],   	pch=23, col="orange", 	bg="orange")
# legend("bottom", c("LSE_CI","STEP-Fo","STEP-Ba","STEP-Bo","LASSO"), pch=c(19,22,24,25,23),
# 	   col=c("maroon4","maroon4","maroon4","maroon4","maroon4"))

#
#	Supplementary Figure 2, Figure 3, Table 2, Table 3
#
#	ROC curves for DREAM 4 Challenge, corresponding to challenges "insilico_size10_1"... "**_size10_5" (fig 2), "insilico_size100_1"... "**_size100_5" (fig 3) and many methods. 
#	Update the variables "iSize" and "iFic" to scan the 5 files and draw the corresponding picture.
#	See fig. 4 for details about the Dream Challenge files.
#

NBFIC		<- 5									# There are 5 files for a given size
Sizes		<- c(10, 100)							# Nbr of nodes of the network
NBS			<- length(Sizes) 

Methods		<-  c("MRA-KD", "MRA-KO", "TLR", "LASSO", "STEP-Fo", "STEP-Ba", "STEP-Bo", "CLR", "ARACNE", "MRNET")
NBM			<- length(Methods)						# Number of méthods

vTOP		<- c(0.2, 0.25)							# "top" of the edges to take into account

NBK			<- 2									# NBK corresponds to the number of perturbation types (KO or KD)

Seuils		<- seq(0, 1, by=0.05)					# To plot the curves ROC (KO, KD, TLR)
													# For LASSO's method, we will multiply every value by KLASSO
NBSeuils	<- length(Seuils)

ScoreTests	<-	c("TP", "TN", "FP", "FN", "Se", "Sp")		# Scores to measure (Se = "Sensibility", Sp = "Specificity")
NBSC		<-  length(ScoreTests)					# Number of scores
Score 		<-  array(0, dim=c(NBSC, NBSeuils, NBM, NBFIC, NBS))		# TP, TN, FP, FN, Se, Sp
dimnames(Score)		<- list(ScoreTests, Seuils, Methods, NULL, Sizes)

ScoreDCTests<-	c("AUC", "AUPR", "PVAUC", "PVAUPR")			# Scores used by DreamChallenge
NBSCDC		<-  length(ScoreDCTests)						# Number of scores
ScoreDC 			<- array(0, dim=c(NBSCDC, NBM, NBFIC, NBS))	# "AUC", "AUPR", "PVAUC", "PVAUPR" 1° dim.
dimnames(ScoreDC)	<- list(ScoreDCTests, Methods, NULL, Sizes)

ScoreDCOut 			<- array(0, dim=c(NBM+1, NBFIC+3, NBS))	# Title and methods 1° dim.  - Out value (AUC) for the NBFIC files + mean, sd  and AUROC Score 2° dim
dimnames(ScoreDCOut)<- list(c("AUC",Methods), NULL, Sizes)

Pts			<- array(0, dim=c(2, NBM, NBFIC, NBS))			# Spec and Sensib corresponding to the points to draw on the ROC curves (not plotted methods)
dimnames(Pts)		<- list(c("Se", "Sp"), Methods, NULL, Sizes)
AuROC		<- array(0, dim=c(NBM, NBFIC, NBS))				# AUROC corresponding to the plotted methods
dimnames(AuROC)		<- list(Methods, NULL, Sizes)
		
cat("START ", as.character(Sys.time()), "\n")	
		
for (iSize in 1:NBS) {
	NBN			<- Sizes[iSize]
	
	Solution	<- array(dim=c(NBN, NBN))			# Solution : the true network structure
	Solution2	<- array(dim=c(NBN, NBN))			# Idem "Solution", but row and column numbers correspond to numebers of "levels" ("G1", "G10", "G2", "G3" etc)
	Solution3	<- array(dim=c(NBN, NBN))			# Intermediary of calculation
	
	XMesNP		<- array(dim=c(NBN, NBK))			# Measured values, non pertubed system
	XMesP		<- array(dim=c(NBN, NBN, NBK))		# Measured values following a perturbation	
		
	MatR		<- array(dim=c(NBN, NBN, NBK))		# Matrix R
	pX   		<- array(dim=c(NBN-1, NBN-1))				# row (iRow), column (iCol)
	pY   		<- array(dim=c(NBN-1))						# row
	gX   		<- array(dim=c(NBK*(NBN-1), NBN-1, NBN))	# row (iRow), column (iCol), matrix number (iMat)
	gY   		<- array(dim=c(NBK*(NBN-1), NBN))			# row,  matrix number
	gXKO 		<- array(dim=c(NBN-1, NBN-1, NBN))			# gX matrix truncated so as to keep the values KO only
	gYKO 		<- array(dim=c(NBN-1, NBN))					# gY matrix truncated so as to keep the values KO only
	gXKD 		<- array(dim=c(NBN-1, NBN-1, NBN))			# gX matrix truncated so as to keep the values KD only
	gYKD 		<- array(dim=c(NBN-1, NBN))					# gY matrix truncated so as to keep the values KD only
		
	MatrCc		<- array(dim=c(NBN, NBN))			# "r" matrix computed (according to method MRA, TLR, STEP or LASSO)
	MatrCcDig	<- array(dim=c(NBN, NBN))			# Computed Connection matrix digitalized (0, 1)
	rij 		<- matrix(nrow=1, ncol=NBN-1)		# Intermediate calculation of rij
	rL 			<- array(dim=c(NBN))				# Result of Lasso method (ie sol. of Yi = Ai * Xi)	
	
	vMatrCc  	<- vector(length = NBN*NBN)			# MatrCc as a vector (to use roc.area)
	vSolution	<- vector(length = NBN*NBN)			# Solution as a vector (to use roc.area)
	
	for (iFic in 1:NBFIC) {

		#	1/ Data reading
		
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_wildtype.tsv", sep="")
		XMesNPLu 	<- fread(val, data.table=F)		# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_wildtype.tsv"
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockdowns.tsv", sep="")
		XMesPLu1 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/insilico_size10_1_knockdowns.tsv"
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/insilico_size", NBN, "_", iFic, "_knockouts.tsv", sep="")
		XMesPLu2 	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1//insilico_size10_1_knockouts.tsv"		
		val			<- paste(vRootDC4, "insilico_size", NBN, "_", iFic, "/DREAM4_GoldStandard_InSilico_Size", NBN, "_", iFic, ".tsv", sep="")
		SolutLu  	<- fread(val, data.table=F) 	# ex : "vRootDC4/insilico_size10_1/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
		SolutLu$V1	<- as.factor(SolutLu$V1)
		SolutLu$V2	<- as.factor(SolutLu$V2)
	
		for (i in 1:NBN) {
			XMesNP[i,1] = XMesNPLu[1,i]
			XMesNP[i,2] = XMesNPLu[1,i]
			
			for (j in 1:NBN) {
				XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation level -50%	
													# ATTENTION Dream Challenge data are transposed vs. my data 
				XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation level -100%
			}
		}	
	
		Solution2[,] = 0
		for (i in 1:(NBN*(NBN-1))) {
			if(SolutLu[i, ]$V3 == 1) {
				Solution2[as.numeric(SolutLu[i, ]$V1), as.numeric(SolutLu[i, ]$V2)] = 1
			} else {
				break								# We suppose that all the "1" are situated at the beginning of the file
			}
		}
	
		numero <- as.numeric(sort(as.character(c(1:NBN))))
		for (i in 1:NBN) {
			Solution3[ ,numero[i]] = Solution2[ ,i]
		}
		for (i in 1:NBN) {
			Solution[ ,numero[i]] = Solution3[i,]	# We find the classical connection matrix (nodes "G1", "G2", ... , "G10")
													# Note that Solution is transposed
		}

		for (iPaq in c(1:NBK)) {					# We take every available data : Knock Down and Knock Out
			#	The "input data", delivered by Dream Challenge, have been loaded before in XMesNP and XMesP
			MatR[ , ,iPaq]  <- 2 * (XMesP[ , ,iPaq]-XMesNP[ ,iPaq]) / (XMesP[ , ,iPaq]+XMesNP[ ,iPaq])
			
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
				
				for (iRow in 1:(NBN-1)) {
					pY[iRow] = MatR[iMat, col[iRow], iPaq]
					
					for (iCol in 1:NBN-1) {
						pX[iRow, iCol] = MatR[col[iCol], col[iRow], iPaq]
					}
				}
			
				gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,iMat]  	<- pX[ , ]
				gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),iMat]   	<- pY[ ]
			}
		}
	
		gXD 	<- gX[1 : (NBN-1), , ] 					# KnockDown
		gYD		<- gY[1 : (NBN-1), ]	
		gXO		<- gX[NBN : ((NBN-1)*2), , ] 			# KnockOut
		gYO		<- gY[NBN : ((NBN-1)*2), ] 
		
		for	(iMeth in 1:NBM) {
			Method  <- Methods[iMeth]
		
			cat("iSize ", iSize, " iFic ", iFic, " Method ", Method, " Heure ", as.character(Sys.time()), "\n")
			
			if (Method %in% c("MRA-KD", "MRA-KO")) {	#	Classical MRA method
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
					
					if (Method == "MRA-KD") {
						Form <- paste("gYD[ ,",iMat,"]~gXD[ , ,",iMat,"]+0")
					} else {
						Form <- paste("gYO[ ,",iMat,"]~gXO[ , ,",iMat,"]+0")
					}	
					
					MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}
				MatrCc[is.na(MatrCc)]  <- 0
			}	# MRA methods
	
			if (Method == "TLR") {
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
			
					Form <- paste("gY[ ,",iMat,"]~gX[ , ,",iMat,"]+0")
					MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}
				MatrCc[is.na(MatrCc)]  <- 0	
			}	# TLR with Linear regression, 2 measures per node
			
			if (Method == "LASSO") {				# Digitalization = 1 if the result != 0	-- automatic choice of lambda
				MatrCc[ , ] <- 0
				for (iMat in 1:NBN) {
					col 	<- 1:NBN
					col 	<- col[-iMat]
	
					cv_model 	<- cv.glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1)	# Fit lasso regression model using k-fold cross-validation
					best_lambda <- cv_model$lambda.min
					best_model 	<- glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1, lambda = best_lambda)		# View coefficients of best model
					rL			<-	coef(best_model)
					MatrCc[iMat, col]	<- rL[2:NBN]		# Note that MatrCc is transposed
				}
			}	# Method LASSO			

			if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
			#	Method "step" (automatic removal of coefficients by optimization of AIC : Akaike's Information Criterion)
			
				MatrCc[ , ] <- 0
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]	
					
					Donnees  <- data.frame(gY[ ,iMat], gX[ , ,iMat])
					#	cat("iMat ", iMat, " Data ", colnames(Donnees), "\n")
					Donnees	 <- rename(Donnees, "Y" ="gY...iMat.")		# The follow-up is clearer like this
					cc 		 <- colnames(Donnees)						# Name of the columns "Data". The first one is "Y"
					cc		 <- cc[-1]									# It remains the name of the coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
					colnames(rij) <- cc									# step gives the column names corresponding to the coefficients that are kept
							
					switch(Method, 
						"STEP-Fo" =
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},	# Forward
						"STEP-Ba" = 
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward
						"STEP-Bo" =
						{intercept_only <- lm(Y ~ 1, data=Donnees)
							all 	<- lm(Y ~ ., data=Donnees)
							forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)},		# Both
						)
					
					ll  	<- length(forward$coefficients)				# Number of coefficients kept by step
					nn		<- names(forward$coefficients)				# Name of these coefficients
			
					rij[1, ]	<- 0
					if (ll >= 2) {
						for (i in 2:ll) {								# nn[1] = "(Intercept)"
							rij[1, nn[i]] <- forward$coefficients[nn[i]]
						}
					}

					#	cat (" iMat ", iMat, " ll ", ll, " nn ", nn, " rij ", rij, "\n")

					MatrCc[iMat, col]	<- rij[1, ]
					MatrCc[iMat, iMat]  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
				}		# loop iMat
			}	# STEP methods
			
			if (Method == "CLR") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))	
				MatrCc  <- clr(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# CLR
					
			if (Method == "ARACNE") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))
				MatrCc  <- aracne(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# ARACNE

			if (Method == "MRNET") {		
				MI <- build.mim(dataset = rbind(MatR[,,1], MatR[,,2]))
				MatrCc  <- mrnet(MI)
				diag(MatrCc)  <- 0			# The diagonal element is set to 0 instead of -1, because the matrix "Solution" does so.
			}	# MRNET			

			#	Digitalization of the result and scoring	
			
			for (iSeuil in (1:NBSeuils)) {
				Seuil	<- Seuils[iSeuil]
				
				if (Method %in% c("MRA-KD", "MRA-KO")) {
					MatrCcDig	<-	Digitalize (MatrCc, "Top", Seuil)
				} else if (Method == "TLR") {
					MatrCcDig	<-	Digitalize (MatrCc, "Threshold", Seuil)					
				} else if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
					MatrCcDig	<- Digitalize (MatrCc, "Threshold", 0)
				}			
				
				Score[ ,iSeuil,iMeth,iFic,iSize]	<- Benchmark(MatrCcDig,Solution)
			}	# Loop on iSeuil	
				
			if (Method %in% c("MRA-KD", "MRA-KO")) {
				AuROC[iMeth,iFic,iSize]	<-  trapz(1-Score["Sp", ,iMeth,iFic,iSize], Score["Se", ,iMeth,iFic,iSize])# AUC ROC for the plotted methods
			} else if (Method %in% c("TLR", "LASSO")) {
				AuROC[iMeth,iFic,iSize]	<- -trapz(1-Score["Sp", ,iMeth,iFic,iSize], Score["Se", ,iMeth,iFic,iSize])# AUC ROC for the plotted methods
			}	# This method is less accurate than roc.area. When possible, use the following one (roc.area) to compute AUROC
			
			#	Computation of Area Under Curve (AUROC) and pValue (PVAUROC)
		
			vMatrCc[]  		<- as.numeric(abs(MatrCc[]))
			vSolution[]		<- as.numeric(Solution[])
			
			rb 	<- roc.area(vSolution, vMatrCc)
			ScoreDC["AUC", 	 iMeth, iFic, iSize]   = rb$A
			ScoreDC["PVAUC", iMeth, iFic, iSize] 	= rb$p.value					# Mann-Whitney U Test (wilcox.test)
		}		# Loop on iMeth

		Pts[ ,"TLR",iFic,iSize]		<-	Score[c("Se","Sp"),6,"TLR",iFic,iSize]		# 6 corresponds to the threshold 0.25	
	}			# Loop on iFic
}				# Loop on iSize

#	To get "Sensitivity" and "Specificity" corresponding to the values described in the article, we must choose :
#	iSeuil = 5 (Seuil = 0.2) if iSize = 1 and Method in ("MRA_KO", "MRA_KD")  (and "STEP" whose values don't depend on iSeuil)
#	iSeuil = 6 (Seuil = 0.25) if Method = "TLR" or if iSize = 2 and Method in ("MRA_KO", "MRA_KD")
#

cat("END ", as.character(Sys.time()), "\n")	

#
#	Supp. Info. Fig 2 a to e and Fig 3 a to e
#

Noms	<- c("Fig_SI2_", "Fig_SI3_")

iSize	<- 1			# 1 : insilico_Size10 (Fig_SI2),  2 : insilico_Size100 (Fig_SI3)
iFic	<- 1			# Scan the integers 1 to 5 to draw the different pictures ("a" to "e")

NomFic 	<- paste(vRoot, Noms[iSize], "_", letters[iFic], ".pdf", sep="")		# Fig_SI2_a.pdf
pdf(file=NomFic, height=7, width=7) 					# Sizes default to 7

plot (1-Score["Sp", ,"MRA-KD",iFic,iSize], Score["Se", ,"MRA-KD",iFic,iSize], type="l", col="blue", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
lines(1-Score["Sp", ,"MRA-KO",iFic,iSize],	Score["Se", ,"MRA-KO",iFic,iSize],	col="black")
lines(1-Score["Sp", ,"TLR",iFic,iSize],   	Score["Se", ,"TLR",iFic,iSize],    col="red")
abline(c(0,0), c(1,1), col="grey", lty="dashed")
points(1-Pts["Sp","TLR",iFic,iSize],   	Pts["Se","TLR",iFic,iSize], pch=23, col="orange", bg="orange")
segments(x0 = -0.03, y0 = 	Pts["Se","TLR",iFic,iSize], 	x1 = 1-Pts["Sp","TLR",iFic,iSize], y1 = Pts["Se","TLR",iFic,iSize], col = "darkgrey", lty="dashed")
text(0.01, Pts["Se","TLR",iFic,iSize], Pts["Se","TLR",iFic,iSize], font=3)				# font : 1=normal, 2=bold, 3=italic, 4=bold+italic   -- y coordinate of TLR, standard threshold
segments(x0 = 1-Pts["Sp","TLR",iFic,iSize], y0 = -0.03,	x1 = 1-Pts["Sp","TLR",iFic,iSize], y1 = Pts["Se","TLR",iFic,iSize], col = "darkgrey", lty="dashed")
text(1-Pts["Sp","TLR",iFic,iSize]+0.01, 0.03, 1-Pts["Sp","TLR",iFic,iSize], font=3)		# font : 1=normal, 2=bold, 3=italic, 4=bold+italic  -- x coordinate of TLR, standard threshold
text(0.03, 0.95, paste(Sizes[iSize], "_", iFic, sep=""), font=2)						# Name of the figure
#text(0.4,  0.03, "1 - Sp")																# x label
legend("bottomright", c("MRA-KD","MRA-KO","TLR"), fill=c("blue","black","red"))

dev.off()


#
#	Supp. Info. Fig 2 f and Fig 3 f + doc. "Answers to the reviewers (2 nd table § 1.1.1)
#

dimnames(ScoreDCOut)<- list(c("AUC",Methods), NULL, Sizes)
iSize	<- 2			# 1 : insilico_Size10 (Fig_SI2),  2 : insilico_Size100 (Fig_SI3)

for (iFic in 1:NBFIC) {
	ScoreDCOut[1,		 iFic,iSize]	<- paste(Sizes[iSize], "-", iFic, sep="")
	ScoreDCOut[2:(NBM+1),iFic,iSize]	<- round(ScoreDC["AUC", ,iFic,iSize], digits=2)
}

ScoreDCOut[1, NBFIC+1, iSize]			<- "Mean"
for (iMeth in 1:NBM) {
	ScoreDCOut[iMeth+1,NBFIC+1,iSize]	<- round(mean(ScoreDC["AUC",iMeth, ,iSize]), digits=2)
}	

ScoreDCOut[1,NBFIC+2,iSize]				<- "sd"
for (iMeth in 1:NBM) {
	ScoreDCOut[iMeth+1,NBFIC+2,iSize]	<- round(var(ScoreDC["AUC",iMeth, ,iSize])**0.5, digits=2)
}	

ScoreDCOut[1,NBFIC+3,iSize]				<- "Score"
for (iMeth in 1:NBM) {
	ScoreDCOut[iMeth+1,NBFIC+3,iSize]	<- round(-0.2*sum(log10(ScoreDC["AUC",iMeth, ,iSize])), digits=3)
}

write.csv(ScoreDCOut[ , ,iSize], file=paste(vRoot, Noms[iSize], "_f.csv", sep=""), col.names=FALSE, append=TRUE)


#
#	Doc. "Answers to the reviewers (ranking of our methods, §1.1.6: DREAM CHALLENGE 4)
#

dimnames(ScoreDCOut)<- list(c("PVAUC",Methods), NULL, Sizes)
iSize	<- 2			# 1 : insilico_Size10 (Fig_SI2),  2 : insilico_Size100 (Fig_SI3)

for (iFic in 1:NBFIC) {
	ScoreDCOut[1,		 iFic,iSize]	<- paste(Sizes[iSize], "-", iFic, sep="")
	ScoreDCOut[2:(NBM+1),iFic,iSize]	<- ScoreDC["PVAUC", ,iFic,iSize]
}

ScoreDCOut[1, NBFIC+1, iSize]			<- "Mean"
for (iMeth in 1:NBM) {
	ScoreDCOut[iMeth+1,NBFIC+1,iSize]	<- round(mean(ScoreDC["PVAUC",iMeth, ,iSize]), digits=4)
}

ScoreDCOut[1,NBFIC+2,iSize]				<- "sd"
for (iMeth in 1:NBM) {
	ScoreDCOut[iMeth+1,NBFIC+2,iSize]	<- round(var(ScoreDC["PVAUC",iMeth, ,iSize])**0.5, digits=4)
}

ScoreDCOut[1,NBFIC+3,iSize]				<- "Score"
for (iMeth in 1:NBM) {
	ScoreDCOut[iMeth+1,NBFIC+3,iSize]	<- round(-0.2*sum(log10(ScoreDC["PVAUC",iMeth, ,iSize])), digits=3)
}

write.csv(ScoreDCOut[ , ,iSize], file=paste(vRoot, Noms[iSize], "_f.csv", sep=""), col.names=FALSE, append=TRUE)


#
#		Reviewer 2 -- Major (1, 2 and 3)
#

#
#	Use of FRANK, a gene network simulator
#	https://m2sb.org/?page=FRANK
#

#
#	Drawing the curve "Quadratic Error" as a function of the size of the network.
#
#	To do this, we use "FRANK Network Generator" : "https://m2sb.org/?page=FRANK"
#	Ref : "Reverse engineering highlights potential principles of large gene regulatory network design and learning",
#	by Clément Caré, André Mas and Gabriel Krouk -- npj Systems Biology and Applications (2017)
#
#	This generator delivers a file (.csv) named ""Frankxxxxnonmodified_network.csv" (where xxxx is a number identifying the file, delivered by FRANK),
#	according to different parameters which are described below.
#	This generator simulates genes behavior : genes are regulated. TF act on other TFs and TAs, but TAs don't act on TFs. 
#	In our program, NBN (number of nodes) = TF+TA (TF represents the "number of Transcription Factors" and TA the "number of targets"). 
#	If TA=0, the file matches with a square matrix [NBN, NBN]. If TA>0, we get a rectangular matrix (NBN rows, TF columns).
#	In this case, we add TA null columns on the right, to get a square matrix [NBN, NBN], which is mandatory to use MRA.
#	Then, we replace the diagonal elements with -1, to comply with MRA requirements, and save this file with the name 
#	"Frank_TFxxx_TAyyy_z_Sol.csv" (xxx, yyy, z see below). This file represents the "Solution" used to score the results.
#	Starting from the "Solution", we compute the "Global Response Matrix R" (called MatR in the program. Read SI to know how to compute MatR).
#	MatR is an array [NBN, NBN, 2]. The first part corresponds to KnockOut perturbations (KO), the second part to KnockDown perturbations (-50%, called KD).
#	We compute "RMoy", the mean value of |MatR|.
#	Starting from MatR, we compute another array MatRN [NBN, NBN, 2, NBT] by adding to each element of MatR a "noise", which is an independant
#	random gaussian value N(0, sd), where sd = k*Rmoy (sd : standard deviation). NBT equals the number of noise levels. Here, NBT=2 and k=0.1 or 0.5.
#	This array is saved in a file named "Frank_TFxxx_TAyyy_z_R.csv".
#	For each size, we compute five independant files (z = 1,2 ..,5) allowing us to plot the boxplots showing the quadratic error as a
#	function of the network size, the method used and the noise level.
#	For each network size, method used and noise level, we compute the "result" (a matrix named MatrCc[NBN, NBN]) and the quadratic error
#	is defined by ErrorQ = sqrt((MatrCc-Solution)**2)/NBN. Dividing by NBN "normalizes" the error.
#
#	FRANK Network Generator parameters :
#		Number of Transcription Factors		 : TF number (xxx in the file name)
#	 	Number of TArget genes				 : TA number (yyy in the file name)
#	 	Number of eigenvalues of the TF matrix on the unit circle	1   (to get a system leading to a steady state, so as to comply with MRA requirements)
#	 	Value of the variance of the noise added to the log-normal observations		0 (default value, meaning no noise added by FRANK)
#	 	Time-serie observations 			 : not concerned. Choose dynamic 
#		nb expce 							 : not concerned. Choose 1
#	 
#	 	Seed of the random value generator	 :  define a different value for each of the 5 files (see below).
#	 	Minimum sparsity (number of nonzero elements per row) of the TF and TA matrices		: set to 0.15*(TF+TA),   to have a sparsity of 15%
#	 	Maximum sparsity (number of nonzero elements per row) of the TF and TA matrices		: set to 0.15*(TF+TA)+1, to have a sparsity of 15%
#	 	Slop for the sparsity of the matrices												: -2  (default value)
#	 	Magnitude of deviation from zero of the nonnull elements of the matrix				: 1   (default value)
#	 	Value of the mean of the log-normal probability distribution of the observations	: 5   (default value)	
#	 	Value of the variance of the log-normal probabilty distribution of the observations	: 3.5 (default value)
#
#	To get five independant files (z = 1,2 ..,5), we use 5 different "Seed" values, chosen at random : 
#	Seeds = c(12345, 72090, 87577, 45648, 16637). These seeds are also used as "seed" for the noise generator.
#
#	We will study these combinations of networks (the number of nodes, NBN, equals TF+TA) :
#	(TF=30, TA=0), (TF=60, TA=0), (TF=100, TA=0), (TF=200, TA=0), (TF=300, TA=0), (TF=500, TA=0), (TF=800, TA=0), (TF=1000, TA=0),
#	(TF=30, TA=30), (TF=50, TA=50), (TF=100, TA=100), (TF=150, TA=150), (TF=250, TA=250), (TF=400, TA=400), (TF=500, TA=500),   
#	and the following methods :
#	- MRA (it means "classical" MRA, described by Kholodenko et al). As MRA cannot use KO and KD data simultaneously, we will compute
#	MatrCc for KO and KD data successively and then compute the average of these two values, to allow a fair comparison of the methods (same number of measurements
#	for all methods).
#	- LSE_CI (MRA computed by means of standard regression method : minimize Least Square Error : LSE). This method uses KO and KD data simultaneously.
#	- 2 methods integrating a variable selection scheme, namely STEP-Fo (Forward) and STEP-Ba (Backward). These methods also use KO and KD data 
#	simultaneously.
#
#	The original data delivered by FRANK Network Generator are archived, in case more computations are needed. Files are named ""Frankxxxxnonmodified_network.csv",
#	where xxxx are :
#	  11,  649, 8363, 5883, 6798	for TF=30,   TA=0, (z = 1,2 ..5 and Seed = Seeds[z]), 
#	4685, 3537, 1797, 7380, 2568	for TF=60,   TA=0,
#	9725, 9163, 2648, 9991, 8602	for TF=100,  TA=0,
#	2776, 5669, 5944, 6979, 6139	for TF=200,  TA=0,
#	7374, 1110, 5563, 8121, 1863	for TF=300,  TA=0,
#	6283, 9679, 2542,  730, 9564	for TF=500,  TA=0,
#	9518, 3919, 8666, 3067, 8099	for TF=800,  TA=0,
#	1619, 9628, 4013, 3416, 3148	for TF=1000, TA=0,
#
#	5191, 5303, 3890, 9658, 5309	for TF=30,   TA=30,
#	7120, 8664,  239, 3835, 4999	for TF=50,	 TA=50,
#	1276, 3865, 7672, 2237, 7584	for TF=100,  TA=100,
#	 920, 1311,  499, 3243, 1876	for TF=150,	 TA=150,
#	4806, 7511, 3494, 3375, 5994	for TF=250,  TA=250,
#	9332  3725, 7218, 5265, 7952	for TF=400,  TA=400,
#	2072, 6819, 6975, 8518, 5660	for TF=500,  TA=500.
#
#	The results (array ErrorQ, for "quadratic error") are saved in the file "Frank_ErrorQ.csv"

#
#	Joint statements
#

vRoot		<- c("/home/jpborg/Documents/Publications/BioInformatics/MRA_Regression_ed1/Vers_Rev1/Donnees/FRANK2/")

NBSets		<-  15				# Nbr of sets of values (TF/TA) tested
NBTrials	<-  5				# Nbr of trials by set of values
NBK			<-	2				# Nbr of perturbations (KO, KD -50%)
Noises		<- c(0.1, 0.5)		# Error coefficients
NBT			<- length(Noises)

Diago		<- c(-2, -2/3)		# KO, KD (-50%)
Seeds 		<- c(12345, 72090, 87577, 45648, 16637)		# To do independant trials

Methods		<- c("MRA", "LSE_CI", "LASSO", "LSE", "STEP-Fo", "CLR", "MRA_CLR", "ARACNE", "MRNET")
NBM			<- length(Methods)								# Nbr of methods to test

ErrorQ		<- array(0, dim=c(NBM,NBT,NBTrials,NBSets))		# Quadratic error, as a function of noise level and method used
dimnames(ErrorQ)	<- list(Methods,Noises,NULL,NULL)

FileErrorQ	<- paste(vRoot, "Resultats/", "Frank_ErrorQ.csv", sep="")
FileScore	<- paste(vRoot, "Resultats/", "Frank_Scores.csv", sep="")
val			<- paste(vRoot, "Donnees/",   "FileNbr.csv", sep="")
FileNbr		<- as.matrix(fread(val, data.table=F), ncol=NBTrials+2)
					# Column name	: TF 	TA    1 	 2 	   3 	 4 	   5
					# First row		: Seed	 	12345 72090 87577 45648 16637
					# Following rows: TF val TA val File nbrs		

#
#	First part : creation of the files "Frank_TFxxx_TAyyy_z_Sol.csv" 
#

for (iSet in 1:NBSets) {
	TF			<- as.numeric(FileNbr[iSet, 1])
	TA			<- as.numeric(FileNbr[iSet, 2])
	NBN			<- TF+TA
	
	for (iTrial in 1:NBTrials) {
		val			<- paste(vRoot, "Reseau/", "Frank", FileNbr[iSet, iTrial+2], "nonmodified_network.csv", sep="")
		cat("iSet ", iSet, " TF ", TF, " TA ", TA, " iTrial ", iTrial, " FicLu ", val, "\n")
		Solution 	<- as.matrix(read.table(val))
									# As there is a header and the first row contains one fewer field than the number of columns, 
									# the first column in the input is used for the row names. 

		if (dim(Solution)[1] > NBN) {
			rowDel		<- c((NBN+1):(dim(Solution)[1]))	# row(s) to delete (added by FRANK by mistake)
			Solution	<- Solution[-rowDel,]
		}
		if (TA > 0) {
			Solut1 		<- array(0, dim=c(NBN, TA))			# TA genes don't act on TF genes
			Solution	<- cbind(Solution, Solut1)
		}
		diag(Solution)	<- -1
		
		val			<- paste(vRoot, "Donnees/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Sol.csv", sep="")
		write.table(Solution, val, sep=",")		# "Frank_TFxxx_TAyyy_z_Sol.csv"
		
#
#	Second part : creation of the files "Frank_TFxxx_TAyyy_z_R.csv" (noisy Global Response Matrix)
#
		
		rm1 		<- as.matrix(solve(Solution), nrow=NBN,ncol=NBN)			# r ** -1
		dgrm1M1		<- array(0, dim=c(NBN, NBN))								# [dg(r**-1)]**-1
		MatR		<- array(dim=c(NBN, NBN, NBK))				# Matrix R (KO, KD)
		MatRN		<- array(dim=c(NBN, NBN, NBK, NBT))			# Matrix R with noise
		
		for (i in 1:NBN) {
			dgrm1M1[i,i] = 1/rm1[i,i]
		}
		
		MatR[ , ,1]		<- (rm1 %*% dgrm1M1)*Diago[1]	# 1 = KO
		MatR[ , ,2]		<- (rm1 %*% dgrm1M1)*Diago[2]	# 2 = KD (-50%)
		RMoy			<- mean(abs(MatR))				# mean value
		
		cat("TF ", TF, " TA ", TA, " Fic ", iTrial, " RMoy ", RMoy, "\n")		
		
		for (iT in 1:NBT) {
			set.seed(Seeds[iTrial])						# So as to generate the same sequences
			Noise	<- Noises[iT]
			
			for (iRow in 1:NBN) {
				for (iCol in 1:NBN) {
					for (iK in 1:NBK) {
						MatRN[iRow, iCol, iK, iT] = MatR[iRow, iCol, iK] + rnorm(1, mean=0, sd=Noise*RMoy)
					}	# Loop on iK
				}		# Loop on iCol
			}			# Loop on iRow
		}				# Loop on iT
		
		val			<- paste(vRoot, "Donnees/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_R.csv", sep="")
		write.table(MatRN, val, sep=",")				# "Frank_TFxxx_TAyyy_z_R.csv" (noisy Global Response Matrix)
	}				# Loop on iTrial
}					# Loop on iSet

#
#	Third part : detection of networks structure and error computation
#	Supplementary Information Table 4
#

for (iSet in 1:NBSets) {
	TF			<- as.numeric(FileNbr[iSet, 1])
	TA			<- as.numeric(FileNbr[iSet, 2])
	NBN			<- TF+TA
	
	MatRN		<- array(dim=c(NBN, NBN, NBK, NBT))			# Matrix R with noise	
	MatRNm1		<- array(dim=c(NBN, NBN, NBK))				# RN ** -1
	dgRNm1M1	<- array(0, dim=c(NBN, NBN, NBK))			# [dg(RN**-1)]**-1
	
	MatrCc		<- array(dim=c(NBN, NBN))					# "r" matrix computed, according to the method used
	rij 		<- matrix(nrow=1, ncol=NBN-1)				# Intermediate calculation of rij
	pX   		<- array(dim=c(NBN-1, NBN-1))				# row (iRow), column (iCol)	
	pY   		<- array(dim=c(NBN-1))						# row
	gX   		<- array(dim=c(NBK*(NBN-1), NBN-1, NBN))	# row (iRow), column (iCol), matrix (iMat)
	gY   		<- array(dim=c(NBK*(NBN-1), NBN))			# row, matrix
	
	rL 			<- array(dim=c(NBN))						# Result of Lasso method (ie sol. of Yi = Ai * Xi)
	
	ValMin		<- array(dim=c(NBN, NBN))					# Min value at IC95 of rij
	ValMax		<- array(dim=c(NBN, NBN))					# Max value at IC95 of rij


	for (iTrial in 1:NBTrials) {
		cat("iSet ", iSet, " TF ", TF, " TA ", TA, " iTrial ", iTrial, "\n")

		val			<- paste(vRoot, "Donnees/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Sol.csv", sep="")
		Solution	<- as.matrix(read.csv(val))				# "Frank_TFxxx_TAyyy_z_Sol.csv"

		val		<- paste(vRoot, "Donnees/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_R.csv", sep="")
		MatRN	<- array(unlist(read.csv(val)),dim=c(NBN, NBN, NBK, NBT))	# "Frank_TFxxx_TAyyy_z_R.csv" (noisy Global Response Matrix)

		for (iT in 1:NBT) {
			for (iK in 1:NBK) {
				MatRNm1[ , ,iK]		 <- solve(MatRN[ , ,iK,iT])
				for (i in 1:NBN) {
					dgRNm1M1[i,i,iK] <- 1/MatRNm1[i,i,iK]
				}
			}
			
			for (iPaq in c(1:NBK)) {					# We take every available data : Knock Down (50%) and Knock Out
				for (iMat in 1:NBN) {
					col <- 1:NBN
					col <- col[-iMat]
					
					for (iRow in 1:(NBN-1)) {
						pY[iRow] = MatRN[iMat, col[iRow], iPaq, iT]
						
						for (iCol in 1:NBN-1) {
							pX[iRow, iCol] = MatRN[col[iCol], col[iRow], iPaq, iT]
						}
					}
				
					gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,iMat]  	<- pX[ , ]
					gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),iMat]   	<- pY[ ]
				}		# Loop on iMat
			}			# Loop on iPaq			
						
			for (iMeth in 1:NBM) {
				Method 	<- Methods[iMeth]
				cat(" iT ", iT, " Method ", Method, " Time ", as.character(Sys.time()), "\n")
							
				if (Method == "MRA") {
					MatrCc	<- -0.5 * (dgRNm1M1[ , ,1] %*% MatRNm1[ , ,1] + dgRNm1M1[ , ,2] %*% MatRNm1[ , ,2])
				}	# Classical MRA

				if (Method == "MRA_CLR") {					# "minet" CLR is applied on the highest values of MatrCc
					rijMRA	<- -0.5 * (dgRNm1M1[ , ,1] %*% MatRNm1[ , ,1] + dgRNm1M1[ , ,2] %*% MatRNm1[ , ,2])					
					Th		<- 0.25*max(abs(rijMRA))
					rijTh	<- ifelse(abs(rijMRA) >= Th, rijMRA, 0)					
					MI 		<- build.mim(dataset = rijTh)			
					MI		<- ifelse(is.na(MI), 0, MI)
					MatrCc  <- clr(MI)
					MatrCc  <- ifelse(is.nan(MatrCc), 0, MatrCc)
				}	# MRA_CLR	

				if (Method == "LSE") {			# CAUTION : this method is called "TLR" in the article. No threshold is applied when computing quadratic error.
					for (iMat in 1:NBN) {
						col 	<- 1:NBN
						col 	<- col[-iMat]
						Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")
						
						MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
					}		# Loop on iMat						
				}	# Linear regression, 2 measures per node
								
				if (Method == "LSE_CI") {
					for (iMat in 1:NBN) {
						col 	<- 1:NBN
						col 	<- col[-iMat]
						Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")
						
						pQ 	 <- lm(as.formula(Form))
						ValMin [iMat,col]  	<- broom::tidy(pQ)$estimate - 1.96*broom::tidy(pQ)$std.error	# 1.96 to get an IC 95%
						ValMax [iMat,col]  	<- broom::tidy(pQ)$estimate + 1.96*broom::tidy(pQ)$std.error
						
						MatrCc[iMat, col] 	<- ifelse(ValMin[iMat,col] > 0 | ValMax[iMat,col] < 0, 1, 0)
					}		# Loop on iMat						
				}	# LSE_CI method
								
				if (Method == "TLR") {
					for (iMat in 1:NBN) {
						Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")	
						MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
					}		# Loop on iMat
					
					Th		<- 0.25*max(abs(MatrCc))
					MatrCc	<- ifelse(abs(MatrCc) >= Th, MatrCc, 0)						
				}	# Linear regression with a threshold, 2 measures per node

				if (Method == "LASSO") {					# Automatic choice of lambda
					for (iMat in 1:NBN) {					
						col 	<- 1:NBN
						col 	<- col[-iMat]
		
						cv_model 	<- cv.glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1)	# Fit lasso regression model using k-fold cross-validation
						best_lambda <- cv_model$lambda.min
						best_model 	<- glmnet(gX[ , ,iMat], gY[ ,iMat], alpha = 1, lambda = best_lambda)		# View coefficients of best model
						rL			<-	coef(best_model)
						MatrCc[iMat, col]	<- rL[2:NBN]		# Note that MatrCc is transposed						
					}
				}	# LASSO method

				if (Method %in% c("STEP-Fo", "STEP-Ba")) {
				#	Method "step" (automatic removal of coefficients by optimization of AIC : Akaike's Information Criterion)
				
					MatrCc[,] <- 0
					for (iMat in 1:NBN) {		
						for (iPaq in c(1:NBK)) {			# We take every available data : Knock Out and Knock Down (50%)	
							col <- 1:NBN
							col <- col[-iMat]
							
							for (iRow in 1:(NBN-1)) {
								pY[iRow] = MatRN[iMat, col[iRow], iPaq, iT]
								
								for (iCol in 1:NBN-1) {
									pX[iRow, iCol] = MatRN[col[iCol], col[iRow], iPaq, iT]
								}
							}
										
							gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,iMat]  	<- pX[ , ]
							gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),iMat]   	<- pY[ ]							
						}	# Loop on iPaq

						Donnees  <- data.frame(gY[ ,iMat], gX[ , ,iMat])
						#	cat("iMat ", iMat, " Data ", colnames(Donnees), "\n")							
						Donnees	 <- rename(Donnees, "Y" ="gY...iMat.")		# The follow-up is clearer like this
						cc 		 <- colnames(Donnees)						# Name of the columns "Data". The first one is "Y"
						cc		 <- cc[-1]									# It remains the name of the coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
						colnames(rij) <- cc									# step gives the column names corresponding to the coefficients that are kept
							
						switch(Method, 
							"STEP-Fo" =
							{intercept_only <- lm(Y ~ 1, data=Donnees)
								all 	<- lm(Y ~ ., data=Donnees)
								forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},	# Forward
							"STEP-Ba" = 
							{intercept_only <- lm(Y ~ 1, data=Donnees)
								all 	<- lm(Y ~ ., data=Donnees)
								forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward
							"STEP-Bo" =
							{intercept_only <- lm(Y ~ 1, data=Donnees)
								all 	<- lm(Y ~ ., data=Donnees)
								forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)},		# Both								
							)
						
						ll  	<- length(forward$coefficients)				# Number of coefficients kept by step
						nn		<- names(forward$coefficients)				# Name of these coefficients
				
						rij[1,]	<- 0
						if (ll >= 2) {
							for (i in 2:ll) {								# nn[1] = "(Intercept)"
								rij[1, nn[i]] <- forward$coefficients[nn[i]]
							}
						}
						
						MatrCc[iMat, col]	<- rij[1,]
					}	# Loop on iMat
				}	# STEP methods					

				if (Method == "ARACNE") {				# Library "minet"
					MI <- build.mim(dataset = rbind(MatRN[,,1,iT], MatRN[,,2,iT]))	
					MatrCc  <- aracne(MI)
				}	# ARACNE method
				
				if (Method == "CLR") {					# Library "minet"
					MI <- build.mim(dataset = rbind(MatRN[,,1,iT], MatRN[,,2,iT]))	
					MatrCc  <- clr(MI)
				}	# CLR method
				
				if (Method == "MRNET") {		
					MI <- build.mim(dataset = rbind(MatRN[,,1,iT], MatRN[,,2,iT]))	
					MatrCc  <- mrnet(MI)
				}	# MRNET	method				
				
				diag(MatrCc)	<- -1
				
				ss	<- sum((MatrCc-Solution)**2)
				ErrorQ[iMeth,iT,iTrial,iSet]	<- sqrt(ss)

				val		<- paste("\n\n\nFrank_TF", TF, "_TA", TA, "_", iTrial, "_Cc", 10*Noises[iT], "_", Method, "\n", sep="")
				write.table(val, FileErrorQ, append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
				write.table(ErrorQ, FileErrorQ, append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")	# Frank_ErrorQ.csv
				# 	This programm may be very time consuming. These partial results are saved in case of interruption.

				val		<- paste(vRoot, "Resultats/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Cc", 10*Noises[iT], "_", Method, ".csv", sep="")
				write.table(MatrCc, val, row.names=FALSE, col.names=FALSE, sep=",")			# Frank_TF30_TA0_1_Cc1_MRA.csv	
			}		# Loop on iMeth		
		}			# Loop on iT
		
#		val		<- paste("\n\n\nFrank_TF", TF, "_TA", TA, "_", iTrial, "\n", sep="")
#		write.table(val, FileErrorQ, append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
#		write.table(ErrorQ[ , ,iTrial,iSet], FileErrorQ, append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")	# Frank_ErrorQ.csv		
	}				# Loop on iTrial
}					# Loop on iSet

write.table(ErrorQ, FileErrorQ, append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")	# Frank_ErrorQ.csv
# If the program reaches the end, the errors are saved here and the temporary values are deleted.

cat("FIN Time ", as.character(Sys.time()), "\n")


#
#	Fourth part : curves drawing (Quadratic Error)
#

NbNodes		<- vector(length=NBSets)								# Number of nodes
Methods		<- c("MRA", "MRA_CLR", "ARACNE", "LSE",  "STEP-Fo")
NBM			<- length(Methods)

RightA		<- array(0, dim=c(NBTrials, NBSets))					# L1 norm : sum(abs(Solution[i,j]))
RightM		<- array(0, dim=c(2,NBSets))							# Statistics about L1 norm : mean, sd
dimnames(RightM)	<- list(c("mean", "sd"), NULL)

val			<- paste(vRoot, "Resultats/", "Frank_ErrorQ.csv", sep="")
ErrorsA		<- array(unlist(fread(val,data.table=F)),dim=c(NBM,NBT,NBTrials,NBSets))		# Errors stored previously (ErrorQ)
dimnames(ErrorsA)	<- list(Methods, Noises, NULL, NULL)
ErrorsM		<- array(0, dim=c(4, NBM, NBT, NBSets))				# Statistics about quadratic errors (mean, sd and log on the different trials)
dimnames(ErrorsM)	<- list(c("mean", "sd", "mean(Log)", "sd(Log)"), Methods, Noises, NULL)

for (iSet in 1:NBSets) {
	TF			<- as.numeric(FileNbr[iSet, 1])
	TA			<- as.numeric(FileNbr[iSet, 2])
	
	for (iTrial in 1:NBTrials) {
		val			<- paste(vRoot, "Donnees/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Sol.csv", sep="")
		Solution	<- as.matrix(read.csv(val))	# "Frank_TFxxx_TAyyy_z_Sol.csv"
		RightA[iTrial,iSet]		<- sum(abs(Solution))
	}		# Loop on iTrial
}			# Loop on iSet	

for (iSet in 1:NBSets) {
	RightM["mean",iSet] 	<- round (mean(RightA[ ,iSet]), digits=3)			# Table RIGHT
	RightM["sd",iSet] 		<- round ((var(RightA[ ,iSet]))**0.5, digits=3)
}			# Loop on iSet

#
#	Mean and standard deviation (computed from the NBTrials results) of quadratic error, as a function of the number of nodes (iSet), the noise level (iT) and the method used (iMeth)
#

for (iSet in 1:NBSets) {
	TF			<- as.numeric(FileNbr[iSet, 1])
	TA			<- as.numeric(FileNbr[iSet, 2])
	NbNodes[iSet]	<- TF+TA
	
	for (iT in 1:NBT)  {
		for (iMeth in 1:NBM) {
			ErrorsM["mean",	iMeth,iT,iSet]  	= round (mean(as.numeric(ErrorsA[iMeth,iT, ,iSet])), 				digits=3)
			ErrorsM["sd", 	iMeth,iT,iSet]  	= round ((var(as.numeric(ErrorsA[iMeth,iT, ,iSet])))**0.5, 			digits=3)
			ErrorsM["mean(Log)",iMeth,iT,iSet] 	= round (mean(log1p(as.numeric(ErrorsA[iMeth,iT, ,iSet]))),  		digits=3)	# log1p : Ln(1+x)
			ErrorsM["sd(Log)",	iMeth,iT,iSet] 	= round ((var(log1p(as.numeric(ErrorsA[iMeth,iT, ,iSet]))))**0.5,	digits=3)			
		}		# Loop on iMeth
	}			# Loop on iT
}				# Loop on iSet

#
#	New Figure 5
#
#	TA = 0		k = 0.1

Err1TA0_MRA.df		<- data.frame(NbNodes=NbNodes[1:8], avg=ErrorsM["mean","MRA",1,1:8], 	 	sd=ErrorsM["sd","MRA",1,1:8], 		Method="MRA")
# Err1TA0_MRACLR.df	<- data.frame(NbNodes=NbNodes[1:8], avg=ErrorsM["mean","MRA_CLR",1,1:8], 	sd=ErrorsM["sd","MRA_CLR",1,1:8], 	Method="MRA-CLR")
# Err1TA0_ARAC.df		<- data.frame(NbNodes=NbNodes[1:8], avg=ErrorsM["mean","ARACNE",1,1:8], 	sd=ErrorsM["sd","ARACNE",1,1:8], 	Method="ARACNE")
Err1TA0_LSE.df		<- data.frame(NbNodes=NbNodes[1:8], avg=ErrorsM["mean","LSE",1,1:8], 	 	sd=ErrorsM["sd","LSE",1,1:8], 		Method="TLR")
Err1TA0_STEP.df		<- data.frame(NbNodes=NbNodes[1:5], avg=ErrorsM["mean","STEP-Fo",1,1:5],	sd=ErrorsM["sd","STEP-Fo",1,1:5], 	Method="STEP")
Err1TA0.df			<- rbind(Err1TA0_MRA.df, Err1TA0_LSE.df, Err1TA0_STEP.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Err1TA0.pdf", sep="")
pdf(file=NomFic, height=6,width=3) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Err1TA0.df, aes(x=NbNodes, y=avg, colour=factor(Method)), show.legend=FALSE)+	# legend removed to get more space (4 figs. with the same caption, added with Inkscape)
	xlab("Number of nodes") + ylab("Average squared error") + ggtitle("Squared Error, k = 10%, TA = 0")   +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "ARACNE"="green", "TLR"="red", "STEP"="black")) +
	geom_ribbon(data=Err1TA0_LSE.df,  aes(x=NbNodes, ymin=avg-sd, ymax=avg+sd), fill="grey70", 	alpha=0.3) +
	geom_ribbon(data=Err1TA0_STEP.df, aes(x=NbNodes, ymin=avg-sd, ymax=avg+sd), fill="olivedrab2", alpha=0.3) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()

#	TF = TA		k = 0.1

Err1TATF_MRA.df		<- data.frame(NbNodes=NbNodes[9:15], avg=ErrorsM["mean","MRA",1,9:15], 	 	sd=ErrorsM["sd","MRA",1,9:15], 		Method="MRA")
# Err1TA0_MRACLR.df	<- data.frame(NbNodes=NbNodes[9:15], avg=ErrorsM["mean","MRA_CLR",1,9:15], 	sd=ErrorsM["sd","MRA_CLR",1,9:15], 	Method="MRA-CLR")
# Err1TATF_ARAC.df	<- data.frame(NbNodes=NbNodes[9:15], avg=ErrorsM["mean","ARACNE",1,9:15], 	sd=ErrorsM["sd","ARACNE",1,9:15], 	Method="ARACNE")
Err1TATF_LSE.df		<- data.frame(NbNodes=NbNodes[9:15], avg=ErrorsM["mean","LSE",1,9:15], 	 	sd=ErrorsM["sd","LSE",1,9:15], 		Method="TLR")
Err1TATF_STEP.df	<- data.frame(NbNodes=NbNodes[9:11], avg=ErrorsM["mean","STEP-Fo",1,9:11], 	sd=ErrorsM["sd","STEP-Fo",1,9:11], 	Method="STEP")
Err1TATF.df			<- rbind(Err1TATF_MRA.df, Err1TATF_LSE.df, Err1TATF_STEP.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Err1TATF.pdf", sep="")
pdf(file=NomFic, height=6,width=3) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Err1TATF.df, aes(x=NbNodes, y=avg, colour=factor(Method)), show.legend=FALSE)+
	xlab("Number of nodes") + ylab("Average squared error") + ggtitle("Squared Error, k = 10%, TA = TF")   +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "ARACNE"="green", "TLR"="red","STEP"="black"))   +
	geom_ribbon(data=Err1TATF_LSE.df, aes(x=NbNodes,  ymin=avg-sd, ymax=avg+sd), fill="grey70", 	 alpha=0.3) +
	geom_ribbon(data=Err1TATF_STEP.df, aes(x=NbNodes, ymin=avg-sd, ymax=avg+sd), fill="olivedrab2", alpha=0.3)	+
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()

#	TA = 0		k = 0.5

Err5TA0_MRA.df		<- data.frame(NbNodes=NbNodes[1:8], avg=ErrorsM["mean(Log)","MRA",2,1:8], 	  	sd=ErrorsM["sd(Log)","MRA",2,1:8], 		Method="MRA")
# Err5TA0_MRACLR.df	<- data.frame(NbNodes=NbNodes[1:8], avg=ErrorsM["mean(Log)","MRA_CLR",2,1:8], 	sd=ErrorsM["sd(Log)","MRA_CLR",2,1:8], 	Method="MRA-CLR")
# Err5TA0_ARAC.df	<- data.frame(NbNodes=NbNodes[1:8], avg=ErrorsM["mean(Log)","ARACNE",2,1:8], 	sd=ErrorsM["sd(Log)","ARACNE",2,1:8], 	Method="ARACNE")
Err5TA0_LSE.df		<- data.frame(NbNodes=NbNodes[1:8], avg=ErrorsM["mean(Log)","LSE",2,1:8], 	  	sd=ErrorsM["sd(Log)","LSE",2,1:8], 		Method="TLR")
Err5TA0_STEP.df		<- data.frame(NbNodes=NbNodes[1:5], avg=ErrorsM["mean(Log)","STEP-Fo",2,1:5],	sd=ErrorsM["sd(Log)","STEP-Fo",2,1:5], 	Method="STEP")
Err5TA0.df			<- rbind(Err5TA0_MRA.df, Err5TA0_LSE.df, Err5TA0_STEP.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Err5TA0.pdf", sep="")
pdf(file=NomFic, height=6,width=3.5) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Err5TA0.df, aes(x=NbNodes, y=avg, colour=factor(Method)), show.legend=FALSE)+
	xlab("Number of nodes") + ylab("Log(1+Average squared error)") + ggtitle("Log(1+Sq. Error), k = 50%, TA = 0") +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "MRA-CLR"="orange", "ARACNE"="green", "TLR"="red", "STEP"="black")) +
	geom_ribbon(data=Err5TA0_LSE.df,  aes(x=NbNodes, ymin=avg-sd, ymax=avg+sd), fill="grey70", 	alpha=0.3) +
	geom_ribbon(data=Err5TA0_STEP.df, aes(x=NbNodes, ymin=avg-sd, ymax=avg+sd), fill="olivedrab2", alpha=0.3) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3.5 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()

#	TF = TA		k = 0.5

Err5TATF_MRA.df		<- data.frame(NbNodes=NbNodes[9:15], avg=ErrorsM["mean(Log)","MRA",2,9:15], 	 	sd=ErrorsM["sd(Log)","MRA",2,9:15], 		Method="MRA")
# Err5TATF_MRACLR.df<- data.frame(NbNodes=NbNodes[9:15], avg=ErrorsM["mean(Log)","MRA_CLR",2,9:15], 	sd=ErrorsM["sd(Log)","MRA_CLR",2,9:15], 	Method="MRA-CLR")
# Err5TATF_ARAC.df	<- data.frame(NbNodes=NbNodes[9:15], avg=ErrorsM["mean(Log)","ARACNE",2,9:15], 		sd=ErrorsM["sd(Log)","ARACNE",2,9:15], 		Method="ARACNE")
Err5TATF_LSE.df		<- data.frame(NbNodes=NbNodes[9:15], avg=ErrorsM["mean(Log)","LSE",2,9:15], 	  	sd=ErrorsM["sd(Log)","LSE",2,9:15], 		Method="TLR")
Err5TATF_STEP.df	<- data.frame(NbNodes=NbNodes[9:11], avg=ErrorsM["mean(Log)","STEP-Fo",2,9:11], 	sd=ErrorsM["sd(Log)","STEP-Fo",2,9:11], 	Method="STEP")
Err5TATF.df			<- rbind(Err5TATF_MRA.df, Err5TATF_LSE.df, Err5TATF_STEP.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Err5TATF.pdf", sep="")
pdf(file=NomFic, height=6,width=3.5) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Err5TATF.df, aes(x=NbNodes, y=avg, colour=factor(Method)), show.legend=FALSE) +
	xlab("Number of nodes") + ylab("Log(1+Average squared error)") + ggtitle("Log(1+Sq. Error), k = 50%, TA = TF") +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "MRA-CLR"="orange", "ARACNE"="green", "TLR"="red", "STEP"="black")) +
	geom_ribbon(data=Err5TATF_LSE.df,  aes(x=NbNodes, ymin=avg-sd, ymax=avg+sd), fill="grey70", alpha=0.3) +
	geom_ribbon(data=Err5TATF_STEP.df, aes(x=NbNodes, ymin=avg-sd, ymax=avg+sd), fill="olivedrab2", alpha=0.3) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3.5 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()


#
#	Fifth part : curves drawing (ROC curves)
#	Supplementary Information Table 5
#

NbNodes			<- vector(length=NBSets)				# Number of nodes

Seuils			<- seq(0, 1, by=0.05)					# To plot the curves ROC (MRA, MRA_CLR, TLR)
NBSeuils		<- length(Seuils)

Noises			<- c(0.1, 0.5)							# Error coefficients
NBT				<- length(Noises)
#Methods			<- c("MRA", "MRA_CLR", "CLR", "ARACNE", "LSE", "TLR", "STEP-Fo")
Methods			<- c("MRA", "LSE")
NBM				<- length(Methods)

ScoreTests		<-	c("TP", "TN", "FP", "FN", "Se", "Sp")		# Scores to measure (Se = "Sensibility", Sp = "Specificity", AUC = "AUROC")
NBSC			<- length(ScoreTests)					# Nbr. of scores tested
Score 			<- array(0, dim=c(NBSC, NBSeuils, NBM, NBT, NBTrials, NBSets))
dimnames(Score)	<- list(ScoreTests, NULL, Methods, Noises, NULL, NULL)
AUROC.tmp		<- array(0, dim=c(NBM, NBT, NBTrials, NBSets))	
AUROC 			<- array(0, dim=c(NBM, NBT, NBSets))
dimnames(AUROC)	<- list(Methods, Noises, NULL)


for (iSet in 1:NBSets) {
	TF			<- as.numeric(FileNbr[iSet, 1])
	TA			<- as.numeric(FileNbr[iSet, 2])
	NBN			<- TF+TA
	NbNodes[iSet]	<- NBN
	
	Solution	<-	array(dim=c(NBN, NBN))				# "Real" value of the solution
	Solut1		<-	array(dim=c(NBN, NBN))				# Digitalized solution (0, 1)
	MatrCc		<-  array(dim=c(NBN, NBN))				# "r" matrix computed, according to the method used
	MatrCc1		<-  array(dim=c(NBN, NBN))				# Digitalized "r" matrix (0, 1)
	
		
	for (iTrial in 1:NBTrials) {
		val			<- paste(vRoot, "Donnees/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Sol.csv", sep="")
		Solution	<- as.matrix(read.csv(val))			# "Frank_TFxxx_TAyyy_z_Sol.csv"
		Solut1		<- Digitalize(Solution, "Top", 0.15)	# We keep the 15% highest absolute values

		for (iT in 1:NBT) {
			Noise	<- 10*Noises[iT]
			
			for (iMeth in 1:NBM) {
				Method	<- Methods[iMeth]
				cat("TF ", TF, " TA ", TA, " Trial ", iTrial, " iT ", iT, " Method ", Method, " Time ", as.character(Sys.time()), "\n")
				
				val		<- paste(vRoot, "Resultats/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Cc", Noise, "_", Method, ".csv", sep="")	
				MatrCc	<- array(unlist(fread(val,data.table=F)),dim=c(NBN,NBN))	# "Frank_TFxxx_TAyyy_z_Cc1_MRA.csv"
#				vMatrCc[]  	<- as.numeric(abs(MatrCc[]))		
				
				for (iSeuil in 1:NBSeuils) {
					MatrCc1	 <- Digitalize(MatrCc, "Threshold", Seuils[iSeuil])

#					if (Method %in% c("MRA", "MRA-CLR", "ARACNE")) {
#						MatrCc1	 <- Digitalize(MatrCc, "Top", Seuils[iSeuil])
#					} else {
#						MatrCc1	 <- Digitalize(MatrCc, "Threshold", Seuils[iSeuil])
#					}
#
					Score[ ,iSeuil,iMeth,iT,iTrial,iSet]	<- Benchmark(MatrCc1,Solut1)
				}	# Loop on iSeuil
					
				AUROC.tmp[iMeth,iT,iTrial,iSet]	 <- abs(trapz((1-Score["Sp", ,iMeth,iT,iTrial,iSet]), Score["Se", ,iMeth,iT,iTrial,iSet]))# AUC ROC for the plotted methods	
			}		# Loop on iMeth
		}			# Loop on iT
	}				# Loop on iTrial
}					# Loop on iSet	
cat("FIN Time ", as.character(Sys.time()), "\n")


for (iSet in 1:NBSets) {
	for (iT in 1:NBT) {
		for (iMeth in 1:NBM) {
			AUROC.tmp[iMeth,iT, ,iSet] 	<- ifelse(is.nan(AUROC.tmp[iMeth,iT, ,iSet]), 0.5, AUROC.tmp[iMeth,iT, ,iSet])
			AUROC[iMeth,iT,iSet]		<- mean(AUROC.tmp[iMeth,iT, ,iSet])
		}			# Loop on iMeth
	}				# Loop on iT
}					# Loop on iSet	


# AUROC curves drawing: new Figure 6

#	TA = 0		k = 0.5

Auroc5TA0_MRA.df		<- data.frame(NbNodes=NbNodes[1:8], auroc=AUROC["MRA",2,1:8],		Method="MRA")
# Auroc5TA0_MRACLR.df	<- data.frame(NbNodes=NbNodes[1:8], auroc=AUROC["MRA_CLR",2,1:8],	Method="MRA_CLR")
Auroc5TA0_LSE.df		<- data.frame(NbNodes=NbNodes[1:8], auroc=AUROC["LSE",2,1:8],		Method="TLR")
Auroc5TA0.df			<- rbind(Auroc5TA0_MRA.df, Auroc5TA0_LSE.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Auroc5TA0.pdf", sep="")
pdf(file=NomFic, height=6,width=3.5) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Auroc5TA0.df, aes(x=NbNodes, y=auroc, colour=factor(Method)), show.legend=FALSE) +		# legend removed to get more space (2 figs. with the same caption, added with Inkscape)
	xlab("Number of nodes") + ylab("AUROC") + ggtitle("Area under curve ROC, k = 50%, TA = 0") +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "TLR"="red")) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3.5 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()


#	TA = TF		k = 0.5

Auroc5TATF_MRA.df		<- data.frame(NbNodes=NbNodes[9:15], auroc=AUROC["MRA",2,9:15],		Method="MRA")
# Auroc5TATF_MRACLR.df	<- data.frame(NbNodes=NbNodes[9:15], auroc=AUROC["MRA_CLR",2,9:15],	Method="MRA_CLR")
Auroc5TATF_LSE.df		<- data.frame(NbNodes=NbNodes[9:15], auroc=AUROC["LSE",2,9:15],		Method="TLR")
Auroc5TATF.df			<- rbind(Auroc5TATF_MRA.df, Auroc5TATF_LSE.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Auroc5TATF.pdf", sep="")
pdf(file=NomFic, height=6,width=3.6) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Auroc5TATF.df, aes(x=NbNodes, y=auroc, colour=factor(Method)), show.legend=FALSE) +
	xlab("Number of nodes") + ylab("AUROC") + ggtitle("Area under curve ROC, k = 50%, TA = TF") +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "TLR"="red")) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3.5 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()


#
#	Sixth part
#	Se and Sp comparison (Reviewer 2 -- Major (1, 2 and 3))
#	Supplementary Information Table 6
#
NbNodes			<- vector(length=NBSets)				# Number of nodes

Noises			<- c(0.1, 0.5)							# Error coefficients
NBT				<- length(Noises)
Methods			<- c("MRA", "LSE_CI", "LASSO", "LSE", "STEP-Fo", "CLR", "MRA_CLR", "ARACNE", "MRNET")			# MRA-CLR
NBM				<- length(Methods)

ScoreTests		<-	c("TP", "TN", "FP", "FN", "Se", "Sp")		# Scores to measure (Se = "Sensibility", Sp = "Specificity", AUC = "AUROC")
NBSC			<- length(ScoreTests)					# Nbr. of scores tested
Score 			<- array(0, dim=c(NBSC, NBM, NBT, NBTrials, NBSets))
dimnames(Score)	<- list(ScoreTests, Methods, Noises, NULL, NULL)
ScoreM 			<- array(0, dim=c(2, NBM, NBT, NBSets))		# Average on the NBTrials files of the same size
dimnames(ScoreM)<- list(c("Se","Sp"), Methods, Noises, NULL)
FileScore		<- paste(vRoot, "Resultats/", "Frank_Scores.csv", sep="")
DistM 			<- array(0, dim=c(NBM, NBT, NBSets))	# Distance to the first bisecting line of the point (Se, 1-Sp)
dimnames(DistM) <- list(Methods, Noises, NULL)


for (iSet in 1:NBSets) {
	TF			<- as.numeric(FileNbr[iSet, 1])
	TA			<- as.numeric(FileNbr[iSet, 2])
	NBN			<- TF+TA
	NbNodes[iSet]	<- NBN
	
	Solution	<-	array(dim=c(NBN, NBN))				# "Real" value of the solution
	Solut1		<-	array(dim=c(NBN, NBN))				# Digitalized solution (0, 1)
	MatrCc		<-  array(dim=c(NBN, NBN))				# "r" matrix computed, according to the method used
	MatrCc1		<-  array(dim=c(NBN, NBN))				# Digitalized "r" matrix (0, 1)
		
	for (iTrial in 1:NBTrials) {
		val			<- paste(vRoot, "Donnees/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Sol.csv", sep="")
		Solution	<- as.matrix(read.csv(val))			# "Frank_TFxxx_TAyyy_z_Sol.csv"
		Solut1		<- Digitalize(Solution, "Top", 0.15)	# We keep the 15% highest absolute values

		for (iT in 1:NBT) {
			Noise	<- 10*Noises[iT]
			cat("TF ", TF, " TA ", TA, " Trial ", iTrial, " iT ", iT)
			
			for (iMeth in 1:NBM) {
				Method	<- Methods[iMeth]
				# cat("TF ", TF, " TA ", TA, " Trial ", iTrial, " iT ", iT, " Method ", Method, " Time ", as.character(Sys.time()))
				
				val		<- paste(vRoot, "Resultats/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Cc", Noise, "_", Method, ".csv", sep="")	
				if (! file.exists(val)) {
					Score[c("Se","Sp"),iMeth,iT,iTrial,iSet]	<- 0
					cat("   SKIP \n")
					next
				}
				
				MatrCc	<- array(unlist(fread(val,data.table=F)),dim=c(NBN,NBN))	# "Frank_TFxxx_TAyyy_z_Cc1_MRA.csv"

				if (Method == "MRA") {
					MatrCc1	<- Digitalize(MatrCc, "Threshold", 0.25)
				} else if (Method == "LSE_CI") {
					MatrCc1	<- Digitalize(MatrCc, "Threshold", 0)				
				} else if (Method == "LASSO") {
					MatrCc1	<- Digitalize(MatrCc, "Threshold", 0)
				} else if (Method == "LSE") {
					MatrCc1	<- Digitalize(MatrCc, "Threshold", 0.25)			# TLR method
				} else if (Method %in% c("STEP-Fo", "STEP-Ba")) {
					MatrCc1	<- Digitalize(MatrCc, "Threshold", 0)
				} else if (Method == "MRA_CLR") {
					MatrCc1	<- Digitalize(MatrCc, "Threshold", 0)
				} else if (Method %in% c("CLR", "ARACNE", "MRNET")) {				
					MatrCc1	<- Digitalize(MatrCc, "Threshold", 0)
				}

				Score[ ,iMeth,iT,iTrial,iSet]	<- Benchmark(MatrCc1,Solut1)
				cat(" Se ", Score["Se",iMeth,iT,iTrial,iSet], " Sp ", Score["Sp",iMeth,iT,iTrial,iSet], "\n")
			}		# Loop on iMeth
		}			# Loop on iT
	}				# Loop on iTrial
}					# Loop on iSet	
cat("FIN Time ", as.character(Sys.time()), "\n")

write.table(Noises[1], FileScore, sep=",")
write.table(Score[ , ,1, ,], FileScore, append=TRUE, sep=",")				# Frank_Scores.csv
write.table(Noises[2], FileScore, append=TRUE, sep=",")
write.table(Score[ , ,2, ,], FileScore, append=TRUE, sep=",")


for (iSet in 1:NBSets) {
	TF			<- as.numeric(FileNbr[iSet, 1])
	TA			<- as.numeric(FileNbr[iSet, 2])
	
	for (iT in 1:NBT) {
		for (iMeth in 1:NBM) {
			Method	<- Methods[iMeth]
			ScoreM["Se",iMeth,iT,iSet]	<- mean (Score["Se",iMeth,iT, ,iSet])
			ScoreM["Sp",iMeth,iT,iSet]	<- mean (Score["Sp",iMeth,iT, ,iSet])
			DistM[iMeth,iT,iSet]		<- -1*(1-ScoreM["Sp",iMeth,iT,iSet]) + 1*ScoreM["Se",iMeth,iT,iSet]	# dot product of (Se, 1-Sp) by (-1,1) = algebric distance to the first bisecting line
			cat("TF ", TF, " TA ", TA, " Noise ", Noises[iT], " Method ", Method, " Se ", ScoreM["Se",iMeth,iT,iSet], " Sp ", ScoreM["Sp",iMeth,iT,iSet], "\n")
		}			# Loop on iMeth
	}				# Loop on iT
}					# Loop on iSet			

write.table("Mean ", FileScore, append=TRUE, sep=",")
write.table(Noises[1], FileScore, append=TRUE, sep=",")
write.table(ScoreM[ , ,1, ], FileScore, append=TRUE, sep=",")				# Frank_Scores.csv
write.table(Noises[2], FileScore, append=TRUE, sep=",")
write.table(ScoreM[ , ,2, ], FileScore, append=TRUE, sep=",")

# DISTANCE curves drawing: new Figure 7

#	TA = 0		k = 0.1

Dist1TA0_MRA.df		<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["MRA",1,1:7],		Method="MRA")
Dist1TA0_LSECI.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["LSE_CI",1,1:7],		Method="LSE_CI")
Dist1TA0_LASSO.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["LASSO",1,1:7],		Method="LASSO")
Dist1TA0_LSE.df		<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["LSE",1,1:7],		Method="TLR")
Dist1TA0_STEPFo.df	<- data.frame(NbNodes=NbNodes[1:5], dist=DistM["STEP-Fo",1,1:5],	Method="STEP-Fo")
Dist1TA0_CLR.df		<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["CLR",1,1:7],		Method="CLR")
#Dist1TA0_MRACLR.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["MRA_CLR",1,1:7],	Method="MRA_CLR")
Dist1TA0_ARACNE.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["ARACNE",1,1:7],		Method="ARACNE")
Dist1TA0_MRNET.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["MRNET",1,1:7],		Method="MRNET")

Dist1TA0.df			<- rbind(Dist1TA0_MRA.df, Dist1TA0_LSECI.df, Dist1TA0_LASSO.df, Dist1TA0_LSE.df, Dist1TA0_STEPFo.df, Dist1TA0_CLR.df, Dist1TA0_ARACNE.df, Dist1TA0_MRNET.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Dist1TA0.pdf", sep="")
pdf(file=NomFic, height=6,width=3.5) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Dist1TA0.df, aes(x=NbNodes, y=dist, colour=factor(Method)), show.legend=FALSE) +		# legend removed to get more space (4 figs. with the same caption, added with Inkscape)
	xlab("Number of nodes") + ylab("Distance") + ggtitle("Dist. to the diagonal, k = 10%, TA = 0") +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "LSE_CI"="black", "LASSO"="green", "TLR"="red", "STEP-Fo"="orange", "CLR"="darkorchid2", "MRA_CLR"="lightpink3", 
												"ARACNE"="paleturquoise3", "MRNET"="grey75")) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3.5 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()


#	TA = 0		k = 0.5

Dist5TA0_MRA.df		<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["MRA",2,1:7],		Method="MRA")
Dist5TA0_LSECI.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["LSE_CI",2,1:7],		Method="LSE_CI")
Dist5TA0_LASSO.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["LASSO",2,1:7],		Method="LASSO")
Dist5TA0_LSE.df		<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["LSE",2,1:7],		Method="TLR")
Dist5TA0_STEPFo.df	<- data.frame(NbNodes=NbNodes[1:5], dist=DistM["STEP-Fo",2,1:5],	Method="STEP-Fo")
Dist5TA0_CLR.df		<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["CLR",2,1:7],		Method="CLR")
#Dist5TA0_MRACLR.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["MRA_CLR",2,1:7],	Method="MRA_CLR")
Dist5TA0_ARACNE.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["ARACNE",2,1:7],		Method="ARACNE")
Dist5TA0_MRNET.df	<- data.frame(NbNodes=NbNodes[1:7], dist=DistM["MRNET",2,1:7],		Method="MRNET")

Dist5TA0.df			<- rbind(Dist5TA0_MRA.df, Dist5TA0_LSECI.df, Dist5TA0_LASSO.df, Dist5TA0_LSE.df, Dist5TA0_STEPFo.df, Dist5TA0_CLR.df, Dist5TA0_ARACNE.df, Dist5TA0_MRNET.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Dist5TA0.pdf", sep="")
pdf(file=NomFic, height=6,width=3.5) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Dist5TA0.df, aes(x=NbNodes, y=dist, colour=factor(Method)), show.legend=FALSE) +
	xlab("Number of nodes") + ylab("Distance") + ggtitle("Dist. to the diagonal, k = 50%, TA = 0") +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "LSE_CI"="black", "LASSO"="green", "TLR"="red", "STEP-Fo"="orange", "CLR"="darkorchid2", "MRA_CLR"="lightpink3", 
												"ARACNE"="paleturquoise3", "MRNET"="grey75")) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3.5 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()

#	TA = TF		k = 0.1

Dist1TATF_MRA.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["MRA",1,9:13],		Method="MRA")
Dist1TATF_LSECI.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["LSE_CI",1,9:13],	Method="LSE_CI")
Dist1TATF_LASSO.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["LASSO",1,9:13],	Method="LASSO")
Dist1TATF_LSE.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["LSE",1,9:13],		Method="TLR")
Dist1TATF_STEPFo.df	<- data.frame(NbNodes=NbNodes[9:11], dist=DistM["STEP-Fo",1,9:11],	Method="STEP-Fo")
Dist1TATF_CLR.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["CLR",1,9:13],		Method="CLR")
# Dist1TATF_MRACLR.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["MRA_CLR",1,9:13],	Method="MRA_CLR")
Dist1TATF_ARACNE.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["ARACNE",1,9:13],	Method="ARACNE")
Dist1TATF_MRNET.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["MRNET",1,9:13],	Method="MRNET")

Dist1TATF.df		<- rbind(Dist1TATF_MRA.df, Dist1TATF_LSECI.df, Dist1TATF_LASSO.df, Dist1TATF_LSE.df, Dist1TATF_STEPFo.df, Dist1TATF_CLR.df, Dist1TATF_ARACNE.df, Dist1TATF_MRNET.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Dist1TATF.pdf", sep="")
pdf(file=NomFic, height=6,width=3.5) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Dist1TATF.df, aes(x=NbNodes, y=dist, colour=factor(Method)), show.legend=FALSE) +
	xlab("Number of nodes") + ylab("Distance") + ggtitle("Dist. to the diagonal, k = 10%, TA = TF") +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "LSE_CI"="black", "LASSO"="green", "TLR"="red", "STEP-Fo"="orange", "CLR"="darkorchid2", "MRA_CLR"="lightpink3", 
												"ARACNE"="paleturquoise3", "MRNET"="grey75")) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3.5 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()

#	TA = TF		k = 0.5

Dist5TATF_MRA.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["MRA",2,9:13],		Method="MRA")
Dist5TATF_LSECI.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["LSE_CI",2,9:13],	Method="LSE_CI")
Dist5TATF_LASSO.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["LASSO",2,9:13],	Method="LASSO")
Dist5TATF_LSE.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["LSE",2,9:13],		Method="TLR")
Dist5TATF_STEPFo.df	<- data.frame(NbNodes=NbNodes[9:11], dist=DistM["STEP-Fo",2,9:11],	Method="STEP-Fo")
Dist5TATF_CLR.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["CLR",2,9:13],		Method="CLR")
# Dist5TATF_MRACLR.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["MRA_CLR",2,9:13],	Method="MRA_CLR")
Dist5TATF_ARACNE.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["ARACNE",2,9:13],	Method="ARACNE")
Dist5TATF_MRNET.df	<- data.frame(NbNodes=NbNodes[9:13], dist=DistM["MRNET",2,9:13],	Method="MRNET")

Dist5TATF.df		<- rbind(Dist5TATF_MRA.df, Dist5TATF_LSECI.df, Dist5TATF_LASSO.df, Dist5TATF_LSE.df, Dist5TATF_STEPFo.df, Dist5TATF_CLR.df, Dist5TATF_ARACNE.df, Dist5TATF_MRNET.df)

NomFic 	<- paste(vRoot, "Resultats/", "FRANK_Dist5TATF.pdf", sep="")
pdf(file=NomFic, height=6,width=3.5) 							# Sizes default to 7

ggplot() + 
	geom_line(data=Dist5TATF.df, aes(x=NbNodes, y=dist, colour=factor(Method)), show.legend=FALSE) +
	xlab("Number of nodes") + ylab("Distance") + ggtitle("Dist. to the diagonal, k = 50%, TA = TF") +
	scale_colour_manual(name="Method", values=c("MRA"="blue", "LSE_CI"="black", "LASSO"="green", "TLR"="red", "STEP-Fo"="orange", "CLR"="darkorchid2", "MRA_CLR"="lightpink3", 
												"ARACNE"="paleturquoise3", "MRNET"="grey75")) +
	theme(axis.text.x = element_text(color="black")) +
	theme(axis.text.y = element_text(color="black")) +
	theme(text = element_text(size = 10))					# width = 3.5 (4x8 cm) vs size = 13 width = 6 (8x8 cm)
	
dev.off()


#
#	Utilities
#

##################################################################################
#
#	Use of FRANK, a gene network simulator  (Reviewer 2 -- Major)
#	https://m2sb.org/?page=FRANK
#
#	The third part of the previous program is used to detect networks structure and compute error, but it's very time consuming.
#	In case we have saved the files "Frank_TFxxx_TAyyy_z_Sol.csv" and "Frank_TFxxx_TAyyy_z_Cc.csv", this utility computes the error (array ErrorQ)
#

vRoot		<- c("/home/jpborg/Documents/Publications/BioInformatics/MRA_Regression_ed1/Vers_Rev1/Donnees/FRANK2/")

#Methods		<- c("MRA", "LSE", "TLR", "STEP-Fo", "ARACNE")
#Methods		<- c("MRA", "MRA-CLR", "LSE", "STEP-Fo", "ARACNE")			# stocké ainsi
Methods		<- c("MRA", "MRA_CLR", "ARACNE", "LSE",  "STEP-Fo")
NBM			<- length(Methods)							# Nbr of methods to test

ErrorQ		<- array(0, dim=c(NBM, NBT, NBTrials, NBSets))				# Quadratic error, as a function of noise level and method used
rownames(ErrorQ)	<- c(Methods)
FileErrorQ	<- paste(vRoot, "Resultats/", "Frank_ErrorQ3.csv", sep="")

cat("DEBUT Time ", as.character(Sys.time()), "\n")

for (iSet in 1:NBSets) {
	TF			<- as.numeric(FileNbr[iSet, 1])
	TA			<- as.numeric(FileNbr[iSet, 2])
	NBN			<- TF+TA

	for (iTrial in 1:NBTrials) {
		cat("iSet ", iSet, " TF ", TF, " TA ", TA, " iTrial ", iTrial, "\n")

		val			<- paste(vRoot, "Donnees/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Sol.csv", sep="")
		Solution	<- as.matrix(read.csv(val))				# "Frank_TFxxx_TAyyy_z_Sol.csv"

		for (iT in 1:NBT) {
			for (iMeth in 1:NBM) {
				Method 	<- Methods[iMeth]
				cat(" iT ", iT, " Method ", Method, " Time ", as.character(Sys.time()), "\n")
				
				val		<- paste(vRoot, "Resultats/", "Frank_TF", TF, "_TA", TA, "_", iTrial, "_Cc", 10*Noises[iT], "_", Method, ".csv", sep="")
				MatrCc		<- array(unlist(fread(val,data.table=F)),dim=c(NBN,NBN))	
				
				ss		<- sum((MatrCc-Solution)**2)
				ErrorQ[iMeth,iT,iTrial,iSet]	<- sqrt(ss)
			
				val		<- paste("\n\n\nFrank_TF", TF, "_TA", TA, "_", iTrial, "_Cc", 10*Noises[iT], "_", Method, "\n", sep="")
				write.table(val, FileErrorQ, append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
				write.table(ErrorQ, FileErrorQ, append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")	# Frank_ErrorQ.csv
			}		# Loop on iMeth		
		}			# Loop on iT
	}				# Loop on iTrial
}					# Loop on iSet

cat("FIN Time ", as.character(Sys.time()), "\n")

RightM1		<- array(0, dim=c(2,15))
rownames(RightM1)=c("mean", "sd")
for (iSet in 9:11) {
	RightM1["mean",iSet] 	<- round (mean(ErrorQ["STEP-Fo",1, ,iSet]), digits=3)
	RightM1["sd",iSet] 		<- round ((var(ErrorQ["STEP-Fo",1, ,iSet]))**0.5, digits=3)
}			# Loop on iSet

##################################################################################
