#
#	Readme : test
#	R software used to write the article to publish. 
#	The files to download concerning figures 4, 5 and SI, along with their URL, are indicated in the code.
#


#
#	Réinitialisation de la session R et des variables d'environnement, pour repartir à 0
#
rm(list = ls())		# Il faut aussi faire Ctrl + Shift + F10, qui décharge les packages en relançant une session R.
					# On peut aussi détacher chaque package individuellement par : detach(package:packagename)

#
# Déclaration des librairies et des fonctions utilisées dans le document
#

library ("deSolve")
library ("rootSolve")
library ("Deriv")
library ("stringr")
library	("forcats")			# Pour la fonction fct_relevel
library ("ggplot2")
library ("glmnet")
library ("data.table")		# Pour la fonction fread
library ("dplyr")			# Pour utiliser rename
library ("PRROC")			# Pour calcul des AUC
library ("ROCR")			# Inutile, chargé automatiquement au démarrage d'une session R
library ("pROC")			# Pour calcul des AUC
library ("verification")	# Pour roc.area
library ("pracma")			# pour trapz
library	("lars")

F <- function(X)
   c(F1 = (KC1*U*X[1]) 	/ ((K11+X[1]+(MKKK0-X[1]-X[2])*K11/K12)*(1+X[6]/Ki)) - ((V4*(MKKK0-X[1]-X[2]))  / (K31+X[2]+(MKKK0-X[1]-X[2])*K31/K32+X[1]*K31/K33)),
	 F2 = (KC2*U*(MKKK0-X[1]-X[2]))    / ((K11+X[1]+(MKKK0-X[1]-X[2])*K11/K12)*(1+X[6]/Ki)) - ((V3*X[2])  / (K31+X[2]+(MKKK0-X[1]-X[2])*K31/K32+X[1]*K31/K33)),
	 F3 = (KC5*X[3]*X[2]) / (K51+X[3]+(MKK0-X[3]-X[4])*K51/K52) - ((V8*(MKK0-X[3]-X[4])*(1+A*X[6]/Kmp))   / ((K71+X[4]+(MKK0-X[3]-X[4])*K71/K72+X[3]*K71/K73)*(1+X[6]/Kmp))),
	 F4 = (KC6*(MKK0-X[3]-X[4])*X[2])  / (K51+X[3]+(MKK0-X[3]-X[4])*K51/K52) - ((V7*X[4]*(1+A*X[6]/Kmp))  / ((K71+X[4]+(MKK0-X[3]-X[4])*K71/K72+X[3]*K71/K73)*(1+X[6]/Kmp))),
	 F5 = (KC9*X[4]*X[5]) / (K91+X[5]+(MAPK0-X[5]-X[6])*K91/K92) - ((V12*(MAPK0-X[5]-X[6])) / (K111+X[6]+(MAPK0-X[5]-X[6])*K111/K112+X[5]*K111/K113)),
	 F6 = (KC10*X[4]*(MAPK0-X[5]-X[6]))/ (K91+X[5]*(MAPK0-X[5]-X[6])*K91/K92) - ((V11*X[6]) / (K111+X[6]+(MAPK0-X[5]-X[6])*K111/K112+X[5]*K111/K113)))

FF <- function(X1,X2,X3,X4,X5,X6)
   c(F1 = (KC1*U*X1) 	/ ((K11+X1+(MKKK0-X1-X2)*K11/K12)*(1+X6/Ki)) - ((V4*(MKKK0-X1-X2))  / (K31+X2+(MKKK0-X1-X2)*K31/K32+X1*K31/K33)),
	 F2 = (KC2*U*(MKKK0-X1-X2))    / ((K11+X1+(MKKK0-X1-X2)*K11/K12)*(1+X6/Ki)) - ((V3*X2)  / (K31+X2+(MKKK0-X1-X2)*K31/K32+X1*K31/K33)),
	 F3 = (KC5*X3*X2) / (K51+X3+(MKK0-X3-X4)*K51/K52) - ((V8*(MKK0-X3-X4)*(1+A*X6/Kmp))   / ((K71+X4+(MKK0-X3-X4)*K71/K72+X3*K71/K73)*(1+X6/Kmp))),
	 F4 = (KC6*(MKK0-X3-X4)*X2)  / (K51+X3+(MKK0-X3-X4)*K51/K52) - ((V7*X4*(1+A*X6/Kmp))  / ((K71+X4+(MKK0-X3-X4)*K71/K72+X3*K71/K73)*(1+X6/Kmp))),
	 F5 = (KC9*X4*X5) / (K91+X5+(MAPK0-X5-X6)*K91/K92) - ((V12*(MAPK0-X5-X6)) / (K111+X6+(MAPK0-X5-X6)*K111/K112+X5*K111/K113)),
	 F6 = (KC10*X4*(MAPK0-X5-X6))/ (K91+X5*(MAPK0-X5-X6)*K91/K92) - ((V11*X6) / (K111+X6+(MAPK0-X5-X6)*K111/K112+X5*K111/K113)))

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


#	Intialisation des calculs
#
st <- c(200,  0, 100,  50, 100,  0)

#
#	Figure 1 (dessin du réseau de kinases étudié)
#

# Constantes liées à la cinétique des réactions
KC1	 <- 1			# Taux catalytiques en s-1
KC2	 <- 15
KC5	 <- 1
KC6	 <- 15
KC9	 <- 1
KC10 <- 15

K11	 <- 300			# Constantes de Michaelis en nM
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

A	 <- 5			# Coeff. sans dimension

V3	 <- 18.8		# Taux d'enzyme maximal en nM.s-1
V4	 <- 16.4
V7	 <- 18.8
V8	 <- 16.4
V11	 <- 8.4
V12	 <- 7.3

MKKK0	<- 200		# X1+X7+X2	: concentration totale de la protéine MKKK
MKK0	<- 180		# X3+X8+X4	: concentration totale de la protéine MKK
MAPK0	<- 360		# X5+X9+X6	: concentration totale de la protéine MAPK
U	<- 20			# Injection de smallGTPase Ras à partir de l'instant 0+


# 1/ Recherche des concentrations à l'état stationnaire (dXi/dt = 0)

ss <- multiroot (f=F, start= st)	# utilise le package rootSolve
XStat <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])		# Concentrations théoriques à l'état stationnaire
XMesNPMoy	<- 	mean(XStat)			# l'écart-type du bruit sera défini par rapport à ce niveau moyen,
									# concentration mesurée dans les 6 noeuds, réseau non perturbé, état stationnaire.

# 2/ Vérifier que l'état stationnaire est atteint pour t = 1000 s

Ini	<- c(X1	 = 199.3489,
		 X2	 = 0.0334,
		 X3	 = 179.8226,
		 X4	 = 0.0091,
		 X5	 = 358.2074,
		 X6	 = 0.1519)
		 
#	Comportement à long terme : recherche de l'état stationnaire
Times	<- seq(0, 1001, by=1)		# temps en s.
Parms	<- 0						# Pas de perturbation

out	<-	ode(Ini, Times, RevSigNet, Parms)
out[1001,]

# 3/ Calcul des dérivées partielles des fi pour t = 1000s
 
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

# 4/ Calcul de la matrice r théorique

matrTh	<- matrix(, nrow=6, ncol=6)

matrTh[1,]	<- -c(X1*dF1dX1, X2*dF1dX2, X3*dF1dX3, X4*dF1dX4, X5*dF1dX5, X6*dF1dX6) / (X1*dF1dX1)
matrTh[2,]	<- -c(X1*dF2dX1, X2*dF2dX2, X3*dF2dX3, X4*dF2dX4, X5*dF2dX5, X6*dF2dX6) / (X2*dF2dX2)
matrTh[3,]	<- -c(X1*dF3dX1, X2*dF3dX2, X3*dF3dX3, X4*dF3dX4, X5*dF3dX5, X6*dF3dX6) / (X3*dF3dX3)
matrTh[4,]	<- -c(X1*dF4dX1, X2*dF4dX2, X3*dF4dX3, X4*dF4dX4, X5*dF4dX5, X6*dF4dX6) / (X4*dF4dX4)
matrTh[5,]	<- -c(X1*dF5dX1, X2*dF5dX2, X3*dF5dX3, X4*dF5dX4, X5*dF5dX5, X6*dF5dX6) / (X5*dF5dX5)
matrTh[6,]	<- -c(X1*dF6dX1, X2*dF6dX2, X3*dF6dX3, X4*dF6dX4, X5*dF6dX5, X6*dF6dX6) / (X6*dF6dX6)
#	Cette matrice est utilisée pour dessiner le réseau par Cytoscape, v3.9.1 (fic MRA3.csv généré manuellement)
#	et pour la matrice "B" (connectivity matrix, exact value)

#	Pour le dessin des deux autres matrices à droite, voir fig. 2 : MatrCc (iTest = 4 : matrice "A", -50% et iTest = 9 : matrice "C", +50%)




#
#	Figure 2  (Influence de l'amplitude des perturbations et du bruit)
#
	
Perturb	<- c(0.2, 0.3, 0.4, 0.5, 0.9, 0.99, 1.01, 1.10, 1.50, 1.6, 1.7, 1.8)	# 12 perturbations
NBN 	<- 6
NBT		<- length(Perturb)
NBK 	<- 1								# Nombre de paquets de mesures. 
		# Un paquet correspond à un ensemble de NBN mesures correspondant à NBN perturbations portant chacune sur chaque noeud indépendamment.
			
MatDif	<- matrix(nrow=2, ncol= NBN*NBN)	# Pour mesurer la distance euclidienne entre deux matrices
XMesNP	<- array(dim=c(NBN, NBK, NBT))		# Valeurs mesurées, système non perturbé
XMesP	<- array(dim=c(NBN, NBN, NBK, NBT))	# Valeurs mesurées, suite à une perturbation
MatrCc	<- array(-1, dim=c(NBN, NBN, NBT))	# Matrice r calculée
Err		<- matrix(nrow=NBT, ncol=1)			# Erreur correspondante (distance euclidienne vs. matrice théorique)

for (iTest in c(1:NBT)) {
	#	Entrée des "données de mesure"
	ss <- multiroot (f=F, start= st)	# Mesure non perturbée  - utilise le package rootSolve
	XMesNP[,1,iTest]	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
										# Dans cet essai, ce calcul aurait pu être sorti de la boucle, mais on le conserve ici pour ne pas modifier l'algorithme
	
	NPerturb	<- Perturb[iTest]		# Niveau de perturbation
	V4 <- NPerturb*V4		# P1
	ss <- multiroot (f=F, start= st)
	XMesP[,1,1,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V4 <- V4/NPerturb		# Retour à la valeur initiale
	
	V3 <- NPerturb*V3		# P2
	ss <- multiroot (f=F, start= st)
	XMesP[,2,1,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V3 <- V3/NPerturb		# Retour à la valeur initiale
	
	V8 <- NPerturb*V8		# P3
	ss <- multiroot (f=F, start= st)
	XMesP[,3,1,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V8 <- V8/NPerturb		# Retour à la valeur initiale
	
	V7 <- NPerturb*V7		# P4
	ss <- multiroot (f=F, start= st)
	XMesP[,4,1,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V7 <- V7/NPerturb		# Retour à la valeur initiale
	
	V12 <- NPerturb*V12		# P5
	ss <- multiroot (f=F, start= st)
	XMesP[,5,1,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V12 <- V12/NPerturb		# Retour à la valeur initiale
	
	V11 <- NPerturb*V11		# P6
	ss <- multiroot (f=F, start= st)
	XMesP[,6,1,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6])
	V11 <- V11/NPerturb		# Retour à la valeur initiale

	#	Calcul de la matrice "r"
	for (iPaq in c(1:NBK)) {
		if(iPaq ==1) {
			MatR  <- 2 * (XMesP[,,iPaq,iTest]-XMesNP[,iPaq,iTest]) / (XMesP[,,iPaq,iTest]+XMesNP[,iPaq,iTest])
		} else {
			MatR  <- rbind(MatR[,], 2 * (XMesP[,,iPaq,iTest]-XMesNP[,iPaq,iTest]) / (XMesP[,,iPaq,iTest]+XMesNP[,iPaq,iTest]))
			#	Dans ce cas, la taille de MatR augmente : [NBN*NBK, NBN] en finale
		}
	}

	for (i in c(1:NBN)) {
		col <- c(1:NBN)
		col <- col[-i]
		
		CY	<- MatR[i,]
		CY	<- CY[-i]				# on enlève la colonne i
		MatR1T	<- MatR
		MatR1T	<- MatR1T[-i,]		# on enlève la ligne i et la colonne i
		MatR1T	<- MatR1T[,-i]
		MatR1T	<- t(MatR1T)		# on transpose la matrice
		MatR1T	<- solve(MatR1T)	# on inverse la matrice transposée
		CX	<- MatR1T %*% CY
		
		MatrCc[i, col, iTest] <- CX
	}
	
	#	Calcul de l'erreur
	for (i in c(1:NBN)) {
		col <- c(((i-1)*NBN +1) : (i*NBN))
		MatDif[1,col] <- matrTh[i,]
		MatDif[2,col] <- MatrCc[i,,iTest]
	}
	Err[iTest, 1] <-  dist(MatDif)
	cat ("iTest", iTest, " Perturb ", NPerturb, " Err ", Err[iTest, 1], "\n")
}

png(filename="MRA3_2_DP2.png")	# fichier png contenant les dessins (influence du niveau de perturbation - en anglais pour publication)
plot (Perturb, Err, xlim=c(0.3, 1.8), ylim=c(0,4), main="Error Cc vs. Th", type="l", col="blue", xlab="Perturbation", ylab="Error")
dev.off()

# Plot de la courbe de gauche (complétée par Inkscape pour tracer les points A, B et C, avec leurs coordonnées)


NBN 	<- 6								# Nombre de noeuds
Noise	<-	c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.05, 0.1, 0.5)	# Niveaux de bruit (coeff. k)
Perturb <-  c(0.99, 1.01, 0.7, 0.9, 1.1, 1.3, 1.5, 0.5)		# -50% ajouté à la fin pour le nouveau dessin de la fig 3, sans toucher aux indices des dessins anciens
PercentPert	<- c("-1", "+1", "-30", "-10", "+10", "+30", "+50", "-50")		# Pourcentages de perturbation

NBT 	<- length(Noise)
NBP		<- length(Perturb)
NBK		<- 100								# 100 tirages de bruit

Form	<- vector(length=NBN)				# Formule décrivant la régression
XMesNP	<- array(dim=c(NBN, NBK, NBT))		# Valeur mesurées, système non perturbé
XMesP	<- array(dim=c(NBN, NBN, NBK, NBT))	# Valeur mesurées, suite à une perturbation
MatR	<- array(dim=c(NBN, NBN, NBK))		# Matrice R
MatrCc	<- array(dim=c(NBN, NBN, NBK, NBT))	# Matrice r calculée

MatDif	<- matrix(nrow=2, ncol= NBN*NBN)	# Pour mesurer la distance euclidienne entre deux matrices
Err		<- array(dim=c(NBK, NBT))			# Erreur correspondante (distance euclidienne entre la matrice "r" calculée et la matrice théorique)
ErrLn	<- array(dim=c(NBK, NBT))			# Ln(1+Err)

rMoy	<- array(dim=c(NBN, NBN, NBT))		# Coefficient moyen
rSd 	<- array(dim=c(NBN, NBN, NBT))		# Ecart-type
ErrMoy	<- array(dim=c(NBT, NBP))			# Erreur moyenne
ErrSd	<- array(dim=c(NBT, NBP))			# Ecart-type
ErrMoyLn<- array(dim=c(NBT, NBP))			# Erreur moyenne de ErrLn
ErrSdLn	<- array(dim=c(NBT, NBP))			# Ecart-type de ErrLn

# Recherche de la formule

for (i in 1:NBN) {
	col <- 1:NBN
	col <- col[-i]
	
	Form[i] <- paste("MatR[", i, ",col,iPaq]-", sep="")			# Formule décrivant la régression
	for (j in col) {
		Form[i] <- paste(Form[i], "MatR[", sep="A")
		Form[i] <- paste(Form[i], j, sep="")
		Form[i] <- paste(Form[i], ",col", sep="")
		Form[i] <- paste(Form[i], ",iPaq]", sep= "")
	}
	
	Form[i] <- paste(Form[i], "+0", sep="")
	Form[i] <- str_replace_all(Form[i], pattern = "-A", replacement = "~")
	Form[i] <- str_replace_all(Form[i], pattern = "A",  replacement = "+")
	# cat("i ", i, " Form ", Form, "\n")
}

#	Calcul de "r" par régression linéaire portant sur UNE expérience (NBK matrices r calculées pour chaque k et chaque perturbation)

for (iPerturb in c(1:NBP)) {
	NPerturb <- Perturb[iPerturb]
	set.seed(12345)								# Pour pouvoir régénérer les mêmes séquences
	
	for (iTest in c(1:NBT)) {
		# Pour chaque niveau de bruit (iTest), on effectue 100 tirages (NBK) et on calcule 100 matrices "r" (chaque matrice correspond à un tirage)
		# Puis on effectue des statistiques sur les erreurs obtenues : moyenne et sd des ri,j, ainsi que l'erreur moyenne et son écart-type
		
		sd <- Noise[iTest]*XMesNPMoy			# écart-type du bruit
	
		for (iPaq in c(1:NBK)) {
			#	Entrée des "données de mesure" bruitées
			ss <- multiroot (f=F, start= st)	# Mesure non perturbée  - utilise le package rootSolve
			XMesNP[,iPaq,iTest]	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			
			V4 <- NPerturb*V4		# P1
			ss <- multiroot (f=F, start= st)
			XMesP[,1,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V4 <- V4/NPerturb		# Retour à la valeur initiale
			
			V3 <- NPerturb*V3		# P2
			ss <- multiroot (f=F, start= st)
			XMesP[,2,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V3 <- V3/NPerturb		# Retour à la valeur initiale
			
			V8 <- NPerturb*V8		# P3
			ss <- multiroot (f=F, start= st)
			XMesP[,3,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V8 <- V8/NPerturb		# Retour à la valeur initiale
			
			V7 <- NPerturb*V7		# P4
			ss <- multiroot (f=F, start= st)
			XMesP[,4,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V7 <- V7/NPerturb		# Retour à la valeur initiale
			
			V12 <- NPerturb*V12		# P5
			ss <- multiroot (f=F, start= st)
			XMesP[,5,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V12 <- V12/NPerturb		# Retour à la valeur initiale
			
			V11 <- NPerturb*V11		# P6
			ss <- multiroot (f=F, start= st)
			XMesP[,6,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
			V11 <- V11/NPerturb		# Retour à la valeur initiale	
		
			for (iPert in c(1:NBN)) {
				MatR[,iPert,iPaq]  <- 2 * (XMesP[,iPert,iPaq,iTest]-XMesNP[,iPaq,iTest]) / (XMesP[,iPert,iPaq,iTest]+XMesNP[,iPaq,iTest])
			}
			
			for (i in 1:NBN) {
				col <- 1:NBN
				col <- col[-i]
						
				MatrCc[i, col, iPaq, iTest]   <- (lm(as.formula(Form[i])))$coefficients
				MatrCc[i, i,   iPaq, iTest]   <- -1
				
				col2 <- c(((i-1)*NBN +1) : (i*NBN))
				MatDif[1,col2] <- matrTh[i,]
				MatDif[2,col2] <- MatrCc[i,, iPaq, iTest]
				
				Err[iPaq, iTest] 	<-  dist(MatDif)	# Err : distance euclidienne entre MatrCc ("r" calculée) et matrTh ("r" théorique)
				ErrLn[iPaq, iTest]	<- 	log1p(Err[iPaq, iTest])		# Ln(1+Err).  La fonction log1p(x) = Ln(1+x) est plus précise (parait-il) que log(1+x) pour |x| << 1
			}		
		}	# Tirages de bruit (NBK)
		
		for (i in 1:NBN) {
			for (j in 1:NBN) {
				rMoy[i, j, iTest] <- mean(MatrCc[i, j, , iTest])
				rSd [i, j, iTest] <- sd  (MatrCc[i, j, , iTest])
			}
		}
		
		ErrMoy  [iTest, iPerturb]	<- mean(Err  [, iTest])
		ErrSd   [iTest, iPerturb]	<- sd  (Err  [, iTest])
		ErrMoyLn[iTest, iPerturb]	<- mean(ErrLn[, iTest])
		ErrSdLn [iTest, iPerturb]	<- sd  (ErrLn[, iTest])
		cat("Perturb ", Perturb[iPerturb], " Coef ", Noise[iTest], " ErrMoy ", ErrMoy[iTest, iPerturb], " ErrSd ", ErrSd[iTest, iPerturb], "\n")
	}		# Boucle sur les niveaux de bruit (NBT)
}			# Boucle sur les perturbations (NBP)

vP		<- c(8, 3, 4, 1,2, 5, 6, 7)		# Perturbations (-50, -30, -10, -1, +1, +10, +30, +50)  --  Pour ne pas modifier l'ordre de Perturb (utilisé dans les premiers dessins fig 3)
matrN2.df	<- data.frame(Pert=1:8, mean=ErrMoyLn[2, vP],  lower=ErrMoyLn[2, vP]-ErrSdLn[2, vP],   upper=ErrMoyLn[2, vP]+ErrSdLn[2, vP])
matrN5.df	<- data.frame(Pert=1:8, mean=ErrMoyLn[5, vP],  lower=ErrMoyLn[5, vP]-ErrSdLn[5, vP],   upper=ErrMoyLn[5, vP]+ErrSdLn[5, vP])
matrN10.df	<- data.frame(Pert=1:8, mean=ErrMoyLn[10, vP], lower=ErrMoyLn[10, vP]-ErrSdLn[10, vP], upper=ErrMoyLn[10, vP]+ErrSdLn[10, vP])

Titre_en	<- "Ln (1+Err) based on perturbation and noise"
png(filename="MRA3_3_ne.png")	# fichier png : Ln(1+Err moyen) en f° de la perturbation, pour 3 niveaux de bruit (moyenne et +- sd de Ln(1+Err)
								# k=0,002  k=0,005  k=0,01
								# Titre en anglais pour publication
ggplot() + geom_linerange(data=matrN2.df,   aes(x=Pert, ymin=0, ymax=mean), size=3, color="blue")+ 
		   geom_errorbar (data=matrN2.df,   aes(x=Pert, ymin=lower, ymax=upper), color = "black", width = 0.05)+
		   geom_linerange(data=matrN5.df,   aes(x=Pert+0.15, ymin=0, ymax=mean), size=3, color="yellow")+ 
		   geom_errorbar (data=matrN5.df,   aes(x=Pert+0.15, ymin=lower, ymax=upper), color = "black", width = 0.05)+
		   geom_linerange(data=matrN10.df,  aes(x=Pert+0.3, ymin=0, ymax=mean), size=3, color="red")+ 
		   geom_errorbar (data=matrN10.df,  aes(x=Pert+0.3, ymin=lower, ymax=upper), color = "black", width = 0.05)+
		   labs(title=Titre_en, x="Perturbation (%)", y="Ln (1+Err)") + scale_x_continuous(breaks=1:8+0.15, labels=PercentPert[vP])	   
dev.off()

# Plot de la courbe de droite, complétée par "Legende.svg" créé par Inkscape.




#
#	Figure 3 : Sensibility and specificity of regression methods
#
NBN 	<- 6								# Nombre de noeuds

Noise		<-	c(0.001, 0.005)
Replic		<-  c(3, 5)
NPerturb 	<-  1.5

NBT 	<- length(Noise)
NBK		<- 100										# 100 tirages de bruit
NBR		<- length(Replic)							# Les essais de réplicats
NBRM	<- max(Replic)								# Nb max de réplications à effectuer
		
Form	<- vector(length=NBN)						# Formule décrivant la régression pour le calcul de MatrCc
XMesNP	<- array(dim=c(NBN, NBK, NBT))				# Valeur mesurées, système non perturbé
XMesP	<- array(dim=c(NBN, NBN, NBK, NBT))			# Valeur mesurées, suite à une perturbation
MatR	<- array(dim=c(NBN, NBN, NBK))				# Matrice R
MatrCc	<- array(dim=c(NBN, NBN, NBK, NBT))			# Matrice r calculée
		
pX 	<- array(dim=c(NBN-1, NBN-1, NBN))				# ligne (iRow), colonne (iCol), n° matrice (iMat)
pY 	<- array(dim=c(NBN-1, NBN))						# ligne, n° matrice
gX 	<- array(dim=c(NBRM*(NBN-1), NBN-1, NBN))		# ligne (iRow), colonne (iCol), n° matrice (iMat)
gY 	<- array(dim=c(NBRM*(NBN-1), NBN))				# ligne, 1, n° matrice
gX1	<- array(dim=c(NBRM*(NBN-1), NBN-1, NBN))		# Copie de gX
gY1	<- array(dim=c(NBRM*(NBN-1), NBN))				# Copie de gY
		
VMean	<- array(dim=c(NBN, NBN, NBR, NBT))			# Moyenne calculée
VMin	<- array(dim=c(NBN, NBN, NBR, NBT))			# Valeur mini. de l'intervalle de confiance à 95%
VMax	<- array(dim=c(NBN, NBN, NBR, NBT))			# Valeur maxi. de l'intervalle de confiance à 95%

ValNom		<- vector(length=(NBN-1)*NBN)			# Nom (n) correspondant à (i,j)
ValVrai		<- vector(length=(NBN-1)*NBN)			# Vraie valeur de  r(i,j)
ValMean		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Moyenne des réplicats
ValMin		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Min IC95 des réplicats (correspond à p_value à 5%)
ValMax		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Max IC95 des réplicats (correspond à p_value à 5%)
ValMeanBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Moyenne des 100 valeurs (cf Boxplot)
ValMinBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Quantile à 2,5%  des 100 valeurs
ValMaxBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Quantile à 97,5% des 100 valeurs

# Recherche de la formule décrivant la régression pour le calcul de MatrCc

for (i in 1:NBN) {
	col <- 1:NBN
	col <- col[-i]
	
	Form[i] <- paste("MatR[", i, ",col,iPaq]-", sep="")			# Formule décrivant la régression
	for (j in col) {
		Form[i] <- paste(Form[i], "MatR[", sep="A")
		Form[i] <- paste(Form[i], j, sep="")
		Form[i] <- paste(Form[i], ",col", sep="")
		Form[i] <- paste(Form[i], ",iPaq]", sep= "")
	}
	
	Form[i] <- paste(Form[i], "+0", sep="")
	Form[i] <- str_replace_all(Form[i], pattern = "-A", replacement = "~")
	Form[i] <- str_replace_all(Form[i], pattern = "A", replacement = "+")
	#	cat("BOOTSTRAP i ", i, " Form ", Form[i], "\n")
}

set.seed(12345)								# Pour pouvoir régénérer les mêmes séquences

for (iTest in c(1:NBT)) {
	# Pour chaque niveau de bruit (iTest), on effectue 100 tirages (NBK) et on calcule 100 matrices "r" (chaque matrice correspond à un tirage)
	# Puis on effectue des statistiques sur les erreurs obtenues : moyenne et sd des ri,j, ainsi que l'erreur moyenne et son écart-type
	# La régression linéaire porte sur une valeur seule. Cette méthode (6 matrices 5*5) est équivalente à l'inversion de la matrice 6*6
	
	sd <- Noise[iTest]*XMesNPMoy			# écart-type du bruit

	for (iPaq in c(1:NBK)) {
		#	Entrée des "données de mesure" bruitées
		ss <- multiroot (f=F, start= st)	# Mesure non perturbée  - utilise le package rootSolve
		XMesNP[,iPaq,iTest]	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		
		V4 <- NPerturb*V4		# P1
		ss <- multiroot (f=F, start= st)
		XMesP[,1,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V4 <- V4/NPerturb		# Retour à la valeur initiale
		
		V3 <- NPerturb*V3		# P2
		ss <- multiroot (f=F, start= st)
		XMesP[,2,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V3 <- V3/NPerturb		# Retour à la valeur initiale
		
		V8 <- NPerturb*V8		# P3
		ss <- multiroot (f=F, start= st)
		XMesP[,3,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V8 <- V8/NPerturb		# Retour à la valeur initiale
		
		V7 <- NPerturb*V7		# P4
		ss <- multiroot (f=F, start= st)
		XMesP[,4,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V7 <- V7/NPerturb		# Retour à la valeur initiale
		
		V12 <- NPerturb*V12		# P5
		ss <- multiroot (f=F, start= st)
		XMesP[,5,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V12 <- V12/NPerturb		# Retour à la valeur initiale
		
		V11 <- NPerturb*V11		# P6
		ss <- multiroot (f=F, start= st)
		XMesP[,6,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V11 <- V11/NPerturb		# Retour à la valeur initiale	
	
		for (iMat in c(1:NBN)) {
			MatR[,iMat,iPaq]  <- 2 * (XMesP[,iMat,iPaq,iTest]-XMesNP[,iPaq,iTest]) / (XMesP[,iMat,iPaq,iTest]+XMesNP[,iPaq,iTest])
		}
		
		for (i in 1:NBN) {
			col <- 1:NBN
			col <- col[-i]
					
			MatrCc[i, col, iPaq, iTest]   <- (lm(as.formula(Form[i])))$coefficients
			MatrCc[i, i,   iPaq, iTest]   <- -1
		}

		if (iPaq <= NBRM) {
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
				
				for (iRow in 1:(NBN-1)) {
					pY[iRow, iMat] = MatR[iMat, col[iRow], iPaq]			
					for (iCol in 1:NBN-1) {
						pX[iRow, iCol, iMat] = MatR[col[iCol], col[iRow], iPaq]
					}
				}
			}
			
			gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,]  <- pX[, ,]
			gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),  ]  <- pY[,  ]
		}
	}			# Essais de bruits (iPaq)
	
	for (iRep in 1:NBR) {
		NRep 	<- Replic[iRep]
		gX1		<- gX
		gY1		<- gY
			
		if (NRep < NBRM) {
			v	<- seq(NRep*(NBN-1)+1, NBRM*(NBN-1), 1)				# Lignes à éliminer
			gX1	<- gX1[-v,,]
			gY1	<- gY1[-v,]
		}
		
		indN	<- 1
		for (iRow in 1:NBN) {
			col <- 1:NBN
			col <- col[-iRow]
			
			Form1 <- paste("gY1[,",iRow,"]~gX1[,,",iRow,"]+0")		# Pour comparaison avec la régression linéaire directe
			# cat(NRep, " REPLICATS iTest ", iTest, " iRow ", iRow, " Form ", Form1, "\n")
				
			pQ <- lm(as.formula(Form1))
			VMean[iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate
			VMin [iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate - 1.96*broom::tidy(pQ)$std.error
			VMax [iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate + 1.96*broom::tidy(pQ)$std.error
			
			for (iCol in col) {	
				ValNom[indN] 	<- 6*(iRow-1)+iCol
				ValVrai[indN]	<- matrTh[iRow, iCol]
				
				#	Réplicats
				ValMean[indN, iRep, iTest] 	<-	VMean[iRow, iCol, iRep, iTest]
				ValMin [indN, iRep, iTest] 	<-	VMin [iRow, iCol, iRep, iTest]
				ValMax [indN, iRep, iTest] 	<-	VMax [iRow, iCol, iRep, iTest]	
				
				#	Calcul des quantiles (bootstrap)
				ValMeanBP[indN, iTest] 			<- mean(MatrCc[iRow,iCol,,iTest])
				qq	<- quantile(MatrCc[iRow,iCol,,iTest], probs=c(0.025,0.975))		# Donne la plage de IC à 95%
				ValMinBP [indN, iTest] 			<- qq[1]
				ValMaxBP [indN, iTest] 			<- qq[2]	
		
				indN = indN +1
			}
		}
	}			# Réplicats
}				# Tests (niveaux de bruit)

ValNom2	<-	c(1:(NBN*NBN))
for(i in NBN:1) {
	ValNom2	<-	ValNom2[-(1 + (i-1)*(NBN+1))]		# Donne n : 2,3,4,5,6, 7,9,10, ... ,35
}
for(i in 1:NBN) {
	ValNom2[((i-1)*(NBN-1)+1) : ((i-1)*(NBN-1)+5)] 	<-	ValNom2[((i-1)*(NBN-1)+1) : ((i-1)*(NBN-1)+5)] - (i-1)*NBN	# Donne j : 2,3,4,5,6, 1,3,4, ... ,5
}
ValNom2[]	<- as.character(ValNom2[])

# k= 0.001

iTest = 1
VMin	<- min(ValMin[,1,iTest], ValMin[,2,iTest], ValMinBP[,iTest])

matr3	<- matrix(VMin, nrow=NBN*NBN, ncol=3)	# 3 réplicats
matr5	<- matrix(VMin, nrow=NBN*NBN, ncol=3)	# 5 réplicats
matrBP	<- matrix(VMin, nrow=NBN*NBN, ncol=3)	# Bootstrap
matrVrai <-matrix(VMin,	nrow=NBN*NBN, ncol=1)	# Valeurs vraies

matNom 		<- vector(length=NBN*NBN)
matNom[]    <- " "
for (i in 1:NBN) {
	matNom[(1+(i-1)*NBN) : (5+(i-1)*NBN)]  		<- ValNom2[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1))]
	
	matr3  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 1]	<- ValMean[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 1, iTest]
	matr3  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 2]	<- ValMin [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 1, iTest]
	matr3  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 3]	<- ValMax [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 1, iTest]
	
	matr5  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 1]	<- ValMean[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 2, iTest]
	matr5  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 2]	<- ValMin [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 2, iTest]
	matr5  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 3]	<- ValMax [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 2, iTest]
	
	matrBP [(1+(i-1)*NBN) : (5+(i-1)*NBN), 1]	<- ValMeanBP[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), iTest]
	matrBP [(1+(i-1)*NBN) : (5+(i-1)*NBN), 2]	<- ValMinBP [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), iTest]
	matrBP [(1+(i-1)*NBN) : (5+(i-1)*NBN), 3]	<- ValMaxBP [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), iTest]
	
	matrVrai [(1+(i-1)*NBN) : (5+(i-1)*NBN), 1]	<- ValVrai[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1))]	
}

matr3.df	<- data.frame(c(1:(NBN*NBN)),      matNom, matr3)
matr5.df	<- data.frame(c(1:(NBN*NBN))+0.25, matNom, matr5)
matrBP.df	<- data.frame(c(1:(NBN*NBN))-0.25, matNom, matrBP)
matrVrai.df	<- data.frame(c(1:(NBN*NBN)),      matNom, matrVrai)

colnames(matr3.df)	<- c("i_j", "nom", "min", "lower", "upper")
colnames(matr5.df)	<- c("i_j", "nom", "min", "lower", "upper")
colnames(matrBP.df)	<- c("i_j", "nom", "min", "lower", "upper")
colnames(matrVrai.df)	<- c("i_j", "nom", "vrai")

Breaks			<-  c(1:(NBN*NBN))
NomFic 			<- 	"MRA3_4c1c.png"
Titre			<- paste("k = ", Noise[iTest])			# Distribution des ri,j
png(filename=NomFic)	# fichier png contenant les IC des "ri,j" (sensibilité au bruit). Perturb +50%	
ggplot() + geom_linerange(data=matr3.df,  aes(x=i_j, ymin=lower, ymax=upper), size=0.9, color="blue")+
		   geom_linerange(data=matr5.df,  aes(x=i_j, ymin=lower, ymax=upper), size=0.9, color="red")+
		   geom_linerange(data=matrBP.df, aes(x=i_j, ymin=lower, ymax=upper), size=0.9, color="black")+
		   geom_hline(yintercept=0, size=1, color="orange")+
		   geom_point(data=matrVrai.df, aes(x=i_j, y=vrai), size=0.7, color="green") +
		   labs(title=Titre, x="", y="CI 95%") + scale_x_continuous(breaks=Breaks, labels=matNom)		# y="IC 95% de ri,j"
dev.off()

# k= 0.005

iTest = 2
VMin	<- min(ValMin[,1,iTest], ValMin[,2,iTest], ValMinBP[,iTest])

matr3	<- matrix(VMin, nrow=NBN*NBN, ncol=3)
matr5	<- matrix(VMin, nrow=NBN*NBN, ncol=3)
matrBP	<- matrix(VMin, nrow=NBN*NBN, ncol=3)
matNom 		<- vector(length=NBN*NBN)
matNom[]    <- " "
for (i in 1:NBN) {
	matNom[(1+(i-1)*NBN) : (5+(i-1)*NBN)]  		<- ValNom2[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1))]
	
	matr3  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 1]	<- ValMean[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 1, iTest]
	matr3  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 2]	<- ValMin [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 1, iTest]
	matr3  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 3]	<- ValMax [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 1, iTest]
	
	matr5  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 1]	<- ValMean[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 2, iTest]
	matr5  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 2]	<- ValMin [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 2, iTest]
	matr5  [(1+(i-1)*NBN) : (5+(i-1)*NBN), 3]	<- ValMax [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), 2, iTest]
	
	matrBP [(1+(i-1)*NBN) : (5+(i-1)*NBN), 1]	<- ValMeanBP[(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), iTest]
	matrBP [(1+(i-1)*NBN) : (5+(i-1)*NBN), 2]	<- ValMinBP [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), iTest]
	matrBP [(1+(i-1)*NBN) : (5+(i-1)*NBN), 3]	<- ValMaxBP [(1+(i-1)*(NBN-1)) : (5+(i-1)*(NBN-1)), iTest]	
}
matr3.df	<- data.frame(c(1:(NBN*NBN)),      matNom, matr3)
matr5.df	<- data.frame(c(1:(NBN*NBN))+0.25, matNom, matr5)
matrBP.df	<- data.frame(c(1:(NBN*NBN))-0.25, matNom, matrBP)
colnames(matr3.df)	<- c("i_j", "nom", "min", "lower", "upper")
colnames(matr5.df)	<- c("i_j", "nom", "min", "lower", "upper")
colnames(matrBP.df)	<- c("i_j", "nom", "min", "lower", "upper")

Breaks			<-  c(1:(NBN*NBN))
NomFic 			<- 	"MRA3_4c2c.png"
Titre			<- paste("k = ", Noise[iTest])		# Distribution des ri,j
png(filename=NomFic)	# fichier png contenant les IC des "ri,j" (sensibilité au bruit). Perturb +50%	
ggplot() + geom_linerange(data=matr3.df,  aes(x=i_j, ymin=lower, ymax=upper), size=0.9, color="blue")+
		   geom_linerange(data=matr5.df,  aes(x=i_j, ymin=lower, ymax=upper), size=0.9, color="red")+
		   geom_linerange(data=matrBP.df, aes(x=i_j, ymin=lower, ymax=upper), size=0.9, color="black")+
		   geom_hline(yintercept=0, size=1, color="orange")+
		   geom_point(data=matrVrai.df, aes(x=i_j, y=vrai), size=0.7, color="green") +		   
		   labs(title=Titre, x="", y="CI 95%") + scale_x_continuous(breaks=Breaks, labels=matNom)		# y="IC 95% de ri,j"
dev.off()

# Dessins des deux figures de gauche


NBN 	<- 6										# Nombre de noeuds
Verbose	<- FALSE									# Si TRUE, afficher certains commentaires

Noise		<-	c(0.001, 0.005, 0.01)				# 3 niveaux de bruit
Replic		<-  c(3, 5)								# 3 et 5 réplicats techniques
NPerturb 	<-  1.5									# Perturbation +50%
Methods		<-  c("Th", "BS IC95", "PVal 5", "Lasso", "RLD", "Step F", "Step B", "Step BO", "Step F_0", "Step B_0", "Step BO_0")		
																						# Méthodes testées
Steps		<- 	c("Step F", "Step B", "Step BO", "Step F_0", "Step B_0", "Step BO_0")	# Méthodes "Step" (forward, backward, both, avec ou sans intercept forcé à 0)																		
ScoreTest	<-	c("T+", "T-", "T0", "F+", "F-", "F0", "Sensib", "Specif")				# Scores mesurés

NBT 	<- length(Noise)
NBK		<- 100										# 100 tirages de bruit
NBR		<- length(Replic)							# Les essais de réplicats
NBRM	<- max(Replic)								# Nb max de réplications à effectuer
NBM		<- length(Methods)							# Nb de méthodes à tester
NBRST	<- length(Steps)							# Nb de méthodes STEP
NBSC	<- length(ScoreTest)						# Nb. de scores testés	
		
Form	<- vector(length=NBN)						# Formule décrivant la régression pour le calcul de MatrCc
XMesNP	<- array(dim=c(NBN, NBK, NBT))				# Valeur mesurées, système non perturbé
XMesP	<- array(dim=c(NBN, NBN, NBK, NBT))			# Valeur mesurées, suite à une perturbation
MatR	<- array(dim=c(NBN, NBN, NBK))				# Matrice R
MatrCc	<- array(dim=c(NBN, NBN, NBK, NBT))			# Matrice r calculée pour chaque valeur de bruit (bootstrap)
		
pX 	<- array(dim=c(NBN-1, NBN-1, NBN))				# ligne (iRow), colonne (iCol), n° matrice (iMat)
pY 	<- array(dim=c(NBN-1, NBN))						# ligne, n° matrice
gX 	<- array(dim=c(NBRM*(NBN-1), NBN-1, NBN))		# ligne (iRow), colonne (iCol), n° matrice (iMat)
gY 	<- array(dim=c(NBRM*(NBN-1), NBN))				# ligne, 1, n° matrice
gX1	<- array(dim=c(NBRM*(NBN-1), NBN-1, NBN))		# Copie de gX, limitée au nombre de tirages correspondant au nb. de réplicats
gY1	<- array(dim=c(NBRM*(NBN-1), NBN))				# Copie de gY, idem

rL 		<- matrix(nrow=NBN, ncol=1)					# Résultat de la méthode de Lasso, NBN-1 * NBN-1 (ie sol. de Yi = Ai * Xi)
rLasso  <- array(dim=c(NBN, NBN, NBR, NBT))			# Intégration des résultats Xi, ajout des -1 sur la diagonale et transpo ==> donne la matrice "r"

MatrStp <- array(dim=c(NBN, NBN, NBR, NBT, NBRST))	# Calcul de MatrCc par les méthodes STEP
rij 	<- matrix(nrow=1, ncol=NBN-1)				# Calcul intermédiaire des rij
		
VMean	<- array(dim=c(NBN, NBN, NBR, NBT))			# Moyenne calculée
VMin	<- array(dim=c(NBN, NBN, NBR, NBT))			# Valeur mini. de l'intervalle de confiance à 95%
VMax	<- array(dim=c(NBN, NBN, NBR, NBT))			# Valeur maxi. de l'intervalle de confiance à 95%
VCc		<- array(dim=c(NBN, NBN, NBR, NBT))			# MatrCc calculée avec les réplicats (iRep)

ValNom		<- vector(length=(NBN-1)*NBN)			# Nom (n) correspondant à (i,j)
ValMean		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Moyenne des réplicats
ValMin		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Min IC95 des réplicats (correspond à p_value à 5%)
ValMax		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Max IC95 des réplicats (correspond à p_value à 5%)
ValCc		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# MatrCc calculée avec les réplicats (iRep) (indices homogènes à ValMean etc ...)
ValStp		<- array(dim=c((NBN-1)*NBN, NBR, NBT, NBRST))	# MatrStp calculée avec les réplicats (iRep) (indices homogènes à ValMean etc ...)
ValMeanBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Moyenne des 100 valeurs (cf Boxplot)
ValMinBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Quantile à 2,5%  des 100 valeurs
ValMaxBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Quantile à 97,5% des 100 valeurs

Result 	<- array(dim=c((NBN-1)*NBN, NBR, NBT, NBM))	# Résultats obtenus
Score 	<- array(dim=c(NBSC, NBR, NBT, NBM-1))		# Vrais +, vrais - vrais 0, faux +, faux - faux 0, Sensib, Specif    --- par comparaison à la matrice "r" théorique
											
dimnames(Result)	<- list(NULL, NULL, NULL, Methods)				# Nom des méthodes correspondant aux résultats
dimnames(Score)		<- list(ScoreTest, NULL, NULL, Methods[2:NBM])	# Nom des scores testés
									
# Recherche de la formule décrivant la régression pour le calcul de MatrCc (Intercept nul)

for (i in 1:NBN) {
	col <- 1:NBN
	col <- col[-i]
	
	Form[i] <- paste("MatR[", i, ",col,iPaq]-", sep="")			# Formule décrivant la régression
	for (j in col) {
		Form[i] <- paste(Form[i], "MatR[", sep="A")
		Form[i] <- paste(Form[i], j, sep="")
		Form[i] <- paste(Form[i], ",col", sep="")
		Form[i] <- paste(Form[i], ",iPaq]", sep= "")
	}
	
	Form[i] <- paste(Form[i], "+0", sep="")
	Form[i] <- str_replace_all(Form[i], pattern = "-A", replacement = "~")
	Form[i] <- str_replace_all(Form[i], pattern = "A", replacement = "+")
	if (Verbose) {
		cat("BOOTSTRAP i ", i, " Form ", Form[i], "\n")
	}
}

set.seed(12345)								# Pour pouvoir régénérer les mêmes séquences

for (iTest in c(1:NBT)) {
	# Pour chaque niveau de bruit (iTest), on effectue 100 tirages (NBK) et on calcule 100 matrices "r" (chaque matrice correspond à un tirage)
	# Puis on effectue des statistiques sur les erreurs obtenues : moyenne et sd des ri,j, ainsi que l'erreur moyenne et son écart-type
	# La régression linéaire porte sur une valeur seule. Cette méthode (6 matrices 5*5) est équivalente à l'inversion de la matrice 6*6
	
	sd <- Noise[iTest]*XMesNPMoy			# écart-type du bruit

	for (iPaq in c(1:NBK)) {
		#	Entrée des "données de mesure" bruitées
		ss <- multiroot (f=F, start= st)	# Mesure non perturbée  - utilise le package rootSolve
		XMesNP[,iPaq,iTest]	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		
		V4 <- NPerturb*V4		# P1
		ss <- multiroot (f=F, start= st)
		XMesP[,1,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V4 <- V4/NPerturb		# Retour à la valeur initiale
		
		V3 <- NPerturb*V3		# P2
		ss <- multiroot (f=F, start= st)
		XMesP[,2,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V3 <- V3/NPerturb		# Retour à la valeur initiale
		
		V8 <- NPerturb*V8		# P3
		ss <- multiroot (f=F, start= st)
		XMesP[,3,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V8 <- V8/NPerturb		# Retour à la valeur initiale
		
		V7 <- NPerturb*V7		# P4
		ss <- multiroot (f=F, start= st)
		XMesP[,4,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V7 <- V7/NPerturb		# Retour à la valeur initiale
		
		V12 <- NPerturb*V12		# P5
		ss <- multiroot (f=F, start= st)
		XMesP[,5,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V12 <- V12/NPerturb		# Retour à la valeur initiale
		
		V11 <- NPerturb*V11		# P6
		ss <- multiroot (f=F, start= st)
		XMesP[,6,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V11 <- V11/NPerturb		# Retour à la valeur initiale	
	
		for (iMat in c(1:NBN)) {
			MatR[,iMat,iPaq]  <- 2 * (XMesP[,iMat,iPaq,iTest]-XMesNP[,iPaq,iTest]) / (XMesP[,iMat,iPaq,iTest]+XMesNP[,iPaq,iTest])
		}
		
		for (i in 1:NBN) {
			col <- 1:NBN
			col <- col[-i]
					
			MatrCc[i, col, iPaq, iTest]   <- (lm(as.formula(Form[i])))$coefficients
			MatrCc[i, i,   iPaq, iTest]   <- -1
		}

		if (iPaq <= NBRM) {
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
				
				for (iRow in 1:(NBN-1)) {
					pY[iRow, iMat] = MatR[iMat, col[iRow], iPaq]			
					for (iCol in 1:NBN-1) {
						pX[iRow, iCol, iMat] = MatR[col[iCol], col[iRow], iPaq]
					}
				}
			}
			
			gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,]  <- pX[, ,]
			gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),  ]  <- pY[,  ]
		}
	}			# Essais de bruits (iPaq)
	
	#	Application des différentes méthodes
	
	for (iRep in 1:NBR) {							# 3 ou 5 réplicats
		NRep 	<- Replic[iRep]						# Nb de réplicats
		gX1		<- gX
		gY1		<- gY
			
		if (NRep < NBRM) {
			v	<- seq(NRep*(NBN-1)+1, NBRM*(NBN-1), 1)				# Lignes à éliminer (cas où il y a moins de réplicats que NBRM)
			gX1	<- gX1[-v,,]
			gY1	<- gY1[-v,]
		}
		
		MatrStp[,, iRep, iTest,] <- 0
				
		indN	<- 1
		for (iRow in 1:NBN) {
			col <- 1:NBN
			col <- col[-iRow]
			
			#	Méthode BS IC95  (ie. MRA de base) et p-Value (5%)
			Form1 <- paste("gY1[,",iRow,"]~gX1[,,",iRow,"]+0")		# Pour comparaison avec la régression linéaire directe
													# Intercept forcé à 0
													# On les valeurs correspondant à Nrep (3 ou 5) tirages
			if (Verbose) {
				cat(NRep, " REPLICATS iTest ", iTest, " iRow ", iRow, " Form ", Form1, "\n")
			}

			pQ <- lm(as.formula(Form1))
			VMean[iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate
			VMin [iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate - 1.96*broom::tidy(pQ)$std.error
			VMax [iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate + 1.96*broom::tidy(pQ)$std.error
			VCc	 [iRow, col, iRep, iTest]  <- lm(pQ)$coefficients
																# Intervalle de confiance à 95%
			
			#	Méthode de LASSO (avec choix automatique de Lambda)
			cv_model 	<- cv.glmnet(gX1[,,iRow], gY1[,iRow], alpha = 1)	# Fit lasso regression model using k-fold cross-validation
			best_lambda <- cv_model$lambda.min
			best_model 	<- glmnet(gX1[,,iRow], gY1[,iRow], alpha = 1, lambda = best_lambda)	# View coefficients of best model
			rL			<-	coef(best_model)
			rLasso[iRow, iRow, iRep, iTest] <- -1						
			rLasso[col,iRow, iRep, iTest]	<- rL[2:NBN,1]
			rLasso[,, iRep, iTest] <- t(rLasso[,,iRep, iTest])	# on trouve ici l'équivalent de la matrice "r"				
			
			#	Méthodes STEP
			Donnees  <- data.frame(gY1[,iRow], gX1[,,iRow])
			Donnees	 <- rename(Donnees, "Y" ="gY1...iRow.")		# Le suivi est plus clair comme cela
			cc 		 <- colnames(Donnees)						# Nom des colonnes de "Données". La première est "Y"
			cc		 <- cc[-1]									# Il reste le nom des coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
			colnames(rij) <- cc									# step indique le nom des colonnes correspondant aux coefficients conservés			
	
			for (iStep in 1:NBRST) {							# méthodes "step" (fonction step() -- package "stats")			
				switch(iStep, 
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},		# Forward
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)},		# Both
					   {intercept_only <- lm(Y ~ 1+0, data=Donnees)
						all 	<- lm(Y ~ .+0, data=Donnees)
						forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},		# Forward Intercept=0
					   {intercept_only <- lm(Y ~ 1+0, data=Donnees)
						all 	<- lm(Y ~ .+0, data=Donnees)
						forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward Intercept=0
					   {intercept_only <- lm(Y ~ 1+0, data=Donnees)
						all 	<- lm(Y ~ .+0, data=Donnees)
						forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)})		# Both Intercept=0
					
				ll  	<- length(forward$coefficients)				# Nombre de coefficients conservés par step
				nn		<- names(forward$coefficients)				# Nom de ces coefficients
				if (Verbose) {
					cat ("STEP ", iStep, " iRow ", iRow, " Coefs. retenus ", nn, " Valeurs ", forward$coefficients[nn[i]], "\n")
				}
	
				rij[1,]	<- 0
				if (ll >= 2) {
					for (i in 2:ll) {								# nn[1] = "(Intercept)"
						rij[1, nn[i]] <- forward$coefficients[nn[i]]
					}
				}
	
				MatrStp[iRow, col,  iRep, iTest, iStep]		<- rij[1,]
				MatrStp[iRow, iRow, iRep, iTest, iStep]  	<- -1
			}
			
			for (iCol in col) {	
				ValNom[indN] 	<- 6*(iRow-1)+iCol
				
				#	Valeurs théoriques
				if (matrTh[iRow, iCol] > 0) {
					Result [indN, iRep, iTest, 1] = 1			# 1° méthode : valeurs théoriques
				}	else if (matrTh[iRow, iCol] == 0) {
					Result [indN, iRep, iTest, 1] = 0	
				}	else {
					Result [indN, iRep, iTest, 1] = -1
				}									
				
				#	Réplicats
				ValMean[indN, iRep, iTest] 	<-	VMean[iRow, iCol, iRep, iTest]
				ValMin [indN, iRep, iTest] 	<-	VMin [iRow, iCol, iRep, iTest]
				ValMax [indN, iRep, iTest] 	<-	VMax [iRow, iCol, iRep, iTest]
				ValCc  [indN, iRep, iTest] 	<-	VCc  [iRow, iCol, iRep, iTest]
				for (iStep in 1:NBRST) {
					ValStp [indN, iRep, iTest, iStep] 	<-	MatrStp [iRow, iCol, iRep, iTest, iStep]
				}
				
				#	Calcul des quantiles (bootstrap)
				ValMeanBP[indN, iTest] 			<- mean(MatrCc[iRow,iCol,,iTest])
				qq	<- quantile(MatrCc[iRow,iCol,,iTest], probs=c(0.025,0.975))		# Donne la plage de IC à 95%
				ValMinBP [indN, iTest] 			<- qq[1]
				ValMaxBP [indN, iTest] 			<- qq[2]					
				
				if (ValMinBP[indN, iTest] > 0) {
					Result [indN, iRep, iTest, 2] = 1			# 2° méthode : bootstrap sur les valeurs min et max de IC95
				}	else if (ValMaxBP[indN, iTest] < 0) {
					Result [indN, iRep, iTest, 2] = -1	
				}	else {
					Result [indN, iRep, iTest, 2] = 0
				}

				#	p-Value à 5%
				if (ValMin[indN, iRep, iTest] > 0) {
					Result [indN, iRep, iTest, 3] = 1			# 3° méthode : p-Value à 5%
				}	else if (ValMax[indN, iRep, iTest] < 0) {
					Result [indN, iRep, iTest, 3] = -1	
				}	else {
					Result [indN, iRep, iTest, 3] = 0
				}
				
				#	Méthode de Lasso
				if (rLasso[iRow, iCol, iRep, iTest] > 0) {
					Result [indN, iRep, iTest, 4] = 1			# 4° méthode : méthode de Lasso
				}	else if (rLasso[iRow, iCol, iRep, iTest] == 0) {
					Result [indN, iRep, iTest, 4] = 0	
				}	else {
					Result [indN, iRep, iTest, 4] = -1
				}
				
				indN <- indN +1
			}	# iCol
		}		# iRow
		
		ThRLD	<- 0.25 * max(ValCc[, iRep, iTest])
		for (ind in 1:(indN-1)) {
			if (ValCc[ind, iRep, iTest] > ThRLD) {
				Result [ind, iRep, iTest, 5] = 1				# 5° méthode : RLD à seuil
			}	else if (ValCc[ind, iRep, iTest] < -ThRLD) {
				Result [ind, iRep, iTest, 5] = -1	
			}	else {
				Result [ind, iRep, iTest, 5] = 0
			}
		}		# ind				
			
		for (iStep in 1:NBRST) {
			for (ind in 1:(indN-1)) {
				if (ValStp[ind, iRep, iTest, iStep] > 0) {
					Result [ind, iRep, iTest, 5+iStep] = 1		# 6° à 11° méthodes : STEP (for, back, both, for0, back0, both0)
				}	else if (ValStp[ind, iRep, iTest, iStep] < 0) {
					Result [ind, iRep, iTest, 5+iStep] = -1	
				}	else {
					Result [ind, iRep, iTest, 5+iStep] = 0
				}
			}		# ind
		}		# iStep
	}			# Réplicats
}				# Tests (niveaux de bruit)

#	Scores obtenus

for (iTest in 1:NBT)		{
	for (iRep in 1:NBR)		{
		for (iMet in 2:NBM)	{						# iMet = 1 correspond au résultat exact (Th éorique)
			Score[,iRep, iTest, iMet-1]	<- 0		# mise à 0 des compteurs
			
			for (ind in 1:(indN-1)) {				# ind balaye tous les couples (i,j), i # j
				if (Result[ind, iRep, iTest, iMet]  > 0) {
					if (Result[ind, iRep, iTest, 1] > 0) {		# le dernier indice (1) fait référence à la valeur exacte
						Score["T+", iRep, iTest, iMet-1] <- Score["T+", iRep, iTest, iMet-1] +1
					} else {
						Score["F+", iRep, iTest, iMet-1] <- Score["F+", iRep, iTest, iMet-1] +1
					}
				}
				if (Result[ind, iRep, iTest, iMet]  < 0) {
					if (Result[ind, iRep, iTest, 1] < 0) {
						Score["T-", iRep, iTest, iMet-1] <- Score["T-", iRep, iTest, iMet-1] +1
					} else {
						Score["F-", iRep, iTest, iMet-1] <- Score["F-", iRep, iTest, iMet-1] +1
					}
				}
				if (Result[ind, iRep, iTest, iMet]  == 0) {
					if (Result[ind, iRep, iTest, 1] == 0) {
						Score["T0", iRep, iTest, iMet-1] <- Score["T0", iRep, iTest, iMet-1] +1
					} else {
						Score["F0", iRep, iTest, iMet-1] <- Score["F0", iRep, iTest, iMet-1] +1
					}
				}				
			}	# ind   (Test des valeurs calculées pour tous les couples i,j)
			
			P	<- sum (abs(Result[, iRep, iTest, 1]))		# Nb. d'elts non nuls (+ et -)
			N 	<- NBN*(NBN-1) - P							# Nb. d'elts nuls
			Score["Sensib", iRep, iTest, iMet-1] <- (Score["T+", iRep, iTest, iMet-1]+Score["T-", iRep, iTest, iMet-1]) / P
			Score["Specif", iRep, iTest, iMet-1] <- (Score["T0", iRep, iTest, iMet-1]) / N
		}		# iMet	(Les différentes méthodes)
	}			# iRep	(Nb. de réplicats techniques)
}				# iTest (Niveaux de bruit)

#	Dessin du tableau de droite





#
#	Figure 4 : Dream Challenge 4
#
#	Pour obtenir les fichiers de données, aller sur "https://www.synapse.org/#!Synapse:syn3049712/wiki/74630"
#	et télécharger :
##	DREAM4_InSilico_Size10.zip
##	DREAM4_InSilico_Size100.zip
##	Télécharger aussi 
##	DREAM4_InSilicoNetworks_GoldStandard.zip, à télécharger par "https://www.synapse.org/#!Synapse:syn3049736"
##	(Ouvrir "DREAM4_Challenge2_GoldStandards", puis "Size10", " Size 10 bonus round" et enfin " "DREAM4_GoldStandard_InSilico_Size10_1.tsv" 
##	etc… pour les autres fichiers).

#
Verbose	<- FALSE
NBFIC	<- 5				# Il y a 5 fichiers pour une taille donnée
indNBN	<- 1				# Valeurs possibles : 1 (pour 10 noeuds), 2 (pour 100 noeuds)

Methods		<-  c("MRA-KD", "MRA-KO", "PVAL", "LASSO", "RLD", "LARS", "STEP-Fo", "STEP-Ba", "STEP-Bo")
											# Méthodes testées pour l'article BioInformatics  + "ARACNE", "CLR", "LARS"
ScoreTest	<-	c("T+-", "T0", "F+-", "F0", "Sensib", "Specif")	# Scores mesurés
ScoreMoyR	<-	c("mean", "sd")									# Scores globaux pour les 5 fichiers d'une taille
ScoreMoyC	<-	c("Sensib", "Specif")							# Valeurs à afficher

NBM		<- length(Methods)										# Nb de méthodes à tester
NBSC	<- length(ScoreTest)									# Nb. de scores testés

vNBN	<- c(10, 100)																			# Nb noeuds
vNBK	<- c(2, 2)																				# Nb de types de perturbations (KD ou KO)
vTOP	<- c(0.2, 0.25)																			# "top" des arêtes à prendre

vRoot	<- c("Documents/DreamChallenge10/", "Documents/DreamChallenge100/")						# Nom de la racine des fichiers
vRSMTxt	<- c("insilico_size10_",  "insilico_size100_")											# Texte à afficher dans les tables de regroupement
vRSMFic	<- c("MRA3_Resume_ZF10",  "MRA3_Resume_ZF100")											# Nom des fichiers de regroupement
vSolTxt	<- c("DREAM4_GoldStandard_InSilico_Size10_", "DREAM4_GoldStandard_InSilico_Size100_")	# Partie du nom des fichiers solution
ParNds	<- data.frame(vNBN, vNBK, vTOP, vRoot, vRSMTxt, vRSMFic, vSolTxt)						# Tableau des données précédentes, correspondant aux nb. de noeuds possibles

NBN = ParNds$vNBN[indNBN]
NBK = ParNds$vNBK[indNBN]					# NBK ne correspond plus à du bruit ajouté, mais à des types de perturbations (KD ou KO)

XMesNP		<- array(dim=c(NBN, NBK))		# Valeur mesurées, système non perturbé
XMesP		<- array(dim=c(NBN, NBN, NBK))	# Valeur mesurées, suite à une perturbation
Solution	<- array(dim=c(NBN, NBN))		# Solution : vraie structure du réseau
Solution2	<- array(dim=c(NBN, NBN))		# Idem "Solution", mais les n° de ligne et colonne correspondent aux n° de "levels" ("G1", "G10", "G2", "G3" etc)
Solution3	<- array(dim=c(NBN, NBN))		# Intermédiaire de calcul

MatR		<- array(dim=c(NBN, NBN, NBK))	# Matrice R
rij 		<- matrix(nrow=1, ncol=NBN-1)	# Calcul intermédiaire des rij
rL 			<- vector(length=NBN)			# Résultat de la méthode de Lasso, NBN-1 * NBN-1 (ie sol. de Yi = Ai * Xi)
MatrCc		<- array(dim=c(NBN, NBN))		# Matrice r calculée #  sauf par méthode PVAL    # selon la méthode MRA, RLD ou STEP ou LASSO
vrCc  		<- vector(length = NBN*NBN)		# MatrCc en vecteur
MatrCc1		<- array(dim=c(NBN, NBN))		# Matrice r résultat de l'algorithme
											# Elle est calculée en tronquant les éléments de MatrCc à 0 et 1 (0 sur la diagonale, comme Solution)

MatPVal		<- array(dim=c(NBN, NBN))		# Calcul direct de la p Value

VMean		<- array(dim=c(NBN, NBN))		# Valeur moyenne
VMin		<- array(dim=c(NBN, NBN))		# Valeur mini. de l'intervalle de confiance à 95%
VMax		<- array(dim=c(NBN, NBN))		# Valeur maxi. de l'intervalle de confiance à 95%

pX   <- array(dim=c(NBN-1, NBN-1, NBN))						# ligne (iRow), colonne (iCol), n° matrice (iMat)
pY   <- array(dim=c(NBN-1, NBN))							# ligne, n° matrice
gX   <- array(dim=c(NBK*(NBN-1), NBN-1, NBN))				# ligne (iRow), colonne (iCol), n° matrice (iMat)
gY   <- array(dim=c(NBK*(NBN-1), NBN))						# ligne, 1, n° matrice
gX1  <- array(dim=c(NBN-1, NBN-1, NBN))						# matrice gX tronquée pour ne prendre que les valeurs KD ou KO
gY1  <- array(dim=c(NBN-1, NBN))							# idem gX1

Score 		<- array(0, dim=c(NBSC, NBM, NBFIC))			# Vrais +- (TP), vrais 0 (TN), faux +- (FP), faux 0 (FN)
dimnames(Score)		<- list(ScoreTest, Methods, NULL)
ScoreMoy	<- array(0, dim=c(length(ScoreMoyR), length(ScoreMoyC), NBM))	# Moyenne des scores d'une série de fichiers
dimnames(ScoreMoy)	<- list(ScoreMoyR, ScoreMoyC, Methods)

NomFic <- paste(ParNds$vRoot[indNBN], "Classement_DC4.txt", sep="")
for(indFIC in 1:NBFIC) {
	valw <- paste("indFic", indFIC, sep=";")
	
	#	1/ Lecture des données
		
	val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_wildtype.tsv", sep="")
	XMesNPLu 	<- fread(val, data.table=F)						# "Documents/DreamChallenge10/insilico_size10_1_wildtype.tsv"
	val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_knockdowns.tsv", sep="")
	XMesPLu1 	<- fread(val, data.table=F) 					# ("Documents/DreamChallenge10/insilico_size10_1_knockdowns.tsv"
	val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_knockouts.tsv", sep="")
	XMesPLu2 	<- fread(val, data.table=F) 					# "Documents/DreamChallenge10/insilico_size10_1_knockouts.tsv"
	
	val			<- paste(ParNds$vRoot[indNBN], ParNds$vSolTxt[indNBN], indFIC, ".tsv", sep="")
	SolutLu  	<- fread(val, data.table=F) 					# "Documents/DreamChallenge10/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
	SolutLu$V1	<- as.factor(SolutLu$V1)
	SolutLu$V2	<- as.factor(SolutLu$V2)
	
	for (i in 1:NBN) {
		XMesNP[i,1] = XMesNPLu[1,i]
		XMesNP[i,2] = XMesNPLu[1,i]
		
		for (j in 1:NBN) {
			XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation à -50%	
												# ATTENTION les données de Dream Cha. sont transposées par rapport à mon écriture
			XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation à -100%
		}
	}
	
	Solution2[,] = 0
	for (i in 1:(NBN*(NBN-1))) {
	#	cat("i ", i, " V3 ", SolutLu[i,]$V3, " V1 ", SolutLu[i,]$V1, " : ", as.numeric(SolutLu[i,]$V1), " V2 ", SolutLu[i,]$V2, " : ", as.numeric(SolutLu[i,]$V2), "\n")
		if(SolutLu[i,]$V3 == 1) {
			Solution2[as.numeric(SolutLu[i,]$V1), as.numeric(SolutLu[i,]$V2)] = 1
		} else {
			break								# On suppose que tous les 1 sont au début du fichier
		}
	}
	
	numero <- as.numeric(sort(as.character(c(1:NBN))))
	
	for (i in 1:NBN) {
		Solution3[,numero[i]] = Solution2[,i]
	}
	for (i in 1:NBN) {
		Solution[numero[i],] = Solution3[i,]	# On retrouve la matrice de connexion classique (noeuds "G1", "G2", ... , "G10")
	}
	Solution <- t(Solution)						# Deuxième configuration 
												# correspond à la 1° interprétation en réalité : erreur sur l'écriture de la matrice
	
	for (iPaq in c(1:NBK)) {					# On prend toutes les données disponibles : Knock Down et Knock Out
		#	Les "données de mesure", fournies par Dream Challenge, ont été chargées dans XMesNP et XMesP précédemment	
		for (iPert in c(1:NBN)) {
			MatR[,iPert,iPaq]  <- 2 * (XMesP[,iPert,iPaq]-XMesNP[,iPaq]) / (XMesP[,iPert,iPaq]+XMesNP[,iPaq])
		}
		
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
			
			for (iRow in 1:(NBN-1)) {
				pY[iRow, iMat] = MatR[iMat, col[iRow], iPaq]
				
				for (iCol in 1:NBN-1) {
					pX[iRow, iCol, iMat] = MatR[col[iCol], col[iRow], iPaq]
				}
			}
		}
		
		gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,]  <- pX[, ,]
		gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ]   <- pY[, ]
	}
	
	for	(iMeth in 1:NBM) {
		Method  <- Methods[iMeth]		
	
		#	Recherche de MatrCc
		
		if (Method %in% c("MRA-KD", "MRA-KO")) {
		#	Méthode MRA classique
			if (Method == "MRA-KD") {
				rSupp <- seq(NBN, NBK*(NBN-1), 1)	# On supprime les lignes 10 à 18 (resp. 100 à 198) pour KD
			} else {
				rSupp <- seq(1, NBN-1, 1)			# On supprime les lignes 1  à 9  (resp. 1 à 99) pour KO
			}
			
			gX1	<- gX[-rSupp,,]
			gY1	<- gY[-rSupp,]
		
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
					
				Form <- paste("gY1[,",iMat,"]~gX1[,,",iMat,"]+0")
				
				MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
				MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
			}
			MatrCc[is.na(MatrCc)]  <- 0
		}	# Méthodes MRA
		
		if (Method %in% c("RLD", "RLD OPTIM", "RLD OPTIM2")) {
		#	Méthode RLD classique (avec seuils)
		
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
					
				Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")
				
				MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
				MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
			}
			MatrCc[is.na(MatrCc)]  <- 0
		}	# Méthode RLD
		
		if (Method == "PVAL") {						#	p-Value à 5%
			SeuilPV	  <- 0.05
			MatrCc[,] <- 0
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
					
				Form1 <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")	# Pour comparaison avec la régression linéaire directe		
				pQ 	  <- lm(as.formula(Form1))
				
				MatPVal[iMat,col] 	<- broom::tidy(pQ)$p.value
				MatrCc [iMat,which(MatPVal[iMat,] < SeuilPV)] = 1
			}		# boucle iMat
		}	# Méthode PVAL
		
		if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
		#	Méthode "step" (suppression automatique de coefficients par optimisation de AIC : Akaike's Information Criterion)
		
			MatrCc[,] <- 0
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
		
				Donnees  <- data.frame(gY[,iMat], gX[,,iMat])
				#	cat("iMat ", iMat, " Donnees ", colnames(Donnees), "\n")
				Donnees	 <- rename(Donnees, "Y" ="gY...iMat.")		# Le suivi est plus clair comme cela
				cc 		 <- colnames(Donnees)						# Nom des colonnes de "Données". La première est "Y"
				cc		 <- cc[-1]									# Il reste le nom des coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
				colnames(rij) <- cc									# step indique le nom des colonnes correspondant aux coefficients conservés
						
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
				
				ll  	<- length(forward$coefficients)				# Nombre de coefficients conservés par step
				nn		<- names(forward$coefficients)				# Nom de ces coefficients
		
				rij[1,]	<- 0
				if (ll >= 2) {
					for (i in 2:ll) {								# nn[1] = "(Intercept)"
						rij[1, nn[i]] <- forward$coefficients[nn[i]]
					}
				}
				if (Verbose) {
					cat (" iMat ", iMat, " ll ", ll, " nn ", nn, " rij ", rij, "\n")
				}
				MatrCc[iMat, col]	<- rij[1,]
				MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
			}		# boucle iMat
		}	# Méthodes STEP
			
		if (Method == "LASSO") {
		#	Méthode "Lasso", choix automatique de Lambda
		
			MatrCc[,] <- 0		# Initialisation
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
		
				cv_model 	<- cv.glmnet(gX[,,iMat], gY[,iMat], alpha = 1)						# Fit lasso regression model using k-fold cross-validation
				best_lambda <- cv_model$lambda.min	
				best_model 	<- glmnet(gX[,,iMat], gY[,iMat], alpha = 1, lambda = best_lambda)	# View coefficients of best model
				rL			<- coef(best_model)								
		
				MatrCc[iMat, col]	<- rL[2:NBN]
			}		
		}	# Méthode  LASSO
		
		if (Method == "LARS") {
			MatrCc[,] <- 0		# Initialisation
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
		
				pp=lars(gX[,,iMat], gY[,iMat], type="lar")
				rL = coef(pp)
		
				MatrCc[iMat, col]	<- rL[2:NBN]
			}		
		}	# Méthode  LARS		
		
		#	Recherche de la solution
		#	Numérisation du résultat		
		
		MatrMax = max(abs(MatrCc[,]))
		if (Method == "RLD") {
			Th = 0.25*MatrMax
		}	else if (Method %in% c("MRA-KD", "MRA-KO")) {	
			NB1 = trunc(ParNds$vTOP[indNBN]*NBN*NBN)
			vrCc[]  <- MatrCc[]
			vrCc1	<- sort(abs(vrCc), decreasing=TRUE)
			Th	= vrCc1[min(NB1+1, NBN*NBN)]
		}	else {
			Th = 0
		}
		
		for (iRow in 1:NBN) {
			for (iCol in 1:NBN) {
				if (abs(MatrCc[iRow, iCol]) > Th) {
					MatrCc1[iRow, iCol] = 1
				} else {
					MatrCc1[iRow, iCol] = 0
				}
			}
			MatrCc1[iRow, iRow] = 0
		}
		
		#	Test de la solution
		
		Score[, iMeth, indFIC] = 0
		for (iRow in 1:NBN) {
			for (iCol in 1:NBN) {
				if (MatrCc1[iRow, iCol] != 0) {
					if (Solution[iRow, iCol] != 0) {
						Score["T+-", iMeth, indFIC] <- Score["T+-", iMeth, indFIC] +1
					} else {
						Score["F+-", iMeth, indFIC] <- Score["F+-", iMeth, indFIC] +1
					}
				} else {
					if (Solution[iRow, iCol] == 0) {
						Score["T0", iMeth, indFIC] <- Score["T0", iMeth, indFIC] +1
					} else {
						Score["F0", iMeth, indFIC] <- Score["F0", iMeth, indFIC] +1
					}
				}		
			}
		}
		
		P = sum(Solution)			# Nb. vrais positifs
		N = NBN*(NBN-1) - P			# Nb. vrais négatifs
		
		Score["Sensib", iMeth, indFIC] = Score["T+-", iMeth, indFIC] / P
									# (Score["T+-", iMeth] + Score["F0", iMeth])
		Score["Specif", iMeth, indFIC] = (Score["T0", iMeth, indFIC] - NBN) / N
									# (Score["T0", iMeth] + Score["F+-", iMeth] - NBN)
		
		Err	<- (Score["F+-",iMeth,indFIC] + Score["F0",iMeth,indFIC])/(NBN*(NBN-1))
	
		cat(as.character(ParNds$vRSMTxt[indNBN]), indFIC, sep="")
		cat(" Méthode ", Method, " Th ", Th, " Score ", Score[, iMeth, indFIC], " Err ", Err, "\n")
	
		valw <- paste(valw, Method, "Th", Th, "Score", paste(Score[,iMeth,indFIC], collapse=";"), "Err", Err,  sep=";")
	}	# boucle sur les méthodes
	
	cat("\n\n")
	write(valw, NomFic, append=TRUE)
}		# boucle sur les fichiers

#	Calcul des moyennes

for	(iMeth in 1:NBM) {
	ScoreMoy["mean", "Sensib", iMeth] = mean(Score["Sensib", iMeth,])
	ScoreMoy["sd",   "Sensib", iMeth] = (var(Score["Sensib", iMeth,]))**0.5	
	ScoreMoy["mean", "Specif", iMeth] = mean(Score["Specif", iMeth,])
	ScoreMoy["sd",   "Specif", iMeth] = (var(Score["Specif", iMeth,]))**0.5
	
	cat ("Méthode ", Methods[iMeth], " Moyennes ", ScoreMoy[,,iMeth], "\n")
}

#	Permet d'afficher le tableau de droite


#	Matrices csv utilisées par Cytoscape pour dessiner "insilico_10_1" (ou 100_1)
#	Il faut supprimer la 1° ligne et la 1° colonne au préalable
#	Refaire passer, au préalable, le code précédent (boucle sur les fichiers), pour indFIC = 1 et Method = "RLD"   ("RLD" se traduit par "TLR" en anglais)

DessSol		<- matrix(nrow=sum(Solution[,])+1, ncol=4)		# La première ligne correspond aux noms ("Source", "Target" etc...)
DessCc1		<- matrix(nrow=sum(MatrCc1 [,]) + Score["F0", "RLD", 1] +1,  ncol=4)		# Cc1 : méthode RLD
nomNoeud	<- vector(length=NBN)			# Nom des noeuds = au numéro, par défaut

for (i in 1:NBN) {
	nomNoeud[i] = as.character(i)
}

DessSol[1,]		<- c("Source", "Target", "Interaction", "Value")
DessCc1[1,]		<- c("Source", "Target", "Interaction", "Value")
	
iSol 	= 2	
iCc1 	= 2
for (iRow in 1:NBN) {
	for (iCol in 1:NBN) {
		if (Solution[iRow, iCol] == 1) {
			DessSol[iSol,] = c(nomNoeud[iCol], nomNoeud[iRow], "d", 1)
			iSol = iSol +1
		}
		if (MatrCc1[iRow, iCol] == 1) {
			if (Solution[iRow, iCol] == 1) {
				DessCc1[iCc1,] = c(nomNoeud[iCol], nomNoeud[iRow], "d", 1)	# TP (True Positive)
				iCc1 = iCc1 +1			
			} else {
				DessCc1[iCc1,] = c(nomNoeud[iCol], nomNoeud[iRow], "fp", 1)	# FP (False Positive)
				iCc1 = iCc1 +1		
			}				
		} else {
			if (Solution[iRow, iCol] == 1) {
				DessCc1[iCc1,] = c(nomNoeud[iCol], nomNoeud[iRow], "fn", 1)	# FN (False Negative)
				iCc1 = iCc1 +1		
			}
		}
		# cat("Row ", iRow, " Col ", iCol, " MatrCc1 ", MatrCc1[iRow, iCol], " Solution ", Solution[iRow, iCol], " iCc1 ", iCc1, "\n") 
	}	# iCol
}		# iRow	

write.csv(DessSol, "insilico_size10_1_Sol.csv")
write.csv(DessCc1, "insilico_size10_1_Cc1.csv")

#	Les fichiers .csv permettent l'affichage des réseaux avec Cytoscape ("styles_MRA.xml")




#
#	Figure 5 : Knowledge impact on regression methods
#
#	Mettre à jour la variable "indNBN" pour balayer les différentes tailles (10 ou 100)
#	Pour une taille donnée, les 5 fichiers (indFIC) sont balayés automatiquement par une boucle for
#

Verbose	<- FALSE
NBFIC	<- 5				# Il y a 5 fichiers pour une taille donnée

indNBN	<- 1				# Valeurs possibles : 1 (pour 10 noeuds), 2 (pour 100 noeuds)

#Methods	<-  c("MRA-KD", "MRA-KO", "PVAL", "LASSO", "RLD", "LARS", "STEP-Fo", "STEP-Ba", "STEP-Bo")
											# Méthodes utilsables pour l'article BioInformatics  + "ARACNE", "CLR""
Methods		<-	c("RLD", "LASSO", "STEP-Fo")
NBM		<- length(Methods)										# Nb de méthodes à tester

vNBN	<- c(10, 100)																			# Nb noeuds
vNBK	<- c(2, 2)																				# Nb de types de perturbations (KD ou KO)
vTOP	<- c(0.2, 0.25)																			# "top" des arêtes à prendre
vNBT	<- c(100, 10)																			# Nb d'essais à faire
vZF		<- list(from=c(0,0), to=c(100,60), by=c(10,15))											# Paramètres (from, to, by) des séq. %ZF

vRoot	<- c("Documents/DreamChallenge10/", "Documents/DreamChallenge100/")						# Nom de la racine des fichiers
vRSMTxt	<- c("insilico_size10_",  "insilico_size100_")											# Texte à afficher dans les tables de regroupement
vRSMFic	<- c("MRA3_Resume_ZF10",  "MRA3_Resume_ZF100")											# Nom des fichiers de regroupement
vSolTxt	<- c("DREAM4_GoldStandard_InSilico_Size10_", "DREAM4_GoldStandard_InSilico_Size100_")	# Partie du nom des fichiers solution
ParNds	<- data.frame(vNBN, vNBK, vTOP, vNBT, vZF, vRoot, vRSMTxt, vRSMFic, vSolTxt)			# Tableau des données précédentes, correspondant aux nb. de noeuds possibles

NBN = ParNds$vNBN[indNBN]
NBK = ParNds$vNBK[indNBN]					# NBK ne correspond plus à du bruit ajouté, mais à des types de perturbations (KD ou KO)
NBT = ParNds$vNBT[indNBN]					# 100 essais pour chaque valeur de %ZF (réseaux à 10 noeuds) et 10 essais (réseaux à 100 noeuds)
ZF  		<- seq(ParNds$from[indNBN], ParNds$to[indNBN], by=ParNds$by[indNBN])				# Pourcentage de zéros forcés (%ZF)
NBZF		<- length(ZF)					# Nb. de pourcentages à tester

XMesNP		<- array(dim=c(NBN, NBK))		# Valeur mesurées, système non perturbé
XMesP		<- array(dim=c(NBN, NBN, NBK))	# Valeur mesurées, suite à une perturbation
Solution	<- array(dim=c(NBN, NBN))		# Solution : vraie structure du réseau
Solution2	<- array(dim=c(NBN, NBN))		# Idem "Solution", mais les n° de ligne et colonne correspondent aux n° de "levels" ("G1", "G10", "G2", "G3" etc)
Solution3	<- array(dim=c(NBN, NBN))		# Intermédiaire de calcul
VZero		<- vector(length=NBN*NBN)		# Les "vrais" zéros de la solution : (n° ligne-1)*NBN + (n° col-1)
szVZero		<- 0							# Nombre de vrais zéros (en excluant les zéros de la diagonale)
VZero0		<- vector(length=NBN*NBN)		# Recopie de VZero à chaque simulation
szVZero0	<- 0							# Recopie de szVZero à chaque simulation

MatR		<- array(dim=c(NBN, NBN, NBK))	# Matrice R
rij 		<- matrix(nrow=1, ncol=NBN-1)	# Calcul intermédiaire des rij
rL 			<- vector(length=NBN)			# Résultat de la méthode de Lasso, NBN-1 * NBN-1 (ie sol. de Yi = Ai * Xi)
MatrCc		<- array(dim=c(NBN, NBN))		# Matrice r calculée #  sauf par méthode PVAL    # selon la méthode MRA, RLD ou STEP ou LASSO
vrCc  		<- vector(length = NBN*NBN)		# MatrCc en vecteur
MatrCc1		<- array(dim=c(NBN, NBN))		# Matrice r résultat de l'algorithme
											# Elle est calculée en tronquant les éléments de MatrCc à 0 et 1 (0 sur la diagonale, comme Solution)

MatPVal		<- array(dim=c(NBN, NBN))		# Calcul direct de la p Value

VMean		<- array(dim=c(NBN, NBN))		# Valeur moyenne
VMin		<- array(dim=c(NBN, NBN))		# Valeur mini. de l'intervalle de confiance à 95%
VMax		<- array(dim=c(NBN, NBN))		# Valeur maxi. de l'intervalle de confiance à 95%

pX   <- array(dim=c(NBN-1, NBN-1, NBN))						# ligne (iRow), colonne (iCol), n° matrice (iMat)
pY   <- array(dim=c(NBN-1, NBN))							# ligne, n° matrice
gX   <- array(dim=c(NBK*(NBN-1), NBN-1, NBN))				# ligne (iRow), colonne (iCol), n° matrice (iMat)
gY   <- array(dim=c(NBK*(NBN-1), NBN))						# ligne, 1, n° matrice
gX1  <- array(dim=c(NBN-1, NBN-1, NBN))						# matrice gX tronquée pour ne prendre que les valeurs KD ou KO
gY1  <- array(dim=c(NBN-1, NBN))							# idem gX1
gXP	 <- array(dim=c(NBK*(NBN-1), NBN-1, NBN, NBT, NBZF))	# gX modifié par la connaissance de certains zéros
gXP1 <- array(dim=c(NBN-1, NBN-1, NBN, NBT, NBZF))			# idem gXP, tronqué comme gX1

Zero <- array(0, dim=c(NBN, NBN, NBT, NBZF))				# Un 1 en position i,j indique que le ri,j correspondant est nul

Moy				<- c("%ZF", "Mean", "sd")					# Moyennes à afficher
matZFSens.df 	<- data.frame(0, 0, 0)						# Data Frame pour utiliser gplot
matZFSpec.df 	<- data.frame(0, 0, 0)						# Data Frame pour utiliser gplot
mZFSens			<- array(dim=c(length(Moy), NBZF, NBM))		# %ZF, moyenne et écart-type de la sensibilité
mZFSpec			<- array(dim=c(length(Moy), NBZF, NBM))		# %ZF, moyenne et écart-type de la spécificité
dimnames(mZFSens)	<- list(Moy, ZF, Methods)
dimnames(mZFSpec)	<- list(Moy, ZF, Methods)

ScoreTest	<-	c("nbZF", "T+-", "T0", "F+-", "F0", "Sensib", "Specif")	# Scores mesurés
NBSC		<- length(ScoreTest)									# Nb. de scores testés
Score 		<- array(0, dim=c(NBSC, NBT, NBZF, NBM, NBFIC))	# Nb. zéros forcés, Vrais +- (TP), vrais 0 (TN), faux +- (FP), faux 0 (FN)
dimnames(Score)		<- list(ScoreTest, NULL, ZF, Methods, NULL)

NomFic <- paste(ParNds$vRoot[indNBN], "ZeroForce.txt", sep="")
for(indFIC in 1:NBFIC) {
	valw <- paste("\n\n\n", ParNds$vRSMTxt[indNBN], indFIC, sep="")
	write(valw, NomFic, append=TRUE)
	
	#	1/ Lecture des données
		
	val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_wildtype.tsv", sep="")
	XMesNPLu 	<- fread(val, data.table=F)						# "Documents/DreamChallenge10/insilico_size10_1_wildtype.tsv"
	val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_knockdowns.tsv", sep="")
	XMesPLu1 	<- fread(val, data.table=F) 					# ("Documents/DreamChallenge10/insilico_size10_1_knockdowns.tsv"
	val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_knockouts.tsv", sep="")
	XMesPLu2 	<- fread(val, data.table=F) 					# "Documents/DreamChallenge10/insilico_size10_1_knockouts.tsv"
	
	val			<- paste(ParNds$vRoot[indNBN], ParNds$vSolTxt[indNBN], indFIC, ".tsv", sep="")
	SolutLu  	<- fread(val, data.table=F) 					# "Documents/DreamChallenge10/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
	SolutLu$V1	<- as.factor(SolutLu$V1)
	SolutLu$V2	<- as.factor(SolutLu$V2)
	
	for (i in 1:NBN) {
		XMesNP[i,1] = XMesNPLu[1,i]
		XMesNP[i,2] = XMesNPLu[1,i]
		
		for (j in 1:NBN) {
			XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation à -50%	
												# ATTENTION les données de Dream Cha. sont transposées par rapport à mon écriture
			XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation à -100%
		}
	}
	
	Solution2[,] = 0
	for (i in 1:(NBN*(NBN-1))) {
	#	cat("i ", i, " V3 ", SolutLu[i,]$V3, " V1 ", SolutLu[i,]$V1, " : ", as.numeric(SolutLu[i,]$V1), " V2 ", SolutLu[i,]$V2, " : ", as.numeric(SolutLu[i,]$V2), "\n")
		if(SolutLu[i,]$V3 == 1) {
			Solution2[as.numeric(SolutLu[i,]$V1), as.numeric(SolutLu[i,]$V2)] = 1
		} else {
			break								# On suppose que tous les 1 sont au début du fichier
		}
	}
	
	numero <- as.numeric(sort(as.character(c(1:NBN))))
	
	for (i in 1:NBN) {
		Solution3[,numero[i]] = Solution2[,i]
	}
	for (i in 1:NBN) {
		Solution[numero[i],] = Solution3[i,]	# On retrouve la matrice de connexion classique (noeuds "G1", "G2", ... , "G10")
	}
	Solution <- t(Solution)						# Deuxième configuration 
												# correspond à la 1° interprétation en réalité : erreur sur l'écriture de la matrice
	
	for (iPaq in c(1:NBK)) {					# On prend toutes les données disponibles : Knock Down et Knock Out
		#	Les "données de mesure", fournies par Dream Challenge, ont été chargées dans XMesNP et XMesP précédemment	
		for (iPert in c(1:NBN)) {
			MatR[,iPert,iPaq]  <- 2 * (XMesP[,iPert,iPaq]-XMesNP[,iPaq]) / (XMesP[,iPert,iPaq]+XMesNP[,iPaq])
		}
		
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
			
			for (iRow in 1:(NBN-1)) {
				pY[iRow, iMat] = MatR[iMat, col[iRow], iPaq]
				
				for (iCol in 1:NBN-1) {
					pX[iRow, iCol, iMat] = MatR[col[iCol], col[iRow], iPaq]
				}
			}
		}
		
		gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,]  <- pX[, ,]
		gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ]   <- pY[, ]
	}
	
	#	2/ Recherche de la position des vrais zéros
	
	VZero[] <- 0
	ind 	<- 1		# indice de stockage dans VZero
	for (i in 1:NBN) {
		for (j in 1:NBN) {
			if ((Solution[i,j] == 0) && (i != j)) {
				VZero[ind] = (i-1)*NBN + (j-1)
				ind  = ind +1
			}
		}
	}
	szVZero = ind -1
	
	#	3/ Boucle sur le % de zéros et les simulations
	
	set.seed(12345)									# Pour pouvoir régénérer les mêmes séquences
	for (iZF in 1:NBZF) {							# Boucle sur les zéros à forcer
		nbZF = round(ZF[iZF] * szVZero / 100)		# Nb de zéros à forcer
		Score["nbZF",,iZF,,] = nbZF
		Nom <- as.character(ZF[iZF])				# Pour gplot : nom qui s'affichera en abscisse de la BP
		
		valw <- paste("\n", "% zéros forcés : ", ZF[iZF], " nb. zéros forcés : ", nbZF, sep="")
		write(valw, NomFic, append=TRUE)
	
		for (iT in 1:NBT) {							# NBT simulations
			nbZF0			<- nbZF
			Zero[,,iT,iZF]	<- 0
			VZero0 			<- VZero
			szVZero0		<- szVZero
			
			#	3.1/ Choix des zéros
			
			while (nbZF0 > 0) {
				indZF = round(runif(1, min=1, max=szVZero0))	# Position DU vrai zéro à prendre dans le tableau
				val = VZero0[indZF]
				if (Verbose) {
					cat("nbZF0", nbZF0, " indZF ", indZF, " szVZero0 ", szVZero0, " val ", val, " Row ", trunc(val/NBN)+1, " Col ", val%%NBN+1, "\n")
				}
				Zero[trunc(val/NBN)+1, val%%NBN+1, iT, iZF] = 1
				
				for (i in indZF:(szVZero0-1)) {		# On supprime la valeur utilisée dans VZero (tirage sans remise)
					VZero0[i] = VZero0[i+1]
				}	# tirage sans remise
				szVZero0 = szVZero0 - 1
				nbZF0 = nbZF0 - 1
			}		# while (nbZF0 > 0)
			
			#	3.2/ Prise en compte des zéros choisis
			
			gXP[,,, iT, iZF]  <- gX[,,]
			for (Row in 1:NBN) {
				for (Col in 1:NBN) {
					if (Zero[Row, Col, iT, iZF] == 1) {			# On sait que r Row,Col = 0
						if (Col > Row) {
							cc = Col -1
						} else {
							cc = Col
						}
						gXP[,cc, Row, iT, iZF] <- 1E6
						if (Verbose) {
							cat("gXP ", cc, Row, iT, iZF, "\n")
						}
					}
				}	# for Col
			}		# for Row
	
			#	3.3/ Recherche de MatrCc (dépend de la méthode choisie)
			
			for	(iMeth in 1:NBM) {
				Method  <- Methods[iMeth]
				
				cat("Fic :", indFIC, " %ZF ", ZF[iZF], " Test ", iT, " Méthode ", Method, "\n")
				
				if (Method %in% c("MRA-KD", "MRA-KO")) {
				#	Méthode MRA classique
					if (Method == "MRA-KD") {
						rSupp <- seq(NBN, NBK*(NBN-1), 1)	# On supprime les lignes 10 à 18 (resp. 100 à 198) pour KD
					} else {
						rSupp <- seq(1, NBN-1, 1)			# On supprime les lignes 1  à 9  (resp. 1 à 99) pour KO
					}
					
					gX1	 <- gX[-rSupp,,]
					gY1	 <- gY[-rSupp,]
					gXP1 <- gXP[-rSupp,,,,]
				
					for (iMat in 1:NBN) {
						col <- 1:NBN
						col <- col[-iMat]
							
						Form <- paste("gY1[,",iMat,"]~gXP1[,,",iMat,", iT, iZF]+0")					
						
						MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
						MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
					}
					MatrCc[is.na(MatrCc)]  <- 0
				}	# Méthodes MRA
				
				if (Method %in% c("RLD", "RLD OPTIM", "RLD OPTIM2")) {
				#	Méthode RLD classique (avec seuils)
				
					for (iMat in 1:NBN) {
						col <- 1:NBN
						col <- col[-iMat]
							
						Form <- paste("gY[,",iMat,"]~gXP[,,",iMat,", iT, iZF]+0")					
						
						MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
						MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
					}
					MatrCc[is.na(MatrCc)]  <- 0
				}	# Méthode RLD
				
				if (Method == "PVAL") {						#	p-Value à 5%
					SeuilPV	  <- 0.05
					MatrCc[,] <- 0
					for (iMat in 1:NBN) {
						col <- 1:NBN
						col <- col[-iMat]
							
						Form  <- paste("gY[,",iMat,"]~gXP[,,",iMat,", iT, iZF]+0")
						pQ 	  <- lm(as.formula(Form))
						
						MatPVal[iMat,col] 	<- broom::tidy(pQ)$p.value
						MatrCc [iMat,which(MatPVal[iMat,] < SeuilPV)] = 1
					}		# boucle iMat
				}	# Méthode PVAL
				
				if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
				#	Méthode "step" (suppression automatique de coefficients par optimisation de AIC : Akaike's Information Criterion)
				
					MatrCc[,] <- 0
					for (iMat in 1:NBN) {
						col <- 1:NBN
						col <- col[-iMat]
				
						Donnees  <- data.frame(gY[,iMat], gXP[,,iMat,iT,iZF])
						if (Verbose) {
							cat("iMat ", iMat, " Donnees ", colnames(Donnees), "\n")
						}
						Donnees	 <- rename(Donnees, "Y" ="gY...iMat.")		# Le suivi est plus clair comme cela
						cc 		 <- colnames(Donnees)						# Nom des colonnes de "Données". La première est "Y"
						cc		 <- cc[-1]									# Il reste le nom des coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
						colnames(rij) <- cc									# step indique le nom des colonnes correspondant aux coefficients conservés
								
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
		
						ll  	<- length(forward$coefficients)				# Nombre de coefficients conservés par step
						nn		<- names(forward$coefficients)				# Nom de ces coefficients
				
						rij[1,]	<- 0
						if (ll >= 2) {
							for (i in 2:ll) {								# nn[1] = "(Intercept)"
								rij[1, nn[i]] <- forward$coefficients[nn[i]]
							}
						}
						if (Verbose) {
							cat (" iZF ", iZF, " iT ", iT, " iMat ", iMat, " ll ", ll, " nn ", nn, " rij ", rij, "\n")
						}
						MatrCc[iMat, col]	<- rij[1,]
						MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
					}		# boucle iMat
				}	# Méthodes STEP
					
				if (Method == "LASSO") {
				#	Méthode "Lasso", choix automatique de Lambda
				
					MatrCc[,] <- 0		# Initialisation
					for (iMat in 1:NBN) {
						col <- 1:NBN
						col <- col[-iMat]
						
						exc  <- vector(length=0)						# colonnes à exclure
						for (i in 1:NBN) {
							if (Zero[iMat, i, iT, iZF] == 1) {			# On sait que r iMat,i = 0
								if (i > iMat) {
									cc = i -1
								} else {
									cc = i
								}
								exc <- cbind(exc, cc)
							}
						}	# for i					
				
						if (length(exc) < (NBN-1)) {
							cv_model 	<- cv.glmnet(gXP[,,iMat,iT,iZF], gY[,iMat], alpha = 1, exclude = exc)			# Fit lasso regression model using k-fold cross-validation
							best_lambda <- cv_model$lambda.min
							best_model 	<- glmnet(gXP[,,iMat,iT,iZF], gY[,iMat], alpha = 1, lambda = best_lambda, exclude = exc)	# View coefficients of best model
							rL			<-	coef(best_model)								
			
							MatrCc[iMat, col]	<- rL[2:NBN]
						} else {
							MatrCc[iMat, col]	<- 0		# tous les éléments de la ligne iMat sontnuls
						}			
					}		
				}	# Méthode  LASSO
				
				if (Method == "LARS") {
					MatrCc[,] <- 0		# Initialisation
					for (iMat in 1:NBN) {
						col <- 1:NBN
						col <- col[-iMat]
				
						pp=lars(gXP[,,iMat,iT,iZF], gY[,iMat], type="lar")
						rL = coef(pp)
				
						MatrCc[iMat, col]	<- rL[2:NBN]
					}		
				}	# Méthode  LARS		
				
				#	3.4/ 	Recherche de la solution
				#	Numérisation du résultat		
				
				MatrMax = max(abs(MatrCc[,]))
				if (Method == "RLD") {
					Th = 0.25*MatrMax
				}	else if (Method %in% c("MRA-KD", "MRA-KO")) {	
					NB1 = trunc(ParNds$vTOP[indNBN]*NBN*NBN)
					vrCc[]  <- MatrCc[]
					vrCc1	<- sort(abs(vrCc), decreasing=TRUE)
					Th	= vrCc1[min(NB1+1, NBN*NBN)]
				}	else {
					Th = 0
				}
					
				for (iRow in 1:NBN) {
					for (iCol in 1:NBN) {
						if (abs(MatrCc[iRow, iCol]) > Th) {
							MatrCc1[iRow, iCol] = 1
						} else {
							MatrCc1[iRow, iCol] = 0
						}
					}
					MatrCc1[iRow, iRow] = 0
				}		
			
				#	3.5/	Test de la solution
				
				Score[2:NBSC,iT,iZF,iMeth,indFIC] = 0			
	
				for (iRow in 1:NBN) {
					for (iCol in 1:NBN) {
						if (MatrCc1[iRow, iCol] != 0) {
							if (Solution[iRow, iCol] != 0) {
								Score["T+-",iT,iZF,iMeth,indFIC] <- Score["T+-",iT,iZF,iMeth,indFIC] +1
							} else {
								Score["F+-",iT,iZF,iMeth,indFIC] <- Score["F+-",iT,iZF,iMeth,indFIC] +1
							}
						} else {
							if (Solution[iRow, iCol] == 0) {
								Score["T0",iT,iZF,iMeth,indFIC] <- Score["T0",iT,iZF,iMeth,indFIC] +1
							} else {
								Score["F0",iT,iZF,iMeth,indFIC] <- Score["F0",iT,iZF,iMeth,indFIC] +1
							}
						}		
					}
				}
				
				P = sum(Solution)			# Nb. vrais positifs
				N = NBN*(NBN-1) - P			# Nb. vrais négatifs
				
				Score["Sensib",iT,iZF,iMeth,indFIC] = Score["T+-",iT,iZF,iMeth,indFIC] / P
											# (Score["T+-", iMeth] + Score["F0", iMeth])
				Score["Specif",iT,iZF,iMeth,indFIC] = (Score["T0",iT,iZF,iMeth,indFIC] - NBN) / N
											# (Score["T0", iMeth] + Score["F+-", iMeth] - NBN)
				
				Err	<- (Score["F+-",iT,iZF,iMeth,indFIC] + Score["F0",iT,iZF,iMeth,indFIC])/(NBN*(NBN-1))
			
				matZFSens.df 	<- rbind(matZFSens.df, 
						c(Score["T+-",iT,iZF,iMeth,indFIC] / (Score["T+-",iT,iZF,iMeth,indFIC]+Score["F0",iT,iZF,iMeth,indFIC]), Method, Nom))
				matZFSpec.df 	<- rbind(matZFSpec.df, 
						c((Score["T0",iT,iZF,iMeth,indFIC]-NBN) / (Score["T0",iT,iZF,iMeth,indFIC]+Score["F+-",iT,iZF,iMeth,indFIC]-NBN), Method, Nom))			
			
				if (Verbose) {
					cat(as.character(ParNds$vRSMTxt[indNBN]), indFIC, sep="")
					cat(" Méthode ", Method, " Th ", Th, " Score ", Score[,iT,iZF,iMeth,indFIC], " Err ", Err, "\n")
				}
			
				valw <- paste("Méthode : ", Method, "Th", Th, " iT ", iT, "Score", paste(Score[,iT,iZF,iMeth,indFIC], collapse=";"), "Err", Err, sep=";")
				write(valw, NomFic, append=TRUE)
			}	# boucle iMeth (méthodes)
		}		# boucle iT  (100 ou 10 essais)
	}			# boucle iZF (% de zéros connus)	
}				# boucle sur les fichiers

for (iZF in 1:NBZF) {
	for (iMeth in 1:NBM) {
		mZFSens["%ZF",iZF, iMeth] 	= ZF[iZF]
		mZFSens["Mean",iZF, iMeth] 	= mean(Score["T+-",,iZF, iMeth,] / (Score["T+-",,iZF, iMeth,]+Score["F0",,iZF, iMeth,]))
		mZFSens["sd",  iZF, iMeth ] = sd  (Score["T+-",,iZF, iMeth,] / (Score["T+-",,iZF, iMeth,]+Score["F0",,iZF, iMeth,])) 
		mZFSpec["%ZF", iZF, iMeth]  = ZF[iZF]
		mZFSpec["Mean",iZF, iMeth]  = mean((Score["T0",, iZF, iMeth,]-NBN) / (Score["T0",, iZF, iMeth,]+(Score["F+-",, iZF, iMeth,]-NBN)))
		mZFSpec["sd",  iZF, iMeth]  = sd  ((Score["T0",, iZF, iMeth,]-NBN) / (Score["T0",, iZF, iMeth,]+(Score["F+-",, iZF, iMeth,]-NBN)))
		
		cat("ZF ", ZF[iZF], " Méthode ", Methods[iMeth], " mSens ", mZFSens["Mean",iZF,iMeth], " sdSens ", mZFSens["sd",iZF,iMeth], 
			" mSpec ", mZFSpec["Mean",iZF,iMeth], " sdSpec ", mZFSpec["sd",iZF,iMeth],"\n")
		valw <- paste("%ZF ", ZF[iZF], " Méthode ", Methods[iMeth], " mSens ", mZFSens["Mean",iZF,iMeth], " sdSens ", mZFSens["sd",iZF,iMeth], 
					  " mSpec ", mZFSpec["Mean",iZF,iMeth], " sdSpec ", mZFSpec["sd",iZF,iMeth], sep=";")
		write(valw, NomFic, append=TRUE)
	}			# boucle sur les méthodes
}				# boucle iZF (% de zéros connus)

matZFSens.df 	 		<- matZFSens.df[-1,,]			# La 1° ligne, créée à la déclaration, vaut (0,0, 0)
names(matZFSens.df)		<- c("Sensib", "Meth", "ZF")
if (ZF[NBZF] == 100) {
	matZFSens.df$ZF		<- fct_relevel(matZFSens.df$ZF, "100", after=Inf)	#  100% est mis à la fin (R suit l'ordre lexicographique par défaut !) 
}
matZFSens.df$Sensib		<- as.numeric(matZFSens.df$Sensib)
matZFSens.df$ZF			<- as.factor (matZFSens.df$ZF)

png(filename="MRA3_5a.png")			# fichier png contenant les boxplots Sensibilité vs %ZF  -- Réseaux 10 noeuds
#png(filename="MRA3_5c.png")			# fichier png contenant les boxplots Sensibilité vs %ZF  -- Réseaux 100 noeuds
ggplot(matZFSens.df, aes(x=ZF, y=Sensib, fill=Meth), ylimit=c(0.35,1)) + 
	geom_boxplot(outlier.colour="red") + 
	labs(title="Sensibility vs. % values forced to 0", x="% values forced to 0", y="Sensibility") +
    scale_colour_manual(name = "Meth",
		values = c("RLD" = "green", "LASSO" = "red", "STEP-Fo" = "orange"),
		labels = c("RLD" = "TLR", "LASSO" = "LASSO", "STEP-Fo" = "STEP-Fo"))
dev.off()	
#		Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL))
#		Erreur de R. Pour la fixer, faire dev.off()
#

matZFSpec.df 	 		<- matZFSpec.df[-1,,]			# La 1° ligne, créée à la déclaration, vaut (0,0, 0)
names(matZFSpec.df)		<- c("Sensib", "Meth", "ZF")
if (ZF[NBZF] == 100) {
	matZFSpec.df$ZF		<- fct_relevel(matZFSpec.df$ZF, "100", after=Inf)	#  100% est mis à la fin (R suit l'ordre lexicographique par défaut !) 
}
matZFSpec.df$Sensib		<- as.numeric(matZFSpec.df$Sensib)
matZFSpec.df$ZF			<- as.factor (matZFSpec.df$ZF)

png(filename="MRA3_5b.png")			# fichier png contenant les boxplots Spécificité vs %ZF  -- Réseaux 10 noeuds
#png(filename="MRA3_5d.png")			# fichier png contenant les boxplots Sensibilité vs %ZF  -- Réseaux 100 noeuds
ggplot(matZFSpec.df, aes(x=ZF, y=Sensib, fill=Meth)) + 
	geom_boxplot(outlier.colour="red") + 
	labs(title="Specificity vs. % values forced to 0", x="% values forced to 0", y="Specificity") +
    scale_colour_manual(name = "Meth",
		values = c("RLD" = "green", "LASSO" = "red", "STEP-Fo" = "orange"),
		labels = c("RLD" = "TLR", "LASSO" = "LASSO", "STEP-Fo" = "STEP-Fo"))
dev.off()	

# 	Pour obtenir les images du bas (réseau 100 noeuds), remplacer indNBN = 1 par indNBN = 2 et permuter les noms de fichiers.





#
#	Sup. Infos 2
#
NBN 	<- 6										# Nombre de noeuds
Verbose	<- FALSE									# Si TRUE, afficher certains commentaires

Noise		<-	c(0.001, 0.005, 0.01)				# 3 niveaux de bruit
Replic		<-  c(3, 5)								# 3 et 5 réplicats techniques
NPerturb 	<-  1.5									# Perturbation +50%

Methods		<-  c("Th", "BS IC95", "PVal 5", "Lasso", "RLD", "Step F", "Step B", "Step BO", "Step F_0", "Step B_0", "Step BO_0")		
																						# Méthodes testées
Steps		<- 	c("Step F", "Step B", "Step BO", "Step F_0", "Step B_0", "Step BO_0")	# Méthodes "Step" (forward, backward, both, avec ou sans intercept forcé à 0)																		
ScoreTest	<-	c("T+", "T-", "T0", "F+", "F-", "F0", "Sensib", "Specif")				# Scores mesurés

NBT 	<- length(Noise)
NBK		<- 100										# 100 tirages de bruit
NBR		<- length(Replic)							# Les essais de réplicats
NBRM	<- max(Replic)								# Nb max de réplications à effectuer
NBM		<- length(Methods)							# Nb de méthodes à tester
NBRST	<- length(Steps)							# Nb de méthodes STEP
NBSC	<- length(ScoreTest)						# Nb. de scores testés	
		
Form	<- vector(length=NBN)						# Formule décrivant la régression pour le calcul de MatrCc
XMesNP	<- array(dim=c(NBN, NBK, NBT))				# Valeur mesurées, système non perturbé
XMesP	<- array(dim=c(NBN, NBN, NBK, NBT))			# Valeur mesurées, suite à une perturbation
MatR	<- array(dim=c(NBN, NBN, NBK))				# Matrice R
MatrCc	<- array(dim=c(NBN, NBN, NBK, NBT))			# Matrice r calculée pour chaque valeur de bruit (bootstrap)
		
pX 	<- array(dim=c(NBN-1, NBN-1, NBN))				# ligne (iRow), colonne (iCol), n° matrice (iMat)
pY 	<- array(dim=c(NBN-1, NBN))						# ligne, n° matrice
gX 	<- array(dim=c(NBRM*(NBN-1), NBN-1, NBN))		# ligne (iRow), colonne (iCol), n° matrice (iMat)
gY 	<- array(dim=c(NBRM*(NBN-1), NBN))				# ligne, 1, n° matrice
gX1	<- array(dim=c(NBRM*(NBN-1), NBN-1, NBN))		# Copie de gX, limitée au nombre de tirages correspondant au nb. de réplicats
gY1	<- array(dim=c(NBRM*(NBN-1), NBN))				# Copie de gY, idem

rL 		<- matrix(nrow=NBN, ncol=1)					# Résultat de la méthode de Lasso, NBN-1 * NBN-1 (ie sol. de Yi = Ai * Xi)
rLasso  <- array(dim=c(NBN, NBN, NBR, NBT))			# Intégration des résultats Xi, ajout des -1 sur la diagonale et transpo ==> donne la matrice "r"

MatrStp <- array(dim=c(NBN, NBN, NBR, NBT, NBRST))	# Calcul de MatrCc par les méthodes STEP
rij 	<- matrix(nrow=1, ncol=NBN-1)				# Calcul intermédiaire des rij
		
VMean	<- array(dim=c(NBN, NBN, NBR, NBT))			# Moyenne calculée
VMin	<- array(dim=c(NBN, NBN, NBR, NBT))			# Valeur mini. de l'intervalle de confiance à 95%
VMax	<- array(dim=c(NBN, NBN, NBR, NBT))			# Valeur maxi. de l'intervalle de confiance à 95%
VCc		<- array(dim=c(NBN, NBN, NBR, NBT))			# MatrCc calculée avec les réplicats (iRep)

ValNom		<- vector(length=(NBN-1)*NBN)			# Nom (n) correspondant à (i,j)
ValMean		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Moyenne des réplicats
ValMin		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Min IC95 des réplicats (correspond à p_value à 5%)
ValMax		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# Max IC95 des réplicats (correspond à p_value à 5%)
ValCc		<- array(dim=c((NBN-1)*NBN, NBR, NBT))	# MatrCc calculée avec les réplicats (iRep) (indices homogènes à ValMean etc ...)
ValStp		<- array(dim=c((NBN-1)*NBN, NBR, NBT, NBRST))	# MatrStp calculée avec les réplicats (iRep) (indices homogènes à ValMean etc ...)
ValMeanBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Moyenne des 100 valeurs (cf Boxplot)
ValMinBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Quantile à 2,5%  des 100 valeurs
ValMaxBP	<- array(dim=c((NBN-1)*NBN, NBT))		# Quantile à 97,5% des 100 valeurs

Result 	<- array(dim=c((NBN-1)*NBN, NBR, NBT, NBM))	# Résultats obtenus
Score 	<- array(dim=c(NBSC, NBR, NBT, NBM-1))		# Vrais +, vrais - vrais 0, faux +, faux - faux 0    --- par comparaison à la matrice "r" théorique
											
dimnames(Result)	<- list(NULL, NULL, NULL, Methods)				# Nom des méthodes correspondant aux résultats
dimnames(Score)		<- list(ScoreTest, NULL, NULL, Methods[2:NBM])	# Nom des scores testés
									
# Recherche de la formule décrivant la régression pour le calcul de MatrCc (Intercept nul)

for (i in 1:NBN) {
	col <- 1:NBN
	col <- col[-i]
	
	Form[i] <- paste("MatR[", i, ",col,iPaq]-", sep="")			# Formule décrivant la régression
	for (j in col) {
		Form[i] <- paste(Form[i], "MatR[", sep="A")
		Form[i] <- paste(Form[i], j, sep="")
		Form[i] <- paste(Form[i], ",col", sep="")
		Form[i] <- paste(Form[i], ",iPaq]", sep= "")
	}
	
	Form[i] <- paste(Form[i], "+0", sep="")
	Form[i] <- str_replace_all(Form[i], pattern = "-A", replacement = "~")
	Form[i] <- str_replace_all(Form[i], pattern = "A", replacement = "+")
	if (Verbose) {
		cat("BOOTSTRAP i ", i, " Form ", Form[i], "\n")
	}
}

set.seed(12345)								# Pour pouvoir régénérer les mêmes séquences

for (iTest in c(1:NBT)) {
	# Pour chaque niveau de bruit (iTest), on effectue 100 tirages (NBK) et on calcule 100 matrices "r" (chaque matrice correspond à un tirage)
	# Puis on effectue des statistiques sur les erreurs obtenues : moyenne et sd des ri,j, ainsi que l'erreur moyenne et son écart-type
	# La régression linéaire porte sur une valeur seule. Cette méthode (6 matrices 5*5) est équivalente à l'inversion de la matrice 6*6
	
	sd <- Noise[iTest]*XMesNPMoy			# écart-type du bruit

	for (iPaq in c(1:NBK)) {
		#	Entrée des "données de mesure" bruitées
		ss <- multiroot (f=F, start= st)	# Mesure non perturbée  - utilise le package rootSolve
		XMesNP[,iPaq,iTest]	<- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		
		V4 <- NPerturb*V4		# P1
		ss <- multiroot (f=F, start= st)
		XMesP[,1,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V4 <- V4/NPerturb		# Retour à la valeur initiale
		
		V3 <- NPerturb*V3		# P2
		ss <- multiroot (f=F, start= st)
		XMesP[,2,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V3 <- V3/NPerturb		# Retour à la valeur initiale
		
		V8 <- NPerturb*V8		# P3
		ss <- multiroot (f=F, start= st)
		XMesP[,3,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V8 <- V8/NPerturb		# Retour à la valeur initiale
		
		V7 <- NPerturb*V7		# P4
		ss <- multiroot (f=F, start= st)
		XMesP[,4,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V7 <- V7/NPerturb		# Retour à la valeur initiale
		
		V12 <- NPerturb*V12		# P5
		ss <- multiroot (f=F, start= st)
		XMesP[,5,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V12 <- V12/NPerturb		# Retour à la valeur initiale
		
		V11 <- NPerturb*V11		# P6
		ss <- multiroot (f=F, start= st)
		XMesP[,6,iPaq,iTest] <- c (ss$root[1], ss$root[2], ss$root[3], ss$root[4], ss$root[5], ss$root[6]) + rnorm(NBN, mean=0, sd=sd)
		V11 <- V11/NPerturb		# Retour à la valeur initiale	
	
		for (iMat in c(1:NBN)) {
			MatR[,iMat,iPaq]  <- 2 * (XMesP[,iMat,iPaq,iTest]-XMesNP[,iPaq,iTest]) / (XMesP[,iMat,iPaq,iTest]+XMesNP[,iPaq,iTest])
		}
		
		for (i in 1:NBN) {
			col <- 1:NBN
			col <- col[-i]
					
			MatrCc[i, col, iPaq, iTest]   <- (lm(as.formula(Form[i])))$coefficients
			MatrCc[i, i,   iPaq, iTest]   <- -1
		}

		if (iPaq <= NBRM) {
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
				
				for (iRow in 1:(NBN-1)) {
					pY[iRow, iMat] = MatR[iMat, col[iRow], iPaq]			
					for (iCol in 1:NBN-1) {
						pX[iRow, iCol, iMat] = MatR[col[iCol], col[iRow], iPaq]
					}
				}
			}
			
			gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,]  <- pX[, ,]
			gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq),  ]  <- pY[,  ]
		}
	}			# Essais de bruits (iPaq)
	
	#	Application des différentes méthodes
	
	for (iRep in 1:NBR) {							# 3 ou 5 réplicats
		NRep 	<- Replic[iRep]						# Nb de réplicats
		gX1		<- gX
		gY1		<- gY
			
		if (NRep < NBRM) {
			v	<- seq(NRep*(NBN-1)+1, NBRM*(NBN-1), 1)				# Lignes à éliminer (cas où il y a moins de réplicats que NBRM)
			gX1	<- gX1[-v,,]
			gY1	<- gY1[-v,]
		}
		
		MatrStp[,, iRep, iTest,] <- 0
				
		indN	<- 1
		for (iRow in 1:NBN) {
			col <- 1:NBN
			col <- col[-iRow]
			
			#	Méthode BS IC95  (ie. MRA de base) et p-Value (5%)
			Form1 <- paste("gY1[,",iRow,"]~gX1[,,",iRow,"]+0")		# Pour comparaison avec la régression linéaire directe
													# Intercept forcé à 0
													# On les valeurs correspondant à Nrep (3 ou 5) tirages
			if (Verbose) {
				cat(NRep, " REPLICATS iTest ", iTest, " iRow ", iRow, " Form ", Form1, "\n")
			}

			pQ <- lm(as.formula(Form1))
			VMean[iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate
			VMin [iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate - 1.96*broom::tidy(pQ)$std.error
			VMax [iRow, col, iRep, iTest]  <- broom::tidy(pQ)$estimate + 1.96*broom::tidy(pQ)$std.error
			VCc	 [iRow, col, iRep, iTest]  <- lm(pQ)$coefficients
																# Intervalle de confiance à 95%
			
			#	Méthode de LASSO (avec choix automatique de Lambda)
			cv_model 	<- cv.glmnet(gX1[,,iRow], gY1[,iRow], alpha = 1)	# Fit lasso regression model using k-fold cross-validation
			best_lambda <- cv_model$lambda.min
			best_model 	<- glmnet(gX1[,,iRow], gY1[,iRow], alpha = 1, lambda = best_lambda)	# View coefficients of best model
			rL			<-	coef(best_model)
			rLasso[iRow, iRow, iRep, iTest] <- -1						
			rLasso[col,iRow, iRep, iTest]	<- rL[2:NBN,1]
			rLasso[,, iRep, iTest] <- t(rLasso[,,iRep, iTest])	# on trouve ici l'équivalent de la matrice "r"				
			
			#	Méthodes STEP
			Donnees  <- data.frame(gY1[,iRow], gX1[,,iRow])
			Donnees	 <- rename(Donnees, "Y" ="gY1...iRow.")		# Le suivi est plus clair comme cela
			cc 		 <- colnames(Donnees)						# Nom des colonnes de "Données". La première est "Y"
			cc		 <- cc[-1]									# Il reste le nom des coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
			colnames(rij) <- cc									# step indique le nom des colonnes correspondant aux coefficients conservés			
	
			for (iStep in 1:NBRST) {							# méthodes "step" (fonction step() -- package "stats")			
				switch(iStep, 
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},		# Forward
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward
					   {intercept_only <- lm(Y ~ 1, data=Donnees)
						all 	<- lm(Y ~ ., data=Donnees)
						forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)},		# Both
					   {intercept_only <- lm(Y ~ 1+0, data=Donnees)
						all 	<- lm(Y ~ .+0, data=Donnees)
						forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)},		# Forward Intercept=0
					   {intercept_only <- lm(Y ~ 1+0, data=Donnees)
						all 	<- lm(Y ~ .+0, data=Donnees)
						forward <- step(all, direction='backward', scope=formula(all), trace=0)},				# Backward Intercept=0
					   {intercept_only <- lm(Y ~ 1+0, data=Donnees)
						all 	<- lm(Y ~ .+0, data=Donnees)
						forward <- step(intercept_only, direction='both', scope=formula(all), trace=0)})		# Both Intercept=0
					
				ll  	<- length(forward$coefficients)				# Nombre de coefficients conservés par step
				nn		<- names(forward$coefficients)				# Nom de ces coefficients
				if (Verbose) {
					cat ("STEP ", iStep, " iRow ", iRow, " Coefs. retenus ", nn, " Valeurs ", forward$coefficients[nn[i]], "\n")
				}
	
				rij[1,]	<- 0
				if (ll >= 2) {
					for (i in 2:ll) {								# nn[1] = "(Intercept)"
						rij[1, nn[i]] <- forward$coefficients[nn[i]]
					}
				}
	
				MatrStp[iRow, col,  iRep, iTest, iStep]		<- rij[1,]
				MatrStp[iRow, iRow, iRep, iTest, iStep]  	<- -1
			}
			
			for (iCol in col) {	
				ValNom[indN] 	<- 6*(iRow-1)+iCol
				
				#	Valeurs théoriques
				if (matrTh[iRow, iCol] > 0) {
					Result [indN, iRep, iTest, 1] = 1				# 1° méthode : valeurs théoriques
				}	else if (matrTh[iRow, iCol] == 0) {
					Result [indN, iRep, iTest, 1] = 0	
				}	else {
					Result [indN, iRep, iTest, 1] = -1
				}									
				
				#	Réplicats
				ValMean[indN, iRep, iTest] 	<-	VMean[iRow, iCol, iRep, iTest]
				ValMin [indN, iRep, iTest] 	<-	VMin [iRow, iCol, iRep, iTest]
				ValMax [indN, iRep, iTest] 	<-	VMax [iRow, iCol, iRep, iTest]
				ValCc  [indN, iRep, iTest] 	<-	VCc  [iRow, iCol, iRep, iTest]
				for (iStep in 1:NBRST) {
					ValStp [indN, iRep, iTest, iStep] 	<-	MatrStp [iRow, iCol, iRep, iTest, iStep]
				}
				
				#	Calcul des quantiles (bootstrap)
				ValMeanBP[indN, iTest] 			<- mean(MatrCc[iRow,iCol,,iTest])
				qq	<- quantile(MatrCc[iRow,iCol,,iTest], probs=c(0.025,0.975))		# Donne la plage de IC à 95%
				ValMinBP [indN, iTest] 			<- qq[1]
				ValMaxBP [indN, iTest] 			<- qq[2]					
				
				if (ValMinBP[indN, iTest] > 0) {
					Result [indN, iRep, iTest, 2] = 1			# 2° méthode : bootstrap sur les valeurs min et max de IC95
				}	else if (ValMaxBP[indN, iTest] < 0) {
					Result [indN, iRep, iTest, 2] = -1	
				}	else {
					Result [indN, iRep, iTest, 2] = 0
				}

				#	p-Value à 5%
				if (ValMin[indN, iRep, iTest] > 0) {
					Result [indN, iRep, iTest, 3] = 1			# 3° méthode : p-Value à 5%
				}	else if (ValMax[indN, iRep, iTest] < 0) {
					Result [indN, iRep, iTest, 3] = -1	
				}	else {
					Result [indN, iRep, iTest, 3] = 0
				}
				
				#	Méthode de Lasso
				if (rLasso[iRow, iCol, iRep, iTest] > 0) {
					Result [indN, iRep, iTest, 4] = 1			# 4° méthode : méthode de Lasso
				}	else if (rLasso[iRow, iCol, iRep, iTest] == 0) {
					Result [indN, iRep, iTest, 4] = 0	
				}	else {
					Result [indN, iRep, iTest, 4] = -1
				}
				
				indN <- indN +1
			}	# iCol
		}		# iRow
		
		ThRLD	<- 0.25 * max(ValCc[, iRep, iTest])
		for (ind in 1:(indN-1)) {
			if (ValCc[ind, iRep, iTest] > ThRLD) {
				Result [ind, iRep, iTest, 5] = 1				# 5° méthode : RLD à seuil
			}	else if (ValCc[ind, iRep, iTest] < -ThRLD) {
				Result [ind, iRep, iTest, 5] = -1	
			}	else {
				Result [ind, iRep, iTest, 5] = 0
			}
		}		# ind				
			
		for (iStep in 1:NBRST) {
			for (ind in 1:(indN-1)) {
				if (ValStp[ind, iRep, iTest, iStep] > 0) {
					Result [ind, iRep, iTest, 5+iStep] = 1		# 6° à 11° méthodes : STEP (for, back, both, for0, back0, both0)
				}	else if (ValStp[ind, iRep, iTest, iStep] < 0) {
					Result [ind, iRep, iTest, 5+iStep] = -1	
				}	else {
					Result [ind, iRep, iTest, 5+iStep] = 0
				}
			}		# ind
		}		# iStep
	}			# Réplicats
}				# Tests (niveaux de bruit)

#	Scores obtenus

for (iTest in 1:NBT)		{
	for (iRep in 1:NBR)		{
		for (iMet in 2:NBM)	{						# iMet = 1 correspond au résultat exact (Th éorique)
			Score[,iRep, iTest, iMet-1]	<- 0		# mise à 0 des compteurs
			
			for (ind in 1:(indN-1)) {				# ind balaye tous les couples (i,j), i # j
				if (Result[ind, iRep, iTest, iMet]  > 0) {
					if (Result[ind, iRep, iTest, 1] > 0) {		# le dernier indice (1) fait référence à la valeur exacte
						Score["T+", iRep, iTest, iMet-1] <- Score["T+", iRep, iTest, iMet-1] +1
					} else {
						Score["F+", iRep, iTest, iMet-1] <- Score["F+", iRep, iTest, iMet-1] +1
					}
				}
				if (Result[ind, iRep, iTest, iMet]  < 0) {
					if (Result[ind, iRep, iTest, 1] < 0) {
						Score["T-", iRep, iTest, iMet-1] <- Score["T-", iRep, iTest, iMet-1] +1
					} else {
						Score["F-", iRep, iTest, iMet-1] <- Score["F-", iRep, iTest, iMet-1] +1
					}
				}
				if (Result[ind, iRep, iTest, iMet]  == 0) {
					if (Result[ind, iRep, iTest, 1] == 0) {
						Score["T0", iRep, iTest, iMet-1] <- Score["T0", iRep, iTest, iMet-1] +1
					} else {
						Score["F0", iRep, iTest, iMet-1] <- Score["F0", iRep, iTest, iMet-1] +1
					}
				}				
			}	# ind   (Test des valeurs calculées pour tous les couples i,j)
			
			P	<- sum (abs(Result[, iRep, iTest, 1]))		# Nb. d'elts non nuls (+ et -)
			N 	<- NBN*(NBN-1) - P							# Nb. d'elts nuls
			Score["Sensib", iRep, iTest, iMet-1] <- (Score["T+", iRep, iTest, iMet-1]+Score["T-", iRep, iTest, iMet-1]) / P
			Score["Specif", iRep, iTest, iMet-1] <- (Score["T0", iRep, iTest, iMet-1]) / N
		}		# iMet	(Les différentes méthodes)
	}			# iRep	(Nb. de réplicats techniques)
}				# iTest (Niveaux de bruit)





#
#	Sup. Infos 3
#
#	Tracé des courbes ROC correspondant aux différents "challenges" et aux différentes méthodes.
#	Mettre à jour les variables "indNBN" et "indFIC" pour balayer les différentes tailles (10 ou 100) et les 5 fichiers.
#

Verbose	<- FALSE
indNBN	<- 1				# Valeurs possibles : 1 (pour 10 noeuds), 2 (pour 100 noeuds)
indFIC	<- 1				# Valeurs possibles : 1 à 5
# Pour balayer tous les fichiers, on ne peut pas utiliser une boucle car le stockage des courbes ne fonctionne pas alors.

ScoreTest	<-	c("T+-", "T0", "F+-", "F0", "Sensib", "Specif")	# Scores mesurés
NBSC	<- length(ScoreTest)									# Nb. de scores testés
Methods	<-  c("MRA-KD", "MRA-KO", "PVAL", "LASSO", "RLD", "STEP-Fo", "STEP-Ba", "STEP-Bo")
NBM		<- length(Methods)

vNBN	<- c(10, 100)																			# Nb noeuds
vNBK	<- c(2, 2)
vTOP	<- c(0.2, 0.25)																			# Nb de types de perturbations (KD ou KO)

vRoot	<- c("Documents/DreamChallenge10/", "Documents/DreamChallenge100/")						# Nom de la racine des fichiers
vRSMTxt	<- c("insilico_size10_",  "insilico_size100_")											# Texte à afficher dans les tables de regroupement
vRSMFic	<- c("MRA3_Resume_ZF10",  "MRA3_Resume_ZF100")											# Nom des fichiers de regroupement
vSolTxt	<- c("DREAM4_GoldStandard_InSilico_Size10_", "DREAM4_GoldStandard_InSilico_Size100_")	# Partie du nom des fichiers solution

ParNds	<- data.frame(vNBN, vNBK, vTOP, vRoot, vRSMTxt, vRSMFic, vSolTxt)						# Tableau des données précédentes, correspondant aux nb. de noeuds possibles

NBN = ParNds$vNBN[indNBN]
NBK = ParNds$vNBK[indNBN]					# NBK ne correspond plus à du bruit ajouté, mais à des types de perturbations (KD ou KO)

XMesNP		<- array(dim=c(NBN, NBK))		# Valeur mesurées, système non perturbé
XMesP		<- array(dim=c(NBN, NBN, NBK))	# Valeur mesurées, suite à une perturbation
Solution	<- array(dim=c(NBN, NBN))		# Solution : vraie structure du réseau
Solution2	<- array(dim=c(NBN, NBN))		# Idem "Solution", mais les n° de ligne et colonne correspondent aux n° de "levels" ("G1", "G10", "G2", "G3" etc)
Solution3	<- array(dim=c(NBN, NBN))		# Intermédiaire de calcul

MatR		<- array(dim=c(NBN, NBN, NBK))	# Matrice R
rij 		<- matrix(nrow=1, ncol=NBN-1)	# Calcul intermédiaire des rij
rL 			<- vector(length=NBN)			# Résultat de la méthode de Lasso, NBN-1 * NBN-1 (ie sol. de Yi = Ai * Xi)
MatrCc		<- array(dim=c(NBN, NBN))		# Matrice r calculée selon la méthode RLD ou STEP ou LASSO
vrCc  		<- vector(length = NBN*NBN)		# MatrCc en vecteur
MatrCc1		<- array(dim=c(NBN, NBN))		# Matrice r résultat de l'algorithme
											# Elle est calculée en tronquant les éléments de MatrCc à 0 et 1 (0 sur la diagonale, comme Solution)

MatPVal		<- array(dim=c(NBN, NBN))		# Calcul direct de la p Value

pX   <- array(dim=c(NBN-1, NBN-1, NBN))				# ligne (iRow), colonne (iCol), n° matrice (iMat)
pY   <- array(dim=c(NBN-1, NBN))					# ligne, n° matrice
gX   <- array(dim=c(NBK*(NBN-1), NBN-1, NBN))		# ligne (iRow), colonne (iCol), n° matrice (iMat)
gY   <- array(dim=c(NBK*(NBN-1), NBN))				# ligne, 1, n° matrice
gX1  <- array(dim=c(NBN-1, NBN-1, NBN))				# matrice gX tronquée pour ne prendre que les valeurs KD ou KO
gY1  <- array(dim=c(NBN-1, NBN))					# idem gX1
													
Seuils		<- seq(0, 1, by=0.05)					# Pour le tracé des courbes rOC (KO, KD, RLD)
													# Pour LASSO, on multipliera chaque valeur par KLASSO
NBS			<- length(Seuils)
KLASSO		<- 5									# choix arbitraire, peut être modifié si besoin
ThFig6		<- vector(length=NBM)					# Valeur du seuil, correspondant à la figure 6, pour marquer le point correspondant sur la courbe ROC
dimnames(ThFig6)	<- list(Methods)

Score 		<- array(0, dim=c(NBSC, NBS, NBM))		# Vrais +- (TP), vrais 0 (TN), faux +- (FP), faux 0 (FN)
dimnames(Score)		<- list(ScoreTest, NULL, Methods)

Pts			<- array(0, dim=c(2, NBM))				# Spec et Sensib correspondant aux points à indiquer sur les courbes ROC (méthodes non tracées)
AuROC		<- array(0, dim=c(1, NBM))				# AUROC correspondant aux méthodes tracées
dimnames(Pts)		<- list(c("Sensib", "Specif"), Methods)
dimnames(AuROC)		<- list(NULL, Methods)

#	1/ Lecture des données
	
val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_wildtype.tsv", sep="")
XMesNPLu 	<- fread(val, data.table=F)						# "Documents/DreamChallenge10/insilico_size10_1_wildtype.tsv"
val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_knockdowns.tsv", sep="")
XMesPLu1 	<- fread(val, data.table=F) 					# ("Documents/DreamChallenge10/insilico_size10_1_knockdowns.tsv"
val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_knockouts.tsv", sep="")
XMesPLu2 	<- fread(val, data.table=F) 					# "Documents/DreamChallenge10/insilico_size10_1_knockouts.tsv"

val			<- paste(ParNds$vRoot[indNBN], ParNds$vSolTxt[indNBN], indFIC, ".tsv", sep="")
SolutLu  	<- fread(val, data.table=F) 					# "Documents/DreamChallenge10/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
SolutLu$V1	<- as.factor(SolutLu$V1)
SolutLu$V2	<- as.factor(SolutLu$V2)

for (i in 1:NBN) {
	XMesNP[i,1] = XMesNPLu[1,i]
	XMesNP[i,2] = XMesNPLu[1,i]
	
	for (j in 1:NBN) {
		XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation à -50%	
											# ATTENTION les données de Dream Cha. sont transposées par rapport à mon écriture
		XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation à -100%
	}
}

Solution2[,] = 0
for (i in 1:(NBN*(NBN-1))) {
#	cat("i ", i, " V3 ", SolutLu[i,]$V3, " V1 ", SolutLu[i,]$V1, " : ", as.numeric(SolutLu[i,]$V1), " V2 ", SolutLu[i,]$V2, " : ", as.numeric(SolutLu[i,]$V2), "\n")
	if(SolutLu[i,]$V3 == 1) {
		Solution2[as.numeric(SolutLu[i,]$V1), as.numeric(SolutLu[i,]$V2)] = 1
	} else {
		break								# On suppose que tous les 1 sont au début du fichier  -- Vérifié pour les 10 fichiers le 11/05/22
	}
}

numero <- as.numeric(sort(as.character(c(1:NBN))))

for (i in 1:NBN) {
	Solution3[,numero[i]] = Solution2[,i]
}
for (i in 1:NBN) {
	Solution[numero[i],] = Solution3[i,]	# On retrouve la matrice de connexion classique (noeuds "G1", "G2", ... , "G10")
}
Solution <- t(Solution)						# Deuxième configuration

for (iPaq in c(1:NBK)) {					# On prend toutes les données disponibles : Knock Down et Knock Out
	#	Les "données de mesure", fournies par Dream Challenge, ont été chargées dans XMesNP et XMesP précédemment	
	for (iPert in c(1:NBN)) {
		MatR[,iPert,iPaq]  <- 2 * (XMesP[,iPert,iPaq]-XMesNP[,iPaq]) / (XMesP[,iPert,iPaq]+XMesNP[,iPaq])
	}
	
	for (iMat in 1:NBN) {
		col <- 1:NBN
		col <- col[-iMat]
		
		for (iRow in 1:(NBN-1)) {
			pY[iRow, iMat] = MatR[iMat, col[iRow], iPaq]
			
			for (iCol in 1:NBN-1) {
				pX[iRow, iCol, iMat] = MatR[col[iCol], col[iRow], iPaq]
			}
		}
	}
	
	gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,]  <- pX[, ,]
	gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ]   <- pY[, ]
}
	
for	(iMeth in 1:NBM) {
	Method  <- Methods[iMeth]
	
	#	Recherche de MatrCc	
	
	if (Method %in% c("MRA-KD", "MRA-KO")) {
	#	Méthode MRA classique
		if (Method == "MRA-KD") {
			rSupp <- seq(NBN, NBK*(NBN-1), 1)	# On supprime les lignes 10 à 18 (resp. 100 à 198) pour KD
		} else {
			rSupp <- seq(1, NBN-1, 1)			# On supprime les lignes 1  à 9  (resp. 1 à 99) pour KO
		}
		
		gX1	<- gX[-rSupp,,]
		gY1	<- gY[-rSupp,]
	
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
				
			Form <- paste("gY1[,",iMat,"]~gX1[,,",iMat,"]+0")
			
			MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
			MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
		}
		MatrCc[is.na(MatrCc)]  <- 0
		
#			NB1 = trunc(ParNds$vTOP[indNBN]*NBN*NBN)
#			vrCc[]  <- MatrCc[]
#			vrCc1	<- sort(abs(vrCc), decreasing=TRUE)
#			ThFig6[iMeth]	<- vrCc1[min(NB1+1, NBN*NBN)]	
	}	# Méthodes MRA
	
	if (Method %in% c("RLD", "RLD OPTIM", "RLD OPTIM2")) {
	#	Méthode RLD classique (avec seuils)
	
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
				
			Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")
			
			MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
			MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
		}
		MatrCc[is.na(MatrCc)]  <- 0
		
#			ThFig6[iMeth] <- 0.25*MatrMax
	}	# Méthode RLD
	
	if (Method == "PVAL") {						#	p-Value à 5%
		SeuilPV	  <- 0.05
		MatrCc[,] <- 0
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
				
			Form1 <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")	# Pour comparaison avec la régression linéaire directe		
			pQ 	  <- lm(as.formula(Form1))
			
			MatPVal[iMat,col] 	<- broom::tidy(pQ)$p.value
			MatrCc [iMat,which(MatPVal[iMat,] < SeuilPV)] = 1
		}		# boucle iMat
	}	# Méthode PVAL
	
	if (Method %in% c("STEP-Fo", "STEP-Ba", "STEP-Bo")) {
	#	Méthode "step" (suppression automatique de coefficients par optimisation de AIC : Akaike's Information Criterion)
	
		MatrCc[,] <- 0
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
	
			Donnees  <- data.frame(gY[,iMat], gX[,,iMat])
			#	cat("iMat ", iMat, " Donnees ", colnames(Donnees), "\n")
			Donnees	 <- rename(Donnees, "Y" ="gY...iMat.")		# Le suivi est plus clair comme cela
			cc 		 <- colnames(Donnees)						# Nom des colonnes de "Données". La première est "Y"
			cc		 <- cc[-1]									# Il reste le nom des coefficients ("X1", "X2", ... "Xn"  -- n = NBN-1)
			colnames(rij) <- cc									# step indique le nom des colonnes correspondant aux coefficients conservés
					
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
			
			ll  	<- length(forward$coefficients)				# Nombre de coefficients conservés par step
			nn		<- names(forward$coefficients)				# Nom de ces coefficients
	
			rij[1,]	<- 0
			if (ll >= 2) {
				for (i in 2:ll) {								# nn[1] = "(Intercept)"
					rij[1, nn[i]] <- forward$coefficients[nn[i]]
				}
			}
			if (Verbose) {
				cat (" iMat ", iMat, " ll ", ll, " nn ", nn, " rij ", rij, "\n")
			}
			MatrCc[iMat, col]	<- rij[1,]
			MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
		}		# boucle iMat
	}	# Méthodes STEP
		
	if (Method == "LASSO") {
	#	Méthode "Lasso", choix automatique de Lambda
	
		MatrCc[,] <- 0		# Initialisation
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
	
			cv_model 	<- cv.glmnet(gX[,,iMat], gY[,iMat], alpha = 1)						# Fit lasso regression model using k-fold cross-validation
			best_lambda <- cv_model$lambda.min	
		}		
	}	# Méthode  LASSO
	
	if (Method == "LARS") {
		MatrCc[,] <- 0		# Initialisation
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
	
			pp=lars(gX[,,iMat], gY[,iMat], type="lar")
			rL = coef(pp)
	
			MatrCc[iMat, col]	<- rL[2:NBN]
		}		
	}	# Méthode  LARS
	
	MatrMax = max(abs(MatrCc[,]))
	P = sum(Solution)			# Nb. vrais positifs
	N = NBN*(NBN-1) - P			# Nb. vrais négatifs
	
	#	Numérisation du résultat et calcul du score
	
	for (iSeuil in (1:NBS)) {	
		if (Method %in% c("MRA-KD", "MRA-KO", "RLD")) {
			Th = Seuils[iSeuil]*MatrMax
			
			for (iRow in 1:NBN) {
				for (iCol in 1:NBN) {
					if (abs(MatrCc[iRow, iCol]) > Th) {
						MatrCc1[iRow, iCol] = 1
					} else {
						MatrCc1[iRow, iCol] = 0
					}
				}
				MatrCc1[iRow, iRow] = 0	
			}	# iRow
		}			# Méthodes "MRA-KD", "MRA-KO", "RLD"
	
		else if (Method == "LASSO") {
			Blambda <- Seuils[iSeuil]*KLASSO*best_lambda
			
			MatrCc1[,] <- 0		# Initialisation
			for (iMat in 1:NBN) {
				col <- 1:NBN
				col <- col[-iMat]
				
				best_model 	<- glmnet(gX[,,iMat], gY[,iMat], alpha = 1, lambda = Blambda)
				rL			<- coef(best_model)								
				MatrCc1[iMat, col]	<- rL[2:NBN]
			}	# iMat
		}		# Méthode LASSO
	
		else {
			MatrCc1 = MatrCc
		}
		
		Score[, iSeuil, iMeth] = 0
		for (iRow in 1:NBN) {
			for (iCol in 1:NBN) {
				if (MatrCc1[iRow, iCol] != 0) {
					if (Solution[iRow, iCol] != 0) {
						Score["T+-", iSeuil, iMeth] <- Score["T+-", iSeuil, iMeth] +1
					} else {
						Score["F+-", iSeuil, iMeth] <- Score["F+-", iSeuil, iMeth] +1
					}
				} else {
					if (Solution[iRow, iCol] == 0) {
						Score["T0", iSeuil, iMeth] <- Score["T0", iSeuil, iMeth] +1
					} else {
						Score["F0", iSeuil, iMeth] <- Score["F0", iSeuil, iMeth] +1
					}
				}		
			}
		}	# iRow
		
		Score["Sensib", iSeuil, iMeth] = Score["T+-", iSeuil, iMeth] / P
		Score["Specif", iSeuil, iMeth] = (Score["T0", iSeuil, iMeth] - NBN) / N
	}	# Boucle sur les seuils
	
	if (Method %in% c("MRA-KD", "MRA-KO", "RLD", "LASSO")) {
		AuROC[iMeth]	<- -trapz(1-Score["Specif",,iMeth], Score["Sensib",,iMeth])			# AUC ROC pour les méthodes tracées
	} else {
		Pts["Specif", iMeth] = Score["Specif",1,iMeth]										# Point à dessiner pour les méthodes non tracées
		Pts["Sensib", iMeth] = Score["Sensib",1,iMeth]
	}
	
#		if (Method %in% c("MRA-KD", "MRA-KO", "RLD") {	
#			i = min(which(abs(ThFig6[iMeth]-Seuils) < 1E-12))
#			Pts["Specif", Method] = Score["Specif",i,Method]								# Point correspondant au seuil choisi dans les autres figures
#			Pts["Sensib", Method] = Score["Sensib",i,Method]
#		}
	if (Method == "LASSO") {	
		Pts["Specif", "LASSO"] = Score["Specif",5,"LASSO"]									# + Best Lambda pour la méthode LASSO
		Pts["Sensib", "LASSO"] = Score["Sensib",5,"LASSO"]
	}
}		# Boucle sur les méthodes

# Tracé de la courbe ROC

# AURoc	<- -trapz(1-Score["Specif",], Score["Sensib",])			# AUC ROC			trapz(x,y) les x décroissent

NomFic 	<- paste(ParNds$vRSMTxt[indNBN], indFIC, "_GROC.png", sep="")
png(filename=NomFic)											# Courbe ROC globale
plot (1-Score["Specif",,"MRA-KD"], Score["Sensib",,"MRA-KD"], 
	main=paste("ROC : ", as.character(ParNds$vRSMTxt[indNBN]), indFIC, sep=""), type="l", col="blue", 
	xlab=paste("1 - Specificity\nAU ROC  KD = ", round(AuROC[1,"MRA-KD"], digits=3), " KO = ", round(AuROC[1,"MRA-KO"], digits=3),
			    " TLR = ", round(AuROC[1,"RLD"], digits=3), " LASSO = ", round(AuROC[1,"LASSO"], digits=3)), 
	ylab="Sensibility", xlim=c(0,1), ylim=c(0,1))		# plot(x,y,...)
lines(1-Score["Specif",,"MRA-KO"], Score["Sensib",,"MRA-KO"], col="black")
lines(1-Score["Specif",,"RLD"],    Score["Sensib",,"RLD"],    col="red")
lines(1-Score["Specif",,"LASSO"],  Score["Sensib",,"LASSO"],  col="green")
abline(c(0,0), c(1,1), col="grey", lty="dashed")
points(1-Pts["Specif","PVAL"], 	  Pts["Sensib","PVAL"],    pch=19, col="maroon4")
points(1-Pts["Specif","STEP-Fo"], Pts["Sensib","STEP-Fo"], pch=22, col="maroon4", bg="maroon4")
points(1-Pts["Specif","STEP-Ba"], Pts["Sensib","STEP-Ba"], pch=24, col="maroon4", bg="maroon4")
points(1-Pts["Specif","STEP-Bo"], Pts["Sensib","STEP-Bo"], pch=25, col="maroon4", bg="maroon4")

#	points(1-Pts["Specif","MRA-KD"],  Pts["Sensib","MRA-KD"],  pch=16, col="orange", bg="orange")
#	points(1-Pts["Specif","MRA-KO"],  Pts["Sensib","MRA-KO"],  pch=17, col="orange", bg="orange")
#	points(1-Pts["Specif","RLD"],     Pts["Sensib","RLD"],     pch=15, col="orange", bg="orange")
points(1-Pts["Specif","LASSO"],   Pts["Sensib","LASSO"],   pch=23, col="orange", bg="orange")
legend("bottomright", c("MRA-KD","MRA-KO","TLR","LASSO"), fill=c("blue","black","red","green"))
legend("bottom", c("REG","STEP-Fo","STEP-Ba","STEP-Bo","LASSO"), pch=c(19,22,24,25,23),
	   col=c("maroon4","maroon4","maroon4","maroon4","orange"))
# legend("bottom", c("REG","STEP-Fo","STEP-Ba","STEP-Bo","LASSO", "MRA-KD", "MRA-KO", "TLR"), pch=c(19,22,24,25,23, 16,17,15),
#		   col=c("maroon4","maroon4","maroon4","maroon4","orange","orange","orange","orange"))	   
dev.off()




#
#	Sup. Infos 4
#
#	Tracé des courbes ROC correspondant aux différents "challenges" et aux différentes méthodes.
#	Mettre à jour les variables "indNBN" et "indFIC" pour balayer les différentes tailles (10 ou 100) et les 5 fichiers.
#

Verbose	<- FALSE
indNBN	<- 1				# Valeurs possibles : 1 (pour 10 noeuds), 2 (pour 100 noeuds)
indFIC	<- 1				# Valeurs possibles : 1 à 5
# Pour balayer tous les fichiers, on ne peut pas utiliser une boucle car le stockage des courbes ne fonctionne pas alors.

ScoreTest	<-	c("T+-", "T0", "F+-", "F0", "Sensib", "Specif")	# Scores mesurés
NBSC	<- length(ScoreTest)									# Nb. de scores testés
Methods	<-  c("MRA-KD", "MRA-KO", "RLD")
NBM		<- length(Methods)

vNBN	<- c(10, 100)																			# Nb noeuds
vNBK	<- c(2, 2)
vTOP	<- c(0.2, 0.25)																			# Nb de types de perturbations (KD ou KO)

vRoot	<- c("Documents/DreamChallenge10/", "Documents/DreamChallenge100/")						# Nom de la racine des fichiers
vRSMTxt	<- c("insilico_size10_",  "insilico_size100_")											# Texte à afficher dans les tables de regroupement
vRSMFic	<- c("MRA3_Resume_ZF10",  "MRA3_Resume_ZF100")											# Nom des fichiers de regroupement
vSolTxt	<- c("DREAM4_GoldStandard_InSilico_Size10_", "DREAM4_GoldStandard_InSilico_Size100_")	# Partie du nom des fichiers solution

ParNds	<- data.frame(vNBN, vNBK, vTOP, vRoot, vRSMTxt, vRSMFic, vSolTxt)						# Tableau des données précédentes, correspondant aux nb. de noeuds possibles

NBN = ParNds$vNBN[indNBN]
NBK = ParNds$vNBK[indNBN]					# NBK ne correspond plus à du bruit ajouté, mais à des types de perturbations (KD ou KO)

XMesNP		<- array(dim=c(NBN, NBK))		# Valeur mesurées, système non perturbé
XMesP		<- array(dim=c(NBN, NBN, NBK))	# Valeur mesurées, suite à une perturbation
Solution	<- array(dim=c(NBN, NBN))		# Solution : vraie structure du réseau
Solution2	<- array(dim=c(NBN, NBN))		# Idem "Solution", mais les n° de ligne et colonne correspondent aux n° de "levels" ("G1", "G10", "G2", "G3" etc)
Solution3	<- array(dim=c(NBN, NBN))		# Intermédiaire de calcul

MatR		<- array(dim=c(NBN, NBN, NBK))	# Matrice R
rij 		<- matrix(nrow=1, ncol=NBN-1)	# Calcul intermédiaire des rij
rL 			<- vector(length=NBN)			# Résultat de la méthode de Lasso, NBN-1 * NBN-1 (ie sol. de Yi = Ai * Xi)
MatrCc		<- array(dim=c(NBN, NBN))		# Matrice r calculée selon la méthode RLD ou STEP ou LASSO
vrCc  		<- vector(length = NBN*NBN)		# MatrCc en vecteur
MatrCc1		<- array(dim=c(NBN, NBN))		# Matrice r résultat de l'algorithme
											# Elle est calculée en tronquant les éléments de MatrCc à 0 et 1 (0 sur la diagonale, comme Solution)

MatPVal		<- array(dim=c(NBN, NBN))		# Calcul direct de la p Value

pX   <- array(dim=c(NBN-1, NBN-1, NBN))				# ligne (iRow), colonne (iCol), n° matrice (iMat)
pY   <- array(dim=c(NBN-1, NBN))					# ligne, n° matrice
gX   <- array(dim=c(NBK*(NBN-1), NBN-1, NBN))		# ligne (iRow), colonne (iCol), n° matrice (iMat)
gY   <- array(dim=c(NBK*(NBN-1), NBN))				# ligne, 1, n° matrice
gX1  <- array(dim=c(NBN-1, NBN-1, NBN))				# matrice gX tronquée pour ne prendre que les valeurs KD ou KO
gY1  <- array(dim=c(NBN-1, NBN))					# idem gX1
													
Seuils		<- seq(0, 1, by=0.05)					# Pour le tracé des courbes rOC (KO, KD, RLD)
													# Pour LASSO, on multipliera chaque valeur par KLASSO
NBS			<- length(Seuils)
KLASSO		<- 5									# choix arbitraire, peut être modifié si besoin
ThFig6		<- vector(length=NBM)					# Valeur du seuil, correspondant à la figure 6, pour marquer le point correspondant sur la courbe ROC
dimnames(ThFig6)	<- list(Methods)

Score 		<- array(0, dim=c(NBSC, NBS, NBM))		# Vrais +- (TP), vrais 0 (TN), faux +- (FP), faux 0 (FN)
dimnames(Score)		<- list(ScoreTest, NULL, Methods)

Pts			<- array(0, dim=c(2, NBM))				# Spec et Sensib correspondant aux points à indiquer sur les courbes ROC (méthodes non tracées)
AuROC		<- array(0, dim=c(1, NBM))				# AUROC correspondant aux méthodes tracées
dimnames(Pts)		<- list(c("Sensib", "Specif"), Methods)
dimnames(AuROC)		<- list(NULL, Methods)

#	1/ Lecture des données
	
val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_wildtype.tsv", sep="")
XMesNPLu 	<- fread(val, data.table=F)						# "Documents/DreamChallenge10/insilico_size10_1_wildtype.tsv"
val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_knockdowns.tsv", sep="")
XMesPLu1 	<- fread(val, data.table=F) 					# ("Documents/DreamChallenge10/insilico_size10_1_knockdowns.tsv"
val			<- paste(ParNds$vRoot[indNBN], ParNds$vRSMTxt[indNBN], indFIC, "_knockouts.tsv", sep="")
XMesPLu2 	<- fread(val, data.table=F) 					# "Documents/DreamChallenge10/insilico_size10_1_knockouts.tsv"

val			<- paste(ParNds$vRoot[indNBN], ParNds$vSolTxt[indNBN], indFIC, ".tsv", sep="")
SolutLu  	<- fread(val, data.table=F) 					# "Documents/DreamChallenge10/DREAM4_GoldStandard_InSilico_Size10_1.tsv"
SolutLu$V1	<- as.factor(SolutLu$V1)
SolutLu$V2	<- as.factor(SolutLu$V2)

for (i in 1:NBN) {
	XMesNP[i,1] = XMesNPLu[1,i]
	XMesNP[i,2] = XMesNPLu[1,i]
	
	for (j in 1:NBN) {
		XMesP[i,j,1]  = XMesPLu1[j,i]		# KnockDown (KD)  ie Perturbation à -50%	
											# ATTENTION les données de Dream Cha. sont transposées par rapport à mon écriture
		XMesP[i,j,2]  = XMesPLu2[j,i]		# KnockOut  (KO)  ie Perturbation à -100%
	}
}

Solution2[,] = 0
for (i in 1:(NBN*(NBN-1))) {
#	cat("i ", i, " V3 ", SolutLu[i,]$V3, " V1 ", SolutLu[i,]$V1, " : ", as.numeric(SolutLu[i,]$V1), " V2 ", SolutLu[i,]$V2, " : ", as.numeric(SolutLu[i,]$V2), "\n")
	if(SolutLu[i,]$V3 == 1) {
		Solution2[as.numeric(SolutLu[i,]$V1), as.numeric(SolutLu[i,]$V2)] = 1
	} else {
		break								# On suppose que tous les 1 sont au début du fichier  -- Vérifié pour les 10 fichiers le 11/05/22
	}
}

numero <- as.numeric(sort(as.character(c(1:NBN))))

for (i in 1:NBN) {
	Solution3[,numero[i]] = Solution2[,i]
}
for (i in 1:NBN) {
	Solution[numero[i],] = Solution3[i,]	# On retrouve la matrice de connexion classique (noeuds "G1", "G2", ... , "G10")
}
Solution <- t(Solution)						# Deuxième configuration

for (iPaq in c(1:NBK)) {					# On prend toutes les données disponibles : Knock Down et Knock Out
	#	Les "données de mesure", fournies par Dream Challenge, ont été chargées dans XMesNP et XMesP précédemment	
	for (iPert in c(1:NBN)) {
		MatR[,iPert,iPaq]  <- 2 * (XMesP[,iPert,iPaq]-XMesNP[,iPaq]) / (XMesP[,iPert,iPaq]+XMesNP[,iPaq])
	}
	
	for (iMat in 1:NBN) {
		col <- 1:NBN
		col <- col[-iMat]
		
		for (iRow in 1:(NBN-1)) {
			pY[iRow, iMat] = MatR[iMat, col[iRow], iPaq]
			
			for (iCol in 1:NBN-1) {
				pX[iRow, iCol, iMat] = MatR[col[iCol], col[iRow], iPaq]
			}
		}
	}
	
	gX[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ,]  <- pX[, ,]
	gY[((NBN-1)*(iPaq-1)+1) : ((NBN-1)*iPaq), ]   <- pY[, ]
}
	
for	(iMeth in 1:NBM) {
	Method  <- Methods[iMeth]
	
	#	Recherche de MatrCc	
	
	if (Method %in% c("MRA-KD", "MRA-KO")) {
	#	Méthode MRA classique
		if (Method == "MRA-KD") {
			rSupp <- seq(NBN, NBK*(NBN-1), 1)	# On supprime les lignes 10 à 18 (resp. 100 à 198) pour KD
		} else {
			rSupp <- seq(1, NBN-1, 1)			# On supprime les lignes 1  à 9  (resp. 1 à 99) pour KO
		}
		
		gX1	<- gX[-rSupp,,]
		gY1	<- gY[-rSupp,]
	
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
				
			Form <- paste("gY1[,",iMat,"]~gX1[,,",iMat,"]+0")
			
			MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
			MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
		}
		MatrCc[is.na(MatrCc)]  <- 0
	}	# Méthodes MRA
	
	if (Method %in% c("RLD", "RLD OPTIM", "RLD OPTIM2")) {
	#	Méthode RLD classique (avec seuils)
	
		for (iMat in 1:NBN) {
			col <- 1:NBN
			col <- col[-iMat]
				
			Form <- paste("gY[,",iMat,"]~gX[,,",iMat,"]+0")
			
			MatrCc[iMat, col] 	<- (lm(as.formula(Form)))$coefficients
			MatrCc[iMat, iMat]  <- 0			# L'élément diagonal est mis à 0 au lieu de -1, car la matrice Solution le met à 0.
		}
		MatrCc[is.na(MatrCc)]  <- 0
	}	# Méthode RLD
	
	MatrMax = max(abs(MatrCc[,]))
	P = sum(Solution)			# Nb. vrais positifs
	N = NBN*(NBN-1) - P			# Nb. vrais négatifs
	
	#	Numérisation du résultat et calcul du score
	
	for (iSeuil in (1:NBS)) {	
		Th = Seuils[iSeuil]*MatrMax
		
		for (iRow in 1:NBN) {
			for (iCol in 1:NBN) {
				if (abs(MatrCc[iRow, iCol]) > Th) {
					MatrCc1[iRow, iCol] = 1
				} else {
					MatrCc1[iRow, iCol] = 0
				}
			}
			MatrCc1[iRow, iRow] = 0	
		}	# iRow
		
		Score[, iSeuil, iMeth] = 0
		for (iRow in 1:NBN) {
			for (iCol in 1:NBN) {
				if (MatrCc1[iRow, iCol] != 0) {
					if (Solution[iRow, iCol] != 0) {
						Score["T+-", iSeuil, iMeth] <- Score["T+-", iSeuil, iMeth] +1
					} else {
						Score["F+-", iSeuil, iMeth] <- Score["F+-", iSeuil, iMeth] +1
					}
				} else {
					if (Solution[iRow, iCol] == 0) {
						Score["T0", iSeuil, iMeth] <- Score["T0", iSeuil, iMeth] +1
					} else {
						Score["F0", iSeuil, iMeth] <- Score["F0", iSeuil, iMeth] +1
					}
				}		
			}
		}	# iRow
		
		Score["Sensib", iSeuil, iMeth] = Score["T+-", iSeuil, iMeth] / P
		Score["Specif", iSeuil, iMeth] = (Score["T0", iSeuil, iMeth] - NBN) / N
	}	# Boucle sur les seuils
	
	AuROC[iMeth]	<- -trapz(1-Score["Specif",,iMeth], Score["Sensib",,iMeth])			# AUC ROC pour les méthodes tracées
}		# Boucle sur les méthodes

# Tracé de la courbe ROC

NomFic 	<- paste(ParNds$vRSMTxt[indNBN], indFIC, "_GROC100_1.png", sep="")
png(filename=NomFic)											# Courbe ROC globale
plot (1-Score["Specif",,"MRA-KD"], Score["Sensib",,"MRA-KD"], 
	main=paste("ROC : ", as.character(ParNds$vRSMTxt[indNBN]), indFIC, sep=""), type="l", col="blue", 
	xlab=paste("1 - Specificity\nAU ROC  KD = ", round(AuROC[1,"MRA-KD"], digits=3), " KO = ", round(AuROC[1,"MRA-KO"], digits=3),
			    " TLR = ", round(AuROC[1,"RLD"], digits=3)), 
	ylab="Sensibility", xlim=c(0,1), ylim=c(0,1))		# plot(x,y,...)
lines(1-Score["Specif",,"MRA-KO"], Score["Sensib",,"MRA-KO"], col="black")
lines(1-Score["Specif",,"RLD"],    Score["Sensib",,"RLD"],    col="red")
abline(c(0,0), c(1,1), col="grey", lty="dashed")
legend("bottomright", c("MRA-KD","MRA-KO","TLR"), fill=c("blue","black","red")) 
dev.off()
