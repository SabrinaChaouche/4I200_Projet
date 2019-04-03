from gurobipy import *
import numpy as np
import fonctions_f as ff
from PIL import Image, ImageDraw


########### fonctions de g Amelioré (Question3) ################

# ENTREE:-x : les valeurs des variables de décision
#		:-n : le nombre de villes
# Affichage : valeurs des variables de décision 

def affichage_valeurs_x_gam(x,n):
	ff.affichage_valeurs_x(x,n,n)
	print('z=',x[n*n+1].x)
	for i in range(n):
		print('y%d'%(i+1), '=', x[n*n+1+i].x)
	print("")



# ENTREE:-nbcont2: le nombre de contraintes
#		:-nbvar2: le nombre de variables
#		:-n : le nombre de villes
#		:-a : la matrice des contraintes
#		:-b : le vecteur du second membre
#		:-c : le vecteur de la fonction objectif
#RETOURNE: - le model construit
#		 : - les valeurs des variables de decision

def modele_gam(n,a,b,c,nbcont2,nbvar2):
	lignes = range(nbcont2)
	colonnes = range(nbvar2)
	
	#m,x=ga.modele_gam(n,a,b,c,nbcont2,nbvar2)	
	m = Model("mogplex")	 
			
	# declaration variables de decision
	x = []
	for i in range(n):
		for j in range(n):
			x.append(m.addVar(vtype=GRB.BINARY, lb=0, name="x%d.%d" % (i+1, j+1)))
	x.append(m.addVar(vtype=GRB.CONTINUOUS, lb=0, name="z" ))
	for i in range(n):
		x.append(m.addVar(vtype=GRB.BINARY, lb=0, name="y%d" % (i+1)))
	
	# maj du modele pour integrer les nouvelles variables
	m.update()


	obj = LinExpr();
	obj = 0
	for j in colonnes:
		obj += c[j] * x[j]
	
	# definition de l'objectif

	m.setObjective(obj,GRB.MINIMIZE)
	
	# Definition des contraintes
	for i in range(n):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) == b[i], "Contrainte%d" % i)
	for i in range(n,2*n):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) <= b[i], "Contrainte%d" % i)
	for i in range(2*n,3*n):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) == b[i], "Contrainte%d" % i)
	for i in range(3*n,3*n+1):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) == b[i], "Contrainte%d" % i)
	for i in range(3*n+1,4*n+1):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) <= b[i], "Contrainte%d" % i)
	for i in range(4*n+1,5*n+1):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) <= b[i], "Contrainte%d" % i)
	
	for i in range(5*n+1,5*n+2):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) == b[i], "Contrainte%d" % i)
	
	
	# Resolution
	m.optimize()
	
	return m,x

# ENTREE:-n : le nombre de villes
# 		:-k : le nombre de secteurs
# 		:-sigma: le nombre maximal d'individus qu'un secteur peut accueillir
# RETOURNE : Le vecteur du second membre de la fonction g amélioré

def second_membre_gam(n,k,sigma):
	b=ff.second_membre(n,n,sigma)
	for i in range(2*n,3*n):
		b[i]=0
	for i in range(n):
		b.append(0)
	for i in range(n):
		b.append(0)	
	b.append(k)

	return b


# ENTREE:-v: la liste des villes
#		:-n : le nombre de villes
# 		:-matd: la matrice N*N de distances
# 		:-p : la liste de la population des villes
#		:-nbcont2: le nombre de contraintes
#		:-nbvar2: le nombre de variables


# RETOURNE : La matrice des contraintes de g amélioré

def matrice_contrainte_gam(n,p,v,matd,nbcont2,nbvar2):

	a = np.zeros((nbcont2,nbvar2))
	
	#### contrainte :chaque ville appartient a un seul et unique secteur
	for i in range(n):
		for j in range(n):
			a[i][i*n+j]=1
		a[i][n*n]=0

	#### contrainte : le nombre maximal de population par secteur 
	for i in range(n,2*n):
		for j in range(n*n):
			if((j%n==0)):
				a[i][(i-n)+j]=p[int(j/n)]	
		a[i][n*n]=0

	### contrainte: chaque secteur représente son propre secteur 
	for i in range(2*n,3*n):
		for j in range(n*n):
			if ((i-(n+n))*n)==j:
				a[i][j+(i-2*n)]=1 
		a[i][n*n+(i-(2*n)+1)]=-1

	### contrainte: toutes les villes affectés a des secteurs
	for i in range(n*n):
		a[3*n][i]=1

	
	### contrainte sur Z
	for i in range(3*n+1,4*n+1):
		for j in range(n):
			a[i][(i-(3*n+1))*n+j]=matd[i-(3*n+1)][j]
		a[i][n*n]=-1

  	### contrainte : le nombre de ville afféctés a un secteur ne depasse pas n
	#			 si une ville n'est pas un secteur alors elle contient 0 villes 
	
	for i in range(4*n+1,5*n+1):
		for j in range(n*n):
			if (j%(n)==0):
				a[i][j+(i-(4*n+1))]=1

	for i in range(4*n+1,5*n+1):
		for j in range(n*n+1,n*n+n+1):
			if(j%(n)==0):
				a[i][j+(i-(5*n))]=(-n)

	### contrainte: le nombre de secteur total est egal a k

	for i in range(n*n+1,n*n+n+1):
		a[5*n+1][i]=1

	return a


# ENTREE:-n : le nombre de villes
# 		:-epsilon : paramètre donné = 10 puissance(-6)
# 		:-matnk: la matrice N*k de distances
# RETOURNE : Le vecteur de la fonction objectif de g amelioré

def fonction_obj_gam(matd,n,epsilon):
	c=[]
	for i in range(n):
		for j in range(n):
			c.append(epsilon*matd[i][j])
	c.append(1)
	for i in range(n):
		c.append(0)
	return c


# ENTREE:-n : le nombre de villes
# 		:-v : la liste des villes
# 		:-x : les valeurs des variables de décision

# RETOURNE : une liste de secteurs choisis par la fonction g amélioré

def secteurs_choisis(x,n,v):
	secteurs=[]
	for i in range(n):
		if x[n*n+1+i].x==1:
			secteurs.append(v[i])
	return secteurs

	
