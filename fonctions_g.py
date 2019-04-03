from gurobipy import *
import numpy as np
import fonctions_f as ff
from PIL import Image, ImageDraw


# ENTREE:-x : les valeurs des variables de décision
#		:-n : le nombre de villes
# 		:-k : le nombre de secteurs
# Affichage : valeurs des variables de décision 

def affichage_valeurs_x_g(x,n,k):
	ff.affichage_valeurs_x(x,n,k)
	print('z=',x[k+1].x)
	print("")


# ENTREE:-nbcont1:le nombre de contraintes
#		:-nbvar1:le nombre de variables
#		:-n : le nombre de villes
# 		:-k : le nombre de secteurs
#		:-a : la matrice des contraintes
#		:-b : le vecteur du second membre
#		:-c : le vecteur de la fonction objectif
#RETOURNE: - le model construit
#		 : - les valeurs des variables de decision

def modele(n,k,a,b,c,nbcont1,nbvar1):
	lignes = range(nbcont1)
	colonnes = range(nbvar1)
	m = Model("mogplex")	 
			
	# declaration variables de decision
	x = []
	for i in range(n):
		for j in range(k):
			x.append(m.addVar(vtype=GRB.BINARY, lb=0, name="x%d.%d" % (i+1, j+1)))
	x.append(m.addVar(vtype=GRB.CONTINUOUS, lb=0, name="z" ))
	   
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
	for i in range(n,k+n):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) <= b[i], "Contrainte%d" % i)
	for i in range(n+k,2*k+n):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) >= b[i], "Contrainte%d" % i)
	
	m.addConstr(quicksum(a[2*k+n][j]*x[j] for j in colonnes) >= b[2*k+n], "Contrainte%d" % int(2*k+n))
	for i in range(2*k+n+1,2*n+2*k+1):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) <= b[i], "Contrainte%d" % i)
	for i in range(2*n+2*k+1,2*n+3*k+1):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) == b[i], "Contrainte%d" % i)
	
	# Resolution
	m.optimize()
	return m,x


# ENTREE:-n : le nombre de villes
# 		:-k : le nombre de secteurs
# 		:-sigma: le nombre maximal d'individus qu'un secteur peut accueillir
# RETOURNE : Le vecteur du second membre de la fonction g

def second_membre_g(n,k,sigma):
	b=ff.second_membre(n,k,sigma)
	for i in range(n):
		b.append(0)
	for i in range(k):
		b.append(1)	
	return b


# ENTREE:-secteur: la liste des villes reprèsentant les secteurs
# 		:-v : la liste des villes
#		:-n : le nombre de villes
# 		:-k : le nombre de secteurs
# 		:-matnk: la matrice N*k de distances
# 		:-p : la liste de la population des villes
#		:-nbcont1: le nombre de contraintes
#		:-nbvar1: le nombre de variables


# RETOURNE : La matrice des contraintes

def matrice_contraintes_g(n,k,p,matnk,secteur,v,nbcont1,nbvar1):

	a = np.zeros((nbcont1,nbvar1))
	#### contrainte :chaque ville appartient a un seul et unique secteur

	for i in range(n):
		for j in range(k):
			a[i][i*k+j]=1
		a[i][k*n]=0


	#### contrainte : le nombre maximal de population par secteur 

	for i in range(n,k+n):
		for j in range(k*n):
			if((j%k==0)):
				a[i][(i-n)+j]=p[int(j/k)]	
		a[i][k*n]=0
	
	#contrainte sur ...
	for i in range(n+k,2*k+n):
		for j in range(k*n):
			if((j%k==0)):
				a[i][(i-(n+k))+j]=1
		a[i][k*n]=0
	
	#### contrainte: toutes les villes affectés a des secteurs
	for i in range(n*k):
		a[2*k+n][i]=1
	a[2*k+n][k*n]=0

	
	#### contrainte sur Z
	for i in range(2*k+n+1,2*n+2*k+1):
		for j in range(k):
			a[i][(i-(2*k+n)-1)*k+j]=matnk[i-(2*k+n)-1][j]
		a[i][n*k]=-1


	#### contrainte: chaque secteur représente son propre secteur 

	for i in range(2*n+2*k+1,2*n+3*k+1):
		for l in range(n):
			if v.index(secteur[i-(2*n+2*k+1)])==l:
				a[i][l*k+(i-(2*n+2*k+1))]=1
	return a 


# ENTREE:-n : le nombre de villes
# 		:-k : le nombre de secteurs
# 		:-matnk: la matrice N*k de distances
# RETOURNE : Le vecteur de la fonction objectif

def fonction_obj_g(matnk,n,k,epsilon):
	c=[]
	for i in range(n):
		for j in range(k):
			c.append(epsilon*matnk[i][j])
	c.append(1)
	return c
