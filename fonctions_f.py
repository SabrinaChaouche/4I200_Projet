import numpy as np
import random
from PIL import Image, ImageDraw
from gurobipy import *

# ENTREE:-secteur: La liste de villes reprèsentant les secteurs
# 		:-v: La liste des villes
#		:-x : Les valeurs des variables de décision
#		:-n : Le nombre de villes
# 		:-k : Le nombre de secteurs
# 		:-matnk: La matrice N*k de distances
# RETOURNE : -maxx:La distance maximale entre une ville et un secteur dans la solution x
# 		   : -moy: La moyenne de satisfaction 
#		   : -ville: La ville la moins bien servie
#		   : -sec : secteur où ville est affectée

def maire_moinsbien_servi(n,k,matnk,v,x,secteur):
	maxx=0
	moy=0
	for i in range(n):
		for j in range(k):
			if(x[i*k+j].x==1):
				moy+=matnk[i][j]
				if(maxx<matnk[i][j]):
					maxx=matnk[i][j]
					ville=v[i]
					sec=secteur[j]
	moy=moy/n
	return maxx,moy,ville,sec


# ENTREE:-x : Les valeurs des variables de décision
#		:-n : Le nombre de villes
# 		:-k : Le nombre de secteurs
# Affichage : valeurs des variables de décision

def affichage_valeurs_x(x,n,k):
	print("")				
	print('Solution optimale:')
	for i in range(n):
		for j in range(k):
			print('x%d.%d'%(i+1,j+1), '=', x[i*k+j].x)


# ENTREE:-secteur: La liste de villes reprèsentant les secteurs
# 		:-v: La liste des villes
#		:-x : Les valeurs des variables de décision
#		:-n : Le nombre de villes
# 		:-k : Le nombre de secteurs
# 		:-cv: La liste des coordonnées des villes
#		:-s : Le nombre de secteurs (pour une utilisation ultérieure dans g amélioré)
#		:-im: L'image sur laquelle nous allons afficher la solution
# 		:-alpha: paramètre donné
# AFFICHAGE : Image de la carte géographique des hauts-de-Seine modélisant la solution trouvée
def affichage_image(im,n,k,x,secteur,cv,v,alpha,s):
	image=Image.open(im)
	draw=ImageDraw.Draw(image)
	couleur=['black','red','navy','purple','green']
	courant=[]
	for j in range(k):
		l=random.randint(0,4)
		color=couleur[l]
		courant.append(l)
		for i in range(n):
			if(x[i*k+j].x==1):
				point=cv[i][0],cv[i][1]
				draw.text((float(cv[i][0]),float(cv[i][1])),v[i], fill='grey')
				draw.line([(float(cv[i][0]),float(cv[i][1])),(float(cv[v.index(secteur[j])][0]), float(cv[v.index(secteur[j])][1]))],fill=color, width=2	)   
	image.save("/home/sou/master1/mogpl/cartes/"+str(s)+"alpha"+str(int(10*alpha))+".png")
	image.show()


# ENTREE:-secteur: La liste de villes reprèsentant les secteurs
# 		:-v: La liste des villes
#		:-n : Le nombre de villes
# 		:-k : Le nombre de secteurs
# RETOURNE : La matrice des contraintes de la fonction f
def matrice_contrainte(n,k,p,secteur,v):
	nbcont=n+2*k+1
	nbvar=n*k
	a = np.zeros((nbcont,nbvar))
	
	#contrainte:chaque ville appartient a un seul et unique secteur
	for i in range(n):
		for j in range(k):
			a[i][i*k+j]=1

	#contrainte :le nombre maximal de population par secteur 
	for i in range(n,k+n):
		for j in range(k*n):
			if((j%k==0)):
				a[i][(i-n)+j]=p[int(j/k)]	
	
	#contrainte: chaque secteur représente son propre secteur 

	for j in range(n+k,2*k+n):
		for i in range(n):
			if v.index(secteur[j-(n+k)])==i:
				a[j][i*k+(j-(n+k))]=1  
	
	#contrainte: toutes les villes affectés a des secteurs
	for i in range(n*k):
		a[2*k+n][i]=1
	return a


# ENTREE:-n : Le nombre de villes
# 		:-k : Le nombre de secteurs
# 		:-matnk: La matrice N*k de distances
# RETOURNE : Le vecteur des coefficients de la fonction objectif 
def fonction_obj(matnk,n,k):
	c=[]
	for i in range(n):
		for j in range(k):
			c.append(matnk[i][j])
	return c
	
	
# ENTREE:-n : Le nombre de villes
# 		:-k : Le nombre de secteurs
# 		:-sigma: Le nombre maximal d'individus qu'un secteur peut accueillir
# RETOURNE : Le vecteur du second membre de la fonction f

def second_membre(n,k,sigma):
	b=[]
	for i in range(n):
		b.append(1)
	for i in range(k):
		b.append(sigma)
	for i in range(k):
		b.append(1)
	b.append(n)
	return b

	
# ENTREE:-nbcont:Le nombre de contraintes
#		:-nbvar:Le nombre de variables
#		:-n : Le nombre de villes
# 		:-k : Le nombre de secteurs
#		:-a : La matrice des contraintes de la fonction f
#		:-b : Le vecteur du second membre de la fonction f
#		:-c : Le vecteur de la fonction objectif
# RETOURNE: Le modèle construit et les valeurs des variables de décision

def modele(n,k,a,b,c,nbcont,nbvar):
	

	# Range of plants and warehouses
	lignes = range(nbcont)
	colonnes = range(nbvar)

	m = Model("mogplex")	 
			
	# declaration variables de decision
	x = []
	for i in range(n):
		for j in range(k):
			x.append(m.addVar(vtype=GRB.BINARY, lb=0, name="x%d.%d" % (i+1, j+1)))
			
	# maj du modele pour integrer les nouvelles variables
	m.update()
	
	obj = LinExpr();
	obj = 0
	for j in colonnes:
		obj += c[j] * x[j]
			
	# definition de l'objectif
	m.setObjective(obj,GRB.MINIMIZE)
	
	########################## Definition des contraintes ######################
	for i in range(n):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) == b[i], "Contrainte%d" % i)
	for i in range(n,k+n):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) <= b[i], "Contrainte%d" % i)
	for i in range(n+k,2*k+n):
		m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) == b[i], "Contrainte%d" % i)
	
	m.addConstr(quicksum(a[2*k+n][j]*x[j] for j in colonnes) >= b[2*k+n], "Contrainte%d" % int(2*k+n))
	
	# Resolution
	m.optimize()

	return m,x

