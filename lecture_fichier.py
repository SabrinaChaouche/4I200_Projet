import numpy as np
from PIL import Image, ImageDraw
from gurobipy import *
from random import randint

############### lecture fichier distance ##############

def lecture_fichier_distance(dist):
	distance = open(dist, "r")
	dis=distance.read()
	distance.close()
	l=[]
	for ligne in dis.split('\n'):
		l.append(ligne)
	return l


############## lecture fichier villes ################

def lecture_fichier_ville(ville):
	villes = open(ville, "r")
	ville=villes.read()
	villes.close()   
	v=[]
	for ligne in ville.split('\n'):
		if ligne:
			v.append(ligne)
	return v


############# lecture fichier population #################

def lecture_fichier_population(pop):
	population = open(pop, "r")
	popu=population.read()
	population.close()   
	p=[]
	for ligne in popu.split('\n'):
		for li in ligne.split(','):
			if(li.isdigit()):
				p.append(int(li))
	return p


############ lecture du fichier coordonnées des villes ################

def lecture_fichier_coorville(coordonnees):
	coordvilles=open("coordvilles92.txt","r")
	coordv=coordvilles.read()
	coordvilles.close()
	cv=[]
	for ligne in coordv.split('\n'):
		if(ligne):
			lala=[]
			for li in ligne.split(','):
				if(li.isdigit()):
					lala.append(li)
			cv.append(lala)
	return cv



################## Matrice N*N distances ##################
# ENTREE:-l: contenu du fichier distance

# RETOURNE :Matrice N*N de distances
def matrice_distances(l,n): 
	matd=np.zeros((n,n))    
	for i in range(n):
		for j in range(n+1):
			if(j != 0):
				matd[i][j-1]=l[i*(n+1)+j]
	return matd


# ENTREE:-k : le nombre de secteurs
# 		:-v : la liste des villes

# RETOURNE : une liste de secteurs tirés aléatoirement dans v
def secteurs_aleatoire(k,v):
	secteurs_alea=[]
	for i in range(k):
		l=randint(0,len(v)-1)
		secteurs_alea.append(v[l])
	return secteurs_alea
