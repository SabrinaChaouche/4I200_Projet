from gurobipy import *
import numpy as np
import lecture_fichier as lf
import fonctions_f as ff
import fonctions_g as fg
import fonctions_gam as ga
from PIL import Image, ImageDraw

############ nom des fichiers ###################
dist="distances92.txt"
ville="villes92.txt"
pop="populations92.txt"
coordonnees="coordvilles92.txt"
im="92.png"

########### lecture des fichiers ################
l=lf.lecture_fichier_distance(dist)
v=lf.lecture_fichier_ville(ville)
p=lf.lecture_fichier_population(pop)
cv=lf.lecture_fichier_coorville(coordonnees)


########### initialisation des données ##########
#secteur=[v[1],v[26],v[33]]
#secteur=['Garches','Le Plessis-Robinson','Sevres']
#secteur=lf.secteurs_aleatoire(3,v)
secteur=['Sevres', 'Antony', 'Suresnes','Le Plessis-Robinson','Garches']
#secteur=[v[4],v[5],v[19]]
alpha=0.2
n=len(v)
k=len(secteur)
sigma=((1+alpha)/k)*sum(p)

############# Matrice N*N distances ##########
matd= lf.matrice_distances(l,n)

############# matrice N*K distances ########## 

matnk=np.column_stack(([matd[:, v.index(secteur[j])] for j in range(k)] ))



# ENTREE:-secteur: La liste de villes reprèsentant les secteurs
# 		:-alpha: Paramètre donné
# AFFICHAGE : Image de la carte géographique des hauts-de-Seine modélisant la solution trouvée
# RETOURNE : -Les valeurs des variables de décision
# 		   : -Le modèle construit 
def solution_optimale_f(secteur,alpha):
	nbcont=n+2*k+1
	nbvar=n*k

	########### Matrice des contraintes ,##########
	a=ff.matrice_contrainte(n,k,p,secteur,v)

	########### Coefficients de la fonction objectif ########
	c=ff.fonction_obj(matnk,n,k)
	
	################ Second membre #################
	b=ff.second_membre(n,k,sigma)
	
	############ variables de decitions ###########
	m,x=ff.modele(n,k,a,b,c,nbcont,nbvar)

	# maire le moins bien servi
	maxx,moy,ville,sec=ff.maire_moinsbien_servi(n,k,matnk,v,x,secteur)

	############# affichage du resultat ####################

	#ff.affichage_valeurs_x(x,n,k)
	print('__________')
	print('Fonction f')
	print('alpha = ',alpha)
	print('Valeur de la fonction objectif :', m.objVal)	
	print('les secteurs choisis sont : ', [secteur[i] for i in range(len(secteur))])	
	print('le maire le moins bien servi est le maire de la ville de', ville,'avec une distance egale a',maxx  )
	print('la satisfaction moyenne des maires est egale a', moy )
	print('__________')
	
	################ affichage de l'image #####################

	ff.affichage_image(im,n,k,x,secteur,cv,v,alpha,k)
	m.update()
	return m,x
	




# ENTREE:-secteur: La liste de villes reprèsentant les secteurs
# 		:-alpha: Paramètre donné
#		:-epsilon :  Paramètre donné = 10 puissance(-6)
# AFFICHAGE : Image de la carte géographique du département des Hauts-de-Seine modélisant la solution trouvée
# RETOURNE : -Les valeurs des variables de decision
# 		   : -Le modèle construit
def solution_optimale_g(secteur,alpha,epsilon):
	nbcont1=2*n+3*k+1
	nbvar1=n*k+1

	# Range of plants and warehouses
	lignes = range(nbcont1)
	colonnes = range(nbvar1)

	############ Coefficients de la fonction objectif ##############
	c=fg.fonction_obj_g(matnk,n,k,epsilon)


	############## matrice des contraintes g ################
	a=fg.matrice_contraintes_g(n,k,p,matnk,secteur,v,nbcont1,nbvar1)

	##################### Second membre ##################

	b=fg.second_membre_g(n,k,sigma)
			
	#modele 
	m,x=fg.modele(n,k,a,b,c,nbcont1,nbvar1)

	# maire le moins bien servi
	maxx,moy,ville,sec=ff.maire_moinsbien_servi(n,k,matnk,v,x,secteur)

	############ Affichage valeur de X et Z ################

	#fg.affichage_valeurs_x_g(x,n,k)
	print('____________')
	print('Fonction g')
	print('alpha = ',alpha)
	print('les secteurs choisis sont : ', [secteur[i] for i in range(len(secteur))])	
	print('le maire le moins bien servi est le maire de la ville de', ville,'avec une distance du secteur',sec,' egale a',maxx  )
	print('la satisfaction moyenne des maires est egale a', moy )
	print('____________')

	#### affichage de l'image
	ff.affichage_image(im,n,k,x,secteur,cv,v,alpha,k)
	m.update()
	return m,x




# ENTREE:-alpha: Paramètre donné
#		:-epsilon :  Paramètre donné= 10 puissance(-6)
# AFFICHAGE : Image de la carte géographique des Hauts-de-Seine modélisant la solution trouvée
# RETOURNE : -Les valeurs des variables de décision
# 		   : -Le modèle construit
def solution_optimale_g_am(alpha,epsilon,k):
	nbcont2=5*n+2
	nbvar2=n*n+n+1
	#n*n reprèsente xij iE{0....n} jE{0,....,n}
	#n represente yi iE{0,....,n}	y[i]=1 si i est un secteur =0 sinon

	# Range of plants and warehouses
	lignes = range(nbcont2)
	colonnes = range(nbvar2)

	############ Coefficients de la fonction objectif ##############
	c=ga.fonction_obj_gam(matd,n,epsilon)

	############## matrice des contraintes g amélioré (question3)################
	a=ga.matrice_contrainte_gam(n,p,v,matd,nbcont2,nbvar2)

	##################### Second membre ##################
	b=ga.second_membre_gam(n,k,sigma)

	#modele 
	m,x=ga.modele_gam(n,a,b,c,nbcont2,nbvar2)	

	# maire le moins bien servi
	maxx,moy,ville,sec=ff.maire_moinsbien_servi(n,n,matd,v,x,v)
	
	############ Affichage valeur de X,Z et Y ################


	#ga.affichage_valeurs_x_gam(x,n)
	print('____________')
	print('Fonction g amélioré')
	print('alpha = ',alpha)
	print('les secteurs choisis sont : ', ga.secteurs_choisis(x,n,v))	
	print('le maire le moins bien servi est le maire de la ville de', ville,'avec une distance du secteur',sec,' egale a',maxx  )
	print('la satisfaction moyenne des maires est egale a', moy )
	print('____________')
	
	########### Affichage de l'image ##############
	ff.affichage_image(im,n,n,x,v,cv,v,alpha,k)
	m.update()
	return x



# ENTREE:-xf:les valeurs des variables de décision calculées par la fonction f
# 		:-xg:les valeurs des variables de décision calculées par la fonction g
#		:-n : le nombre de villes
# 		:-k : le nombre de secteurs
# 		:-matnk: la matrice N*N de distances
# RETOURNE : le prix de l'équité PE
def prix_equite(xf,xg,n,k,matnk):
	fxf=0
	fxg=0
	for i in range(n):
		for j in range(k):
			fxf+=matnk[i][j]*xf[i*k+j].x
			fxg+=matnk[i][j]*xg[i*k+j].x
	return 1-fxf/fxg




############# test de la fonction f:
#mf,xf=solution_optimale_f(secteur,alpha)


############# test de la fonction g:
#mg,xg=solution_optimale_g(secteur,alpha,10**(-6))



############# test de la fonction g amélioré :
#(alpha,10**(-6),5)


############# prix de l'equité ############
#print(prix_equite(xf,xg,n,k,matnk))























