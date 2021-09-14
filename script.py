
"""@package structure
Secondary structures assignment
    1. helices (Helice)
    2. feuillet β parallèle  (B para)
    3. feuillet β anti-parallèle (B anti)
    4. n-turns (Turn)
    5. coil (Coil)
"""

__author__ = "Rita Eid"

import warnings
warnings.filterwarnings("ignore") 
#pour ne pas imprimer les warnings reçus lors de l'utilisation de certains pdb 

# Biopython :
from Bio.PDB import PDBParser
import sys, os
import numpy as np
import datetime as dt
import shutil

try :
	pdb = sys.argv[1] #c'est la pdb en question qui est donnée par l'utilisateur
	chainechoisi= sys.argv[2] # c'est la chain en question dont l'utilisateur veut avoir la structure secondaire
except:
	print ("votre entrée n'est pas la bonne, je vous pris de revoir le README pour pouvoir poursuivre ")
	sys.exit() # si l'utilisateur ne donne pas le bon input il reçoit ce message

pro_type= pdb [-3:]
if (pro_type != "pdb"):
	print("le fichier que vous avait mis en entrée n'est pas un fichier pdb je vous pris de revoir le README pour pouvoir poursuivre ")
	sys.exit() # le fichier de la proteine à structurer doit etre un pdb

parser = PDBParser()
prot_id = pdb [5:9] # c'est le nom du fichier pdb
prot_file = pdb #c'est tout le non du fichier avec .pdb
structure = parser.get_structure(prot_id, prot_file) # on parse la proteine
model=structure[0]
chain=model[chainechoisi]# c'est la chain choisie
res_start= 1
Hbond=[]#liste des liaisons hydrogénes
atome_name=[]#liste du nom des atome
with open(pdb, "r") as pdb_file:
	for line in pdb_file:
		if (line.startswith("ATOM") and line[13:15] =="H "):
			atome_name.append(line[13:15])	
#on s'assure que le fichier pdb contient des hydrogénes soit l'utilisateur doit avoir bien ajouter des hydrogénes dans son fichier pdb grace à reduce (l'étape decrite dans le README) 		
if (len(atome_name)==0):
	print("votre fichier pdb input ne contient pas des Hydrogénes donc cet algorithme ne peut pas etre calculer, je vous pris de revoir le README pour pouvoir poursuivre surtout la partie de l'hydratation")
	sys.exit()
# si le fichier ne contient pas des hydrogenes le message qui resort est le suivant	
	
for i in range (0,len(chain)):
	for j in range(i+1,len(chain)):
		try:
			idx=res_start+i
			idj=res_start+j
			res=chain[idj]
			if (res.resname== "PRO"): # la proline ne possède pas d'hydrogéne sur sa chaine principales et ne peut pas former de liaison hydrogéne.
				continue
			disON = np.linalg.norm(chain[idj]["N"].coord-chain[idx]["O"].coord)#distance ON
			disCH = np.linalg.norm(chain[idj]["H"].coord-chain[idx]["C"].coord)#distance CH
			disOH = np.linalg.norm(chain[idj]["H"].coord-chain[idx]["O"].coord)#distance OH
			disCN = np.linalg.norm(chain[idj]["N"].coord-chain[idx]["C"].coord)#distance CN
			q1=0.41
			q2=0.20
			F=332
			E=q1*q2*(1/disON + 1/disCH - 1/disOH - 1/disCN)*F
#calcule de l'énergie pour determiner s'il y a une liaison hydrogéne entre les residuts ou non
			if (E < -0.5):
				Hbond.append([E,i,j])
#liste de tout les liaisons hydrogènes possibles dans la protéines
		except KeyError:
			pass
 
Hbond_clean=[]#liste de liaison nettoyer
for n in range (0,len(Hbond)):
	if not (Hbond[n][2] ==Hbond[n][1]+1 or Hbond[n][1] ==Hbond[n][2]+1 or Hbond[n][2] ==Hbond[n][1]+2 or Hbond[n][1] ==Hbond[n][2]+2 or Hbond[n][2] ==Hbond[n][1]): 
#si la liste de liaison en indique une entre 2 acide animés qui se suivent on ne les prend pas en compte
		Hbond_clean.append(Hbond[n])
# on nettoie la liste de liaisons hydrogénes 

resultat=[]#liste des resultats à compléter à fur et à mesure
for idx in range(1,len(chain)+1):
	try:
		res=chain[idx]
		st="Coil"
		resultat.append([res.resname, idx, st])
# list des acide animés , leur mon, indice et strucutre secondaire
	except KeyError:
		continue
#ici on cré une liste qui contient tout les acide animés qui sont present dans la pdb. On assigne la structure secondaire C qui est un coil  pour tout les acide amimés cela sera changer avec les conditions des differentes structures secondaires.

indx1=list(range(0,len(chain)+1))
for t in range(0,len(Hbond_clean)) :
	if (Hbond_clean[t][2]==Hbond_clean[t][1]+3 or Hbond_clean[t][2]==Hbond_clean[t][1]+4 or Hbond_clean[t][2]==Hbond_clean[t][1]+5 ) : 
#condition pour les turns: S'il y a une liaison hydrogéne entre i+3, i+4 ou i+5 	
		aa1=int(Hbond_clean[t][1])# 1ere acide animés de la liaison
		aa2=int(Hbond_clean[t][2])# 2eme acide animés de la liaison
		for idx2 in range(len(resultat)):
			if (resultat[idx2][1]== aa1 or resultat[idx2][1]== aa2):
				resultat[idx2][2]="Turn" # on change la structure secondaire des residus en question dans la liste resultat en T qui et un Turn 


for idx3 in range(0,len(resultat)):
	if (resultat[idx3][2]=='Turn' and resultat[idx3+1][2]=='Turn'):
		resultat[idx3][2]="Helice"
for idx4 in range(0,len(resultat)):
	if (resultat[idx4][2]=='Turn' and resultat[idx4-1][2]=='Helice'):
		resultat[idx4][2]="Helice"
for idx5 in range(1,len(resultat)-1):
	if (resultat[idx5-1][2]=='Helice' and resultat[idx5+1][2]=='Helice'):
		resultat[idx5][2]="Helice"
#les differents condition pour les helices et le combrage des structures secondaire                 

#condition des Brins
diff=[]#list des différences entre des residus d'une liaison hydrogéne
dif=0
for f in range(len(Hbond_clean)):
	br1=int(Hbond_clean[f][1])
	br2=int(Hbond_clean[f][2])
	dif=br1-br2#calcule de la différence entre des residus d'une liaison hydrogéne
	diff.append([dif,br1,br2])
#nous avons une liste avec les différences entre des residus

for d in range(len(diff)):
	if(diff[d][0]<-6):#si la différence est inférieure à -6 soit inferieur q'une helice
		if (diff[d][0]==diff[d-1][0]):# si 2 différences sont égales soit entre 2 liaison hydrogénes avec une égale différence. c'est des brin Beta parallèle.
			brini1=int(diff[d][1])
			brini2=int(diff[d][2])
			for idx7 in range(len(resultat)):
				if (resultat[idx7][1]== brini1 or resultat[idx7][1]== brini2):
					resultat[idx7][2]="B para"

		elif (diff[d-1][2]-diff[d][2]< 4 ):#c'est la différence entre 2 résidus des liaisons qui se suivent. les 2 liaison sont asset proche pour appatenir au meme feulliet beta.
			brini3=int(diff[d][1])
			brini4=int(diff[d][2])
			for idx8 in range(len(resultat)):
				if (resultat[idx8][1]== brini3 or resultat[idx8][1]== brini4 ):
					resultat[idx8][2]="B anti"
		elif (diff[d][0]>-8): 
			brini7=int(diff[d][1])
			brini8=int(diff[d][2])
			for brin2 in range(brini7,brini8):
				for idx9 in range(len(resultat)):
                    			if (resultat[idx9][1]== brin2):
                        			if not (resultat[idx9][2]== 'Turn'):
                            				resultat[idx9][2]="B anti"
#remplissage des residus dont la structure secondaire est connue mais non retrouvée dans la liste des liaison hydrogene
for idx10 in range(len(resultat)-1):
    if (resultat[idx10-1][2]== 'B para' and resultat[idx10+1][2]== 'B para'):
        resultat[idx10][2]="B para"
for idx11 in range(len(resultat)-1):
    if (resultat[idx11-1][2]== 'B anti' and resultat[idx11+1][2]== 'Coil'and resultat[idx11][2] =="Helice"):
        resultat[idx11][2]=resultat[idx11-1][2]
for idx12 in range(len(resultat)-1):
    if (resultat[idx12-1][2]== 'B anti' and resultat[idx12+1][2]== 'B anti'):
        resultat[idx12][2]="B anti"



#sortie des resultats
header = "==== Secondary Structure Assignment using DSSP method ====\nDATE\t\t{}\n".format(dt.date.today())
header += "REFERENCE\tW. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637\n"
header += "PROTEINE\t{}\n".format(prot_id)
header += "CHAINE\t\t{}\n".format(chainechoisi)
print (header+"\n")
assignement="Secondary structures assignment\n 1.helices (Helice) \n 2. feuillet β parallèle  (B para) \n 3. feuillet β anti-parallèle (B anti) \n 4. n-turns (Turn) \n 5. coil (Coil) \n "
print(assignement)
print ("RESIDUE AA STRUCTURE")
for l in range(0,len(resultat)):
	print('{}\t{}\t{}'.format(resultat[l][0],resultat[l][1],resultat[l][2]))

outputfile=prot_id+"res.txt"
wd=os.getcwd()
with open (outputfile, "w") as output:
	output.write(header+"\n")
	output.write(assignement+ "\n")
	output.write("RESIDUE AA STRUCTURE"+ "\n")
	for l2 in range(0,len(resultat)):
		output.write('{}\t{}\t{} \n'.format(resultat[l2][0],resultat[l2][1],resultat[l2][2]))	

source =wd+"/"+outputfile
destination= wd+"/res/"+outputfile
shutil.move(source,destination)
