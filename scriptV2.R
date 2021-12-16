library(ape)
library(xlsx)

#vide l'environnement actuel
rm(list=ls())

#A mettre dans la variable path : le chemin d'accès aux arbres à étudier
path = "C:/Users/Louis/OneDrive/Bureau/projet_orthomam/trees/"

#Les durées de chaque branche
duree_homo = 30
duree_macaca = 30
duree_rattus = 15
duree_mus = 15
duree_rodentia = 75
duree_primate = 60

#Les vecteurs qui stockeront les taux d'évolution des six branches concernées
R_homo <- c()
R_macaca <- c()
R_rattus<- c()
R_mus <- c()
R_primate <- c()
R_rodentia <- c()

#création d'une liste qui contient le nom de fichier de chaque arbre
dir <- list.files(path)
for (file in dir) {
  #print(file)
  tree <- read.tree(paste(path, file, sep="")) # paste(path, file, sep="") permet de créer le chemin d'accès complet à un arbre i : C:/Users/.../Arbre_i
                                               # ensuite on stock l'arbre dans une variable avec read.tree
  
  #Calcul du plus récent ancêtre commun aux primates
  mrca_prim = getMRCA(phy=tree, tip = c("Homo_sapiens", "Macaca_fascicularis"))
  
  #Calcul du plus récent ancêtre commun aux rongeurs 
  mrca_rodent = getMRCA(phy=tree, tip=c("Rattus_norvegicus", "Mus_musculus"))
  
  #Pour obtenir l'ancêtre commun de rodentia et primate (donc Euarchontoglires)
  mrca_rodent_primate = getMRCA(phy=tree, tip = c("Homo_sapiens", "Macaca_fascicularis", "Mus_musculus", "Rattus_norvegicus"))
  
  
  #On récupère l'index de homo, macaca, rattus et mus
  index_homo = which(tree$tip == "Homo_sapiens")
  index_macaca = which(tree$tip == "Macaca_fascicularis")
  index_rattus = which(tree$tip == "Rattus_norvegicus")
  index_mus = which(tree$tip == "Mus_musculus")
  
  #Stock les distances entre toutes les paires de noeuds de l'arbre
  distance = dist.nodes(tree)
  
  #Les distances pour la branche de l'homo, du macaca, du rattus de la mus, de primate et de rodentia
  dist_branche_homo = distance[mrca_prim,index_homo]
  dist_branche_macaca = distance[mrca_prim,index_macaca]
  dist_branche_rattus = distance[mrca_rodent,index_rattus]
  dist_branche_mus = distance[mrca_rodent,index_mus]
  dist_branche_primate = distance[mrca_rodent_primate, mrca_prim]
  dist_branche_rodent = distance[mrca_rodent_primate, mrca_rodent]
  
  #Ici on calcul le taux d'évolution de chaque branche et on le stock dans le vecteur approprié
  Rhomo = dist_branche_homo/duree_homo
  R_homo <- append(R_homo,Rhomo)
  
  Rmacaca = dist_branche_macaca/duree_macaca
  R_macaca <- append(R_macaca, Rmacaca)
  
  Rrattus = dist_branche_rattus/duree_rattus
  R_rattus <- append(R_rattus, Rrattus)
  
  Rmus = dist_branche_mus/duree_mus
  R_mus <- append(R_mus, Rmus)
  
  Rprimate = dist_branche_primate/duree_primate
  R_primate <- append(R_primate, Rprimate)
  
  Rrodentia = dist_branche_rodent/duree_rodentia
  R_rodentia <- append(R_rodentia, Rrodentia)
  
  
}

#Calcul de la médiane pour chacun des 6 taux d'évolution pour l'ensemble des exons
median_homo <- median(x = R_homo)
median_macaca <- median(x = R_macaca)
median_rattus <- median(x = R_rattus)
median_mus <- median(x = R_mus)
median_primate <- median(x = R_primate)
median_rodentia <- median(x = R_rodentia)

#Calcul des quantiles à 5% et 95%
quantiles_homo = quantile(x = R_homo, probs = seq(0.05, 1, 0.90), names = F)
quantiles_macaca = quantile(x = R_macaca, probs = seq(0.05, 1, 0.90), names = F)
quantiles_rattus = quantile(x = R_rattus, probs = seq(0.05, 1, 0.90), names = F)
quantiles_mus = quantile(x = R_mus, probs = seq(0.05, 1, 0.90), names = F)
quantiles_primate = quantile(x = R_primate, probs = seq(0.05, 1, 0.90), names = F)
quantiles_rodentia = quantile(x = R_rodentia, probs = seq(0.05, 1, 0.90), names = F)

# 1.B : calcul des la mediane et des quantiles des branches
dist_homo <- c(quantile(R_homo,probs = 0.05),median(R_homo),quantile(R_homo,probs = 0.95))
dist_macaca <- c(quantile(R_macaca,probs = 0.05),median(R_macaca),quantile(R_macaca,probs = 0.95))
dist_rattus <- c(quantile(R_rattus,probs = 0.05),median(R_rattus),quantile(R_rattus,probs = 0.95))
dist_mus <- c(quantile(R_mus,probs = 0.05),median(R_mus),quantile(R_mus,probs = 0.95))
dist_primate <- c(quantile(R_primate,probs = 0.05),median(R_primate),quantile(R_primate,probs = 0.95))
dist_rodentia <- c(quantile(R_rodentia,probs = 0.05),median(R_rodentia),quantile(R_rodentia,probs = 0.95))

#Analyse des données de tous les groupes de la promo : 
#Dans monDataset on stock les valeurs récupérées par la promo
monDataset <- read.xlsx("C:/Users/Louis/OneDrive/Bureau/projet_orthomam/projet_phylogenie.xlsx", sheetName = 2)

#Pour afficher sur un seul graphe les valeurs de la souris et de l'homme en fonction du %G-C3 min il faut construire
#un dataframe pour chaque catégorie de valeurs (médiane, quantile 5%, quantile 95%). Chaque dataframe aura donc 
#26 lignes et 3 colonnes 
#dans le cas de la médiane : chaque valeur de médiane sera associé à un %G-C3 et à un label "homme" ou "souris"
GC3_min = vector(length = 26)
for(i in 1:26){
  if(i <14){
    GC3_min[i] = monDataset[i, 2]
  }
  else{
    GC3_min[i] = monDataset[i-13,2]
  }
}

label = vector(length = 26)
for(i in 1:26){
  if(i<14){
    label[i] = "souris"
  }
  else{
    label[i] = "homme"
  }
}

Mediane_Taux_Rs = vector(length = 26)
for(i in 1:26){
  if(i < 13){
    Mediane_Taux_Rs[i] = monDataset[i,4]
  }
  else{
    Mediane_Taux_Rs[i] = monDataset[i-13,7]
  }
}
df1 = data.frame(GC3_min, Mediane_Taux_Rs, label)
#On utilise ggplot2 pour plot nos valeurs ensemble.
plot1 <- ggplot(df1, aes(x = GC3_min, y = Mediane_Taux_Rs, colour = factor(label))) + geom_point(size=2.5)


Quantile_5_Taux_Rs = vector(length = 26)
for(i in 1:26){
  if(i < 14){
    Quantile_5_Taux_Rs[i] = monDataset[i,5]
  }
  else{
    Quantile_5_Taux_Rs[i] = monDataset[i-13, 8]
  }
}
df2 = data.frame(GC3_min, Quantile_5_Taux_Rs, label)
plot2 <- ggplot(df2, aes(x = GC3_min, y = Quantile_5_Taux_Rs, colour = factor(label))) + geom_point(size=2.5)



Quantile_95_Taux_Rs = vector(length = 26)
for(i in 1:26){
  if(i < 14){
    Quantile_95_Taux_Rs[i] = monDataset[i,6]
  }
  else{
    Quantile_95_Taux_Rs[i] = monDataset[i-13, 9]
  }
}
df3 = data.frame(GC3_min, Quantile_95_Taux_Rs, label)
plot3 <- ggplot(df2, aes(x = GC3_min, y = Quantile_95_Taux_Rs, colour = factor(label))) + geom_point(size=2.5)




#On enlève de l'environnement les variables qui ne servent plus
rm(list = c("distance", 
            "dist_branche_homo", "dist_branche_macaca", "dist_branche_mus", "dist_branche_primate", "dist_branche_rattus", "dist_branche_rodent", 
            "duree_homo", "duree_macaca", "duree_mus", "duree_rattus", "duree_primate", "duree_rodentia",
            "file", "dir", "path",
            "index_homo", "index_macaca", "index_mus", "index_rattus",
            "mrca_prim", "mrca_rodent", "mrca_rodent_primate",
            "Rhomo", "Rmacaca", "Rmus", "Rrattus", "Rprimate", "Rrodentia"))