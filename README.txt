Descriptif du répertoire data fusion with SEnsorS for Air quality Mapping (SESAM), d'après Gressent et al., (2020)

Ce répertoire contient les scripts permettant de réaliser la fusion de données pour la cartographie de la qualité de l'air à l'échelle urbaine
Dans Gressent et al., 2020, les cartographies ont été réalisées pour les PM10 à Nantes à partir des données de capteurs de la société AtmoTrack, des observations des stations de référence d'Air Pays de la Loire et des sorties du modèle ADMS-Urban utilisées comme dérive pour le krigeage en dérive externe.

Les scripts sont écrits en language R. 
Le krigage est réalisé à partir de la librairie RGeostat (un package pour les applications geostatistiques, MINES ParisTech / ARMINES (2019), RGeostats: The Geostatistical R Package. Version: 11.2.3, Free download from:http://cg.ensmp.fr/rgeostats)

1) 01_Stat:
Analyse exploratoire des données de capteurs, comparaison capteurs fixes / stations de référence
Calcul du seuil de répétabilité pour filtrer le jeu de données

2) 02_Preproc
Application d'un prétraitement sur les données de capteurs
Calcul de la corrélation entre les données de capteurs et la dérive utilisée pour la fusion de données (modèle à l'échelle urbaine)

3) 03_Mapping
Cartographie à l'échelle urbaine à partir d'un krigeage en dérive externe

4) 04_Postproc
Post-traitement des résultats de cartographie

5) INPUTS
Répertoire des fichiers d'entrée

6) OUTPUTS
Répertoire des fichiers de sortie
