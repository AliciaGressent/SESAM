Descriptif des scripts du repertoire 02_Preproc:

1) formate_drift.r:
objectif: créer un fichier csv à partir du fichier .asc ADMS-Urban fourni par Air Pays de la Loire
-création d'un raster
-plot du raster et visualisation 
-sauvegarde du fichier au format csv. 

2) preproc_MS.r:
objectif: nettoyer et corriger les données de capteurs mobiles
-Eliminer les données en-dessous du seuil de répétabilité calculé dans 01_Analyses/repeatability.r
-Lissage des données si nécessaire => trop de bruit: moyenne (pas nécessaire pour les données mobiles)
-Correction de la variation journalière des concentrations de fond: Cf <- (Ci- BG_DATA) + BG_REF avec BG_DATA = la concentration médiane claculée sur la fenêtre glissante de 15 min et BG_REF = la moyenne journalière des observations des stations de référence du réseau
-La correction est réalisée pour chaque run de mesures <=> série continue de mesures (max 5 minutes entre deux mesures)

3) preproc_FS.r:
objectif: nettoyer et corriger les données de capteurs fixes
-Eliminer les données en-dessous du seuil de répétabilité calculé dans 01_Analyses/repeatability.r
-Lissage des données si nécessaire => trop de bruit: moyenne (pas nécessaire pour les données mobiles)
-Correction de la variation journalière des concentrations de fond: Cf <-
(Ci-BG_DATA) + BG_REF avec BG_DATA = la concentration médiane claculée sur la
fenêtre glissante de 15 min et BG_REF = la moyenne journalière des observations de la station de référence la plus proche
-La correction est réalisée pour chaque run de mesures <=> série continue de mesures (max 5 minutes entre deux mesures)

4) compare_data_model.r:
objectif: comparer les données de capteurs et de stations de référence aux sorties du modèle urbain
-comparaison des données de capteurs brutes et corrigées, des données de
référence au modèle à la résolution horaire sur le mois entier de données (période d'échantillonnage disponible)
-tracer les séries temporelles et les pdf sur un même plot

5) correlation_data_drift.r:
objectif: vérifier la corrélation entre les données de capteurs et la dérive (ici sortie de modèle urbain)
-sélection point de grille le plus proche du point de mesure capteur
-calcul de la corrélation


