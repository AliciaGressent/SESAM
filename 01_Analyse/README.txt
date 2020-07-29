Descriptif des scripts du repertoire 01_Analyse:

1) trace_data.r:
objectif: rendre compte de l'échantillonnage des capteurs fixes positionnés à proximité d'une station de référence
-tracer les séries temporelles de la moyenne des trois capteurs en AVG 15min uniquement
-tracer les courbes d'étalonnage aux stations choisies

2) histograms_data.r:
objectifs: visualiser la distribution des concentrations du polluant mesurées aux stations de référence choisies, capteurs fixes et capteurs mobiles
-plot PDF aux stations, moyenne quart horaire des observations de capteurs pour comparer aux stations.
-plot PDF des observations pour la période d'échantillonnage
-plot de l'échantillonnage temporelle pour les stations et les capteurs

3) repeatability.r:
objectifs: calculer la répétabilité moyenne des micro-capteurs à la station de référence choisie afin d'éliminer les valeurs trop faibles qui pourraient apporter de l'erreur dans la fusion de données
-Application de la formule issue du document: http://publications.jrc.ec.europa.eu/repository/handle/JRC83791

