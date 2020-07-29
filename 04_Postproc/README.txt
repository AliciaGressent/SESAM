Descriptif des scripts du repertoire 04_Postproc:

1) compare_pred_obs.r: 
objectif: comparer aux stations de référence l'estimation par krigeage en dérive externe et les observations (résolution horaire et journalière)
-plot la comparaison à chaque station du réseau de l'AASQA: observations référence, fusion de données et modèle urbain
-plot carte de la moyenne journalière

2) mapping_drift_obs_pred_vme_stdev.r:
objectif: plot pour une échéance horaire choisie les résultats de la fusion de données
-plot observations, VME, estimation et erreur associée

3) mapping_hourly_pred.r: 
objectif: plot pour toutes les échéances horaires l'estimation par fusion de données
-unique plot avec toutes les échéances entre 7h et 19h

4) mapping_hourly_vme.r: 
objectif: plot pour toutes les échéances horaires la VME
-unique plot avec toutes les échéances entre 7h et 19h

5) mapping_hourly_error.r: 
objectif: plot pour toutes les échéances horaires l'erreur de krigeage
-unique plot avec toutes les échéances entre 7h et 19h

6) mapping_hourly_obs.r:
objectif: plot pour toutes les échéances horaires les observations capteurs
-unique plot avec toutes les échéances entre 7h et 19h  
