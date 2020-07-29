Descriptif des scripts du repertoire 03_Mapping:

1) mapping_urban_scale_EDK.r:
objectif: réaliser une cartographie du polluant choisi à partir des observations capteurs fixes et mobiles à l'échelle urbaine
- application d'un krigeage en dérive externe
- définition de la variance d'erreur de mesure (Variance of Measurement Error, VME, en anglais) qui permet de pondérer les points de mesures en fonction de la dispersion des observations et de l'incertitude de mesures attribuée aux capteurs
- calcul des performance du krigeage
- stockaqe des fichiers de sorties pour le postprocessing
