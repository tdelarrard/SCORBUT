# -*- coding: utf-8 -*-
"""Les constantes du problème.

"""


livreExcel = "Données.xlsm"

codenames_sheets = {'Discretisation' : "Discretisation du problème",
                    'CondsInitiales' : "Conditions initiales",
                    'CondsLimites'   : "Conditions limites",
                    'Temperature'    : "Temperature",
                    'PropBeton'      : "Propiétés du béton", 
                    'Solution'       : "Solution"}

noms_donnees = ["Longueur", 
                "Nombre de divisions", 
                "Temps d'analyse", 
                "Nombre de pas d'analyse", 
                "Erreur maximum absolu relative", 
                "Humidité relative initiale", 
                "Humidité relative supposée pour le premier pas", 
                "Type de conditions limites", 
                "Conditions limites du nœud à gauche", 
                "Conditions limites du nœud à droite", 
                "Température", 
                "Masse volumique de l'eau", 
                "Viscosité dynamique de l'eau", 
                "Pression de vapeur saturante", 
                "Masse molaire de l'eau", 
                "Constante des gaz parfaits", 
                "Porosité", 
                "Coefficients du polynôme du taux de la saturation du béton", 
                "Perméabilité à l'eau liquide", 
                "Paramètre de Van Genuchten", 
                "Coefficient de diffusion de la vapeur d'eau hors milieu poreux", 
                "Coefficient A", 
                "Coefficient B"
               ]
M_l = 0.018
R = 8.31