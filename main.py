# -*- coding: utf-8 -*-
"""Résoudre le problème de séchage du béton



Notes
-----
    Agregar descripción de este módulo.

"""

import os
from numpy import linspace, savetxt, shape
from discretisation import *
from sechage import resSec
from xlwings import Workbook, Range
from constantes import codenames_sheets, livreExcel
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Résoudre le problème d'humidité et séchage
    h_all_all = resSec(L, J, T, N, erreur_max, h_all_0, h_all_1_0, condLim, 
                       condLim_0_all, condLim_L_all, temp, p_l, n_l, Pvsat, M_l, R, phi, 
                       Sr, k_l, m, D_0, A, B)
    
    # Enregistrer la solution
    savetxt('solution.txt', h_all_all, fmt = '%10.4E', delimiter = '\t', 
            newline =  os.linesep)

    # Début : Plot des solutions
    plt.figure()
    x = linspace(0, L, J + 1) # Coordonnées de l'axe x 
    for i in range(1, int(shape(h_all_all)[1])): # pour chaque pas du temps
        plt.plot(x, h_all_all[:,i], linewidth = 0.5)
    
    # Enregistrer la solution en Solution
    wb = Workbook.active()
    Range(codenames_sheets['Solution'], "C3").value = h_all_all
    
    # title, legends, etc ...
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel(r'x \text{[m]}')
    plt.xlim(0,L)
    plt.ylabel(r'h \text{[-]}')
    plt.title(u"Evolution de l'humidite relative")
    plt.grid(True)
    plt.show()
    # Fin : Plot des solutions