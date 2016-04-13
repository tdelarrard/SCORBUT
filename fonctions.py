"""Fonctions pour lire des fichier d'Excel.

Notes
-----
    Écrire des notes par rapport à ce module.

"""
from numpy import shape, zeros

def table(tetes, donnees, formatDonnees, sep = '', larguer_terminal = 80):
    """Faire une table avec un format personnalisé.
    
    Faire une table avec des têtes et données avec un format personnalisé.
    
    Parameters
    ----------
    tetes : list
        Les têtes de la table.
    donnees : ndarray
        Les données de la table.
    formatDonnees : list
        Le format des données de la table.
    sep : string
        Séparateur des colonnes.
    larguer_terminal : int, optional
        La largeur de la console.
    
    Returns
    -------
    table_formate : list
        La table formaté.
    
    Notes
    -----
    1. Encabezados opcionales.
    2. Cambiar listas y ndarrays por objetos iterables.
    3. Que pasa con encabezados más grandes que el ancho de las columnas
    
    Examples
    --------
    
    """
    # Declaration des variables
    tetes_formates   = []
    donnees_formates = []
    larguerDonnees   = []
    centerDonnees    = []
    table_formate    = []
    
    # Trouver la longuer de chaque colonne
    rem = set('{:+def}')
    
    for j in range(len(formatDonnees)):
        larguerDonnees += [len(formatDonnees[j].format(0))]
    
    # Center les donnes
    for j in range(len(tetes)):
        if len(tetes[j]) < larguerDonnees[j]: centerDonnees += [larguerDonnees[j]]
        else: centerDonnees += [len(tetes[j])]
    
    # Formater les têtes
    for j in range(len(tetes)):
        tetes_formates += [tetes[j].center(centerDonnees[j])]
    
    tetes_formates = [sep.join(tete for tete in tetes_formates).center(larguer_terminal)]
    
    # Formater les données
    for i in range(shape(donnees)[0]):
        dFFilI = []
        for j in range(shape(donnees)[1]):
            dFFilI += [formatDonnees[j].format(donnees[i,j]).center(centerDonnees[j])]
        donnees_formates += [sep.join(donnee for donnee in dFFilI).center(larguer_terminal)]
    
    # Faire la table
    table_formate = tetes_formates + donnees_formates
    
    return table_formate

def masseVolumiqueLEau(T):
    """
    Masse volumique de l'eau liquide en fonction de la temperature.
    
    entre 273.15 K et 2731.5 K en function de la 
    temperature avec la pression externe égal à 1 atm = 101325 Pa.
    
    Paramètres
    ----------
    T : scalaire ou ndarray
        Valeur(s) de la temperature de l'eau [K].
    
    Sorties
    -------
    return : scalaire ou ndarray
        Valeur(s) de la masse volumique de l'eau.
    
    Notes
    -----
    
    """
    
    return - 1.48229e-11 * T ** 6 + 3.00891e-08 * T ** 5 - 2.54839e-05 * T ** 4 \
           + 1.15362e-02 * T ** 3 - 2.94901e+00 * T ** 2 + 4.04138e+02 * T \
           - 2.22010e+04

def presionVapeurSaturante(T):
    """
    Pression de vapeur saturante en fonction de la temperature.
    
    entre 273.15 K et 323.15 K en function de la 
    temperature.
    
    Paramètres
    ----------
    T : scalaire ou ndarray
        Valeur(s) de la temperature de l'eau [K].
    
    Sorties
    -------
    return : scalaire ou ndarray
        Valeur(s) de la pression de vapeur saturante.
    
    Notes
    -----
    
    """
    
    return 4.59689e-09 * T ** 6 - 5.38565e-06 * T ** 5 + 2.51270e-03 * T ** 4 \
         - 5.74269e-01 * T ** 3 + 6.09096e+01 * T ** 2 - 1.56087e+03 * T \
         - 1.21691e+05

def viscositeDynamique(T):
    """
    Viscosité dynamique de l'eau en fonction de la temperature.
    
    entre 273.15 K et 373.15 K en fonction de la 
    temperature.
    
    Paramètres
    ----------
    T : scalaire ou ndarray
        Valeur(s) de la temperature de l'eau [K].
    
    Sorties
    -------
    return : scalaire ou ndarray
        Valeur(s) de la viscosité dynamique de l'eau [Pa·s]
    
    Notes
    -----
    
    """
    
    return 4.50121E-15 * T ** 6 - 9.10566E-12 * T ** 5 + 7.67737E-09 * T ** 4 \
         - 3.45438E-06 * T ** 3 + 8.75136E-04 * T ** 2 - 1.18422E-01 * T \
         + 6.69220E+00

def pol(a, x_i):
    """
    Polynôme de la façon y(x) = sum(a_i x^i) évalué en x = x_i.
    
    Parameters
    ----------
    a : ndarray
        Coefficients du polynôme.
    x_i : float
        Valeur qui prend la variable x.
    
    Returns
    -------
    y : float
        Function y(x) évalué en `x_i`.
    
    Notes
    -----
    --> Vérifié les données d'entrée.
    
    """
    
    y = 0
    for i in range(shape(a)[0]):
        y += a[i] * x_i**i
        
    return y

def solMatTri(A, y):
    """
    Résoudre un système d'équations linéaires tridiagonal de la forme A x = y.
        
    Paramètres
    ----------
    A : ndarray
        Matrice tridiagonal du système.
    y : ndarray
        Matrice colonne second membre connu.
        
    Soties
    ------
    x : ndarray
        Matrice colonne solution.
        
    Notes
    -----
    Impementer une fonction d'error dans le cas de division par zéro.
    """
    
    n = shape(A)[0] # nombre des files
    a, b, c, x = zeros(n), zeros(n), zeros(n), zeros(n) # vectors de zéros
    
    # stoker les valeurs de A en a, b et c
    b[0], b[n - 1] = A[0,0], A[n - 1,n - 1] # diagonal principale
    c[0] = A[0,1] # diagonal au-dessus
    
    for i in range(1, n - 1):
        b[i] = A[i,i] 
        a[i] = A[i,i - 1] # diagonal dessous
        c[i] = A[i,i + 1]
    
    a[n - 1] = A[n - 1,n - 2]
    
    # premier balayage
    D = b[0]
    b[0] = y[0] / D
    a[0] = - c[0] / D
    
    for i in range(1, n):
        D = a[i] * a[i - 1] + b[i]
        b[i] = (y[i] - a[i] * b[i - 1]) / D
        a[i] = - c[i]/D
    
    # second balayage
    x[n - 1] = b[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = a[i] * x[i + 1] + b[i]
    
    return x