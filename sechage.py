"""Fonctions pour résoudre le problème d'humidité et séchage.

Ce module a les fonctions pour résoudre le problème d'humidité et séchage.

Misc Functions
--------------
- `k_rl`   -- Calculer la perméabilité relative à l'eau liquide.
- `Rdg`    -- Calculer la résistance du milieu poreux à la diffusion des gaz.
- `assSys` -- Assembler le système des équations du problème de séchage-humidité.
- `resSec` -- Résoudre le problème de séchage-humidité.

References
----------
.. [1] Jonathan MAI-NHU. Corrosion des armatures du béton: couplage 
       carbonatation/chlorures en présence de cycles hydriques, 2013.
.. [2] Cristian Ramirez. Stage de fin d'études. Modélisation numérique du 
       couplage chloration/carbonatation sous conditions atmosphériques en 
       milieu poreux non saturé, 2016.

"""

import matplotlib.pyplot as plt
from numpy import absolute, append, array, linspace, ones, savetxt, shape, \
                  zeros
from fonctions import pol, solMatTri

__all__ = ['k_rl', 'Rdg', 'assSys', 'resSec']

def k_rl(Sr, m):
    """
    Calculer la perméabilité relative à l'eau liquide (équation 1.82 de [1]).
    
    Parameters
    ----------
    Sr : float
        Taux de la saturation du béton [-].
    m : float
        Paramètre de Van Genuchten [-].
        
    Returns
    -------
    return : float
        Perméabilité relative à l'eau liquide [-].
        
    References
    ----------
    .. [1] Jonathan MAI-NHU. Corrosion des armatures du béton: couplage 
           carbonatation/chlorures en présence de cycles hydriques, 2013.
    
    """
    
    return Sr**0.5 * (1 - (1 - Sr**(1/m))**m)**2

def Rdg(phi, Sr, A, B):
    """
    Calculer la résistance du milieu poreux à la diffusion des gaz (éq. 1.85 [1]).
    
    Parameters
    ----------
    phi : float
        Porosité [-].
    Sr : ndarray
        Taux de la saturation du béton [-].
    A : scalaire
        Coefficient A [-].
    B : scalaire
        Coefficient B [-].
    
    Returns
    -------
    return : float
        Résistance du milieu poreux à la diffusion des gaz [-].
    
    References
    ----------
    .. [1] Jonathan MAI-NHU. Corrosion des armatures du béton: couplage 
           carbonatation/chlorures en présence de cycles hydriques, 2013.
    
    """
    return phi**A * (1 - Sr)**B

def _A(p_l, phi, dSrdh, Pvsat, M_l, R, T, Sr, h):
    """
    Variable A(h, a) de l'équation 1.4 de [1].
    
    Parameters
    ----------
    p_l : float
        Masse volumique de l'eau [kg/m³].
    phi : float
        Porosité [-].
    dSrdh : float
        Devivée de la taux de la saturation du béton [-].
    Pvsat : float
        Pression de vapeur saturante [Pa].
    M_l : float
        Masse molaire de l'eau [kg/mol].
    R : float
        Constante des gaz parfaits [J/mol·K].
    T : float
        Température [K].
    Sr : float
        Taux de la saturation du béton [-].
    h : float
        Humidité relative [-].
    
    Returns
    ------
    return : float
        Coefficient A(h, a) [kg/m³].
    
    References
    ----------
    .. [1] Cristian Ramirez. Stage de fin d'études. Modélisation numérique du 
           couplage chloration/carbonatation sous conditions atmosphériques en 
           milieu poreux non saturé, 2016.
    
    """
    return p_l * phi * dSrdh + (Pvsat * M_l)/(R * T) * phi * (1 - Sr - h * dSrdh)
    
def _B(p_l, k_l, n_l, k_rl, R, T, M_l, h):
    """
    Variable B(h, b) de l'équation 1.4 de [1].
    
    Parameters
    ----------
    p_l : float
        Masse volumique de l'eau [kg/m³].
    k_l : float
        Perméabilité à l'eau liquide [m²].
    n_l : float
        Viscosité dynamique de l'eau [Pa·s].
    k_rl : float
        Perméabilité relative à l'eau liquide [-].
    R : float
        Constante des gaz parfaits [J/mol·K].
    T : float
        Temperature [K].
    M_l : float
        Masse molaire de l'eau [kg/mol].
    h : float
        Humidité relative [-].
        
    Returns
    -------
    return : scalaire
        Coefficient B(h, b) [kg/m·s].
    
    References
    ----------
    .. [1] Cristian Ramirez. Stage de fin d'études. Modélisation numérique du 
           couplage chloration/carbonatation sous conditions atmosphériques en 
           milieu poreux non saturé, 2016.
    
    """
    return -p_l**2 * (k_l/n_l) * k_rl * (R * T/M_l) * (1/(h**2))
    
def _C(p_l, k_l, n_l, k_rl, R, T, M_l, h, Pvsat, R_dg, D_0):
    """
    Variable C(h, c) de l'équation 1.4 de [1].
    
    Parameters
    ----------
    p_l : float
        Masse volumique de l'eau [kg/m³].
    k_l : float
        Perméabilité à l'eau liquide [m²].
    n_l : float
        Viscosité dynamique de l'eau [Pa·s].
    k_rl : float
        Perméabilité relative à l'eau liquide [-].
    R : float
        Constante des gaz parfaits [J/mol·K].
    T : float
        Température [K].
    M_l : float
        Masse molaire de l'eau [kg/mol].
    h : float
        Humidité relative [-].
    Pvsat : float
        Pression de vapeur saturante [Pa].
    R_dg : float
        Résistance du milieu poreux à la diffusion des gaz [-].
    D_0 : float
        Coefficient de diffusion de la vapeur d'eau hors milieu poreux [m²/s].
        
    Returns
    -------
    return : float
       Coefficient C(h,c) [kg/m·s].
    
    References
    ----------
    .. [1] Cristian Ramirez. Stage de fin d'études. Modélisation numérique du 
           couplage chloration/carbonatation sous conditions atmosphériques en 
           milieu poreux non saturé, 2016.
    
    """
    return p_l**2 * (k_l/n_l) * k_rl * ((R * T) / M_l) * (1/h) + \
           Pvsat * (M_l / (R * T)) * R_dg * D_0

def B_inf_sui(dt, dx, B_dem_act, h_inf_sui, C_dem_act):
    """
    Coefficient B_{i - 1}^{j + 1} de l'équation 1.7 de [1].
    
    Parameters
    ----------
    dt : float
        Taille du pas du temps [s].
    dx : float
        Taille du pas du space [m].
    B_dem_act : float
        Coefficient B(h, b) évalué en h(i, j) [kg/m·s].
    h_inf_sui : float
        Humidité relative en (i - 1, j + 1) [-].
    C_dem_act :
        Coefficient C(h, c) évalué en h(i, j) [kg/m·s].
    
    Returns
    -------
    return : float
        Coefficient B_{i - 1}^{j + 1} [kg/m³].
    
    References
    ----------
    .. [1] Cristian Ramirez. Stage de fin d'études. Modélisation numérique du 
           couplage chloration/carbonatation sous conditions atmosphériques en 
           milieu poreux non saturé, 2016.
    
    """
    return -(dt / (4*dx**2)) * (B_dem_act * h_inf_sui + 4*C_dem_act)

def B_dem_sui(A_dem_act, dt, C_dem_act, dx):
    """
    Coefficient B_i^{j + 1} de l'équation 1.7 de [1].
    
    Parameters
    ----------
    A_dem_act : float
        Coefficient A(h, b) évalué en h(i, j) [kg/m·s].
    dt : float
        Taille du pas du temps [s].
    C_dem_act : float
        Coefficient C(h, c) évalué en h(i, j) [kg/m·s].
    dx : float
        Taille du pas du space [m].
    
    Returns
    -------
    return : float
        Coefficient B_i^{j + 1} [kg/m³].
    
    References
    ----------
    .. [1] Cristian Ramirez. Stage de fin d'études. Modélisation numérique du 
           couplage chloration/carbonatation sous conditions atmosphériques en 
           milieu poreux non saturé, 2016.
    
    """
    return A_dem_act + 2*(dt/(dx**2)) * C_dem_act

def B_sup_sui(dt, dx, B_dem_act, h_sup_sui, h_inf_sui, C_dem_act):
    """    
    Coefficient B_{i + 1}^{j + 1} de l'équation 1.7 de [1].
    
    Parameters
    ----------
    dt : float
        Taille du pas du temps [s].
    dx : float
        Taille du pas du space [m].
    B_dem_act : float
        Coefficient B(h, b) évalué en h(i, j) [kg/m·s].
    h_sup_sui : float
        Humidité relative en (i + 1, j + 1) [-].
    h_inf_sui : float
        Humidité relative en (i - 1, j + 1) [-].
    C_dem_act : float
        Coefficient C(h, c) évalué en h(i, j) [kg/m·s].
    
    Returns
    -------
    return : float
        Coefficient B_{i + 1}^{j + 1} [kg/m³].
    
    References
    ----------
    .. [1] Cristian Ramirez. Stage de fin d'études. Modélisation numérique du 
           couplage chloration/carbonatation sous conditions atmosphériques en 
           milieu poreux non saturé, 2016.
    
    """
    return -(dt / (4*dx**2)) * (B_dem_act * (h_sup_sui - 2*h_inf_sui) + 4*C_dem_act)

def assSys(L, J, dt, h_all_act, h_all_sui_0, condLim, condLim_0_sui, condLim_L_sui, 
           temp, p_l, n_l, Pvsat, M_l, R, phi, Sr, k_l, m, D_0, A, B):
    """
    Assembler le système des équations du problème de séchage-humidité 1.3 de [1].
    
    Le système des équations est donné pour l'équation (1.10) de [1] où :
    
	\beta_{i-1}^{j+1} h_{i-1}^{j+1} + 
    \beta_i^{j+1} h_i^{j+1} + 
    \beta_{i+1}^{j+1} h_{i+1}^{j+1} = 
    A h_i^j
    
    Pour assembler le système il faut supposer la solution de l'humidité relative 
    pour le temps suivant.
        
    Parameters
    ----------
    L : float
        Longueur [m].
    J : int
        Nombre de divisions [-].
    dt : float
        Taille du pas du temps [s].
    h_all_act : ndarray
        Humidité relative actuel [-].
    h_all_sui_0 : ndarray
        Humidité relative suivante supposée [-].
    condLim : ndarray
        Type de conditions limites [-].
        a, b = 0 : Condition de frontière type Dirichlet.
        a, b = 1 : Condition de frontière type Neumann.
            (a, b) :
            a : int
                Condition limite en x = 0.
            b : int
                Condition limite en x = L.
    condLim_0_sui : float
        Condition limite suivante à gauche [-] ou [m⁻¹].
    condLim_L_sui : float
        Condition limite suivante à droite [-] ou [m⁻¹].
    temp : float
        Température [K].
    p_l : float
        Masse volumique de l'eau [kg/m³].
    n_l : float
        Viscosité dynamique de l'eau [Pa·s].
    Pvsat : float
        Pression de vapeur saturante [Pa].
    M_l : float
        Masse molaire de l'eau [kg/mol].
    R : float
        Constante des gaz parfaits [J/mol·K].
    phi : float
        Porosité [-].
    Sr : ndarray
        Coefficients du polynôme du taux de la saturation du béton [-].
    k_l : float
        Perméabilité à l'eau liquide [m²].
    m : float
        Paramètre de Van Genuchten [-].
    D_0 : float
        Coefficient de diffusion de la vapeur d'eau hors milieu poreux [m²/s].
    A : float
        Coefficient A [-].
    B : float
        Coefficient B [-].
    
    Returns
    -------
    K : ndarray
        Matrice des coefficients [kg/m³].
    F : ndarray
        Vecteur colonne des coefficients [kg/m³].
    
    Notes
    -----
    --> 
            
    References
    ----------
    .. [1] Cristian Ramirez. Stage de fin d'études. Modélisation numérique du 
           couplage chloration/carbonatation sous conditions atmosphériques en 
           milieu poreux non saturé, 2016.
    
    """
    # Définir des autres variables
    
    N = int((J - 1) + condLim[0] + condLim[1]) # : taille de la matrice des coeff.
    K = zeros((N, N)) # : matrice des coefficients
    F = zeros((N, 1)) # : vecteur colonne des coefficients
    
    dx = L / J # : taille du pas dans l'espace

    dSrdh = array([ ]) # : Coefficients de la derivée du taux de saturation
    for i in range(1, shape(Sr)[0]):
        dSrdh = append(dSrdh, i * array([Sr[i]]))
    
    # Remplir la matrice des coefficients
    
    # Remplir l'extreme gauche.
    if condLim[0] == 0:
        # if condLim[0] == 0: k = 1
        k = 1 # Définir le noeud suivant
        # Valeurs dans le noeud `1`
        h_1_act     = h_all_act[k] # : humidité relative
        Sr_1_act    = pol(Sr, h_1_act) # : taux de saturation du béton
        dSrdh_1_act = pol(dSrdh, h_1_act) # : dérivée du taux de saturation du béton
        k_rl_1_act  = k_rl(Sr_1_act, m) # : pérmeabilité relative à l'eau liquide
        Rdg_1_act   = Rdg(phi, Sr_1_act, A, B) # : résistance du milieu poreux ...
        # Coefficients de la équation 1.4 de [1]
        A_1_act = _A(p_l, phi, dSrdh_1_act, Pvsat, M_l, R, temp, Sr_1_act, h_1_act)
        B_1_act = _B(p_l, k_l, n_l, k_rl_1_act, R, temp, M_l, h_1_act)
        C_1_act = _C(p_l, k_l, n_l, k_rl_1_act, R, temp, M_l, h_1_act, Pvsat, 
                     Rdg_1_act, D_0)
        # Humidité relative suivante supposée dans les noeuds `0` et `2`
        h_0_sui = condLim_0_sui
        h_2_sui = h_all_sui_0[k + 1]
        # Coefficients de l'équation 1.7 de [1]
        B_0_sui = B_inf_sui(dt, dx, B_1_act, h_0_sui, C_1_act)
        B_1_sui = B_dem_sui(A_1_act, dt, C_1_act, dx)
        B_2_sui = B_sup_sui(dt, dx, B_1_act, h_2_sui, h_0_sui, C_1_act)
        # Valeurs des éléments `(0, 0)` et `(0, 1)` de la matrice
        K[0, 0] = B_1_sui
        K[0, 1] = B_2_sui
        # Valeur de l'élément (0,0) du vecteur
        F[0, 0] = A_1_act * h_1_act - B_0_sui * h_0_sui
        
    elif condLim[0] == 1:
        input("condLim[0] == 1 (True) : Je ne suis pas sûr de cette partie :p")
        input("Faire attention !")
        # # valeurs dans le noeud 0
        # h_0_act     = h_all_act[0] # humidité relative
        # Sr_0_act    = pol(Sr, h_0_act) # degré de saturation
        # dSrdh_0_act = pol(dSrdh, h_0_act) # dérivée humidité relative 
        # k_rl_0_act  = k_rl(Sr_0_act, m) # pérmeabilité relative à la phase liquide
        # # coefficients de la équation 3.10
        # A_0_act = A(p_l, phi, dSrdh_0_act, Pvsat, M_l, R, temp, Sr_0_act, h_0_act)
        # B_0_act = B(p_l, k_l, n_l, k_rl_0_act, R, temp, M_l, h_0_act)
        # C_0_act = C(p_l, k_l, n_l, k_rl_0_act, R, temp, M_l, h_0_act, R_dg, D_0)
        # # valeurs supposées dans les noeuds 0 et 1
        # h__1_sui = h_all_sui[1]
        # h_1_sui  = h_all_sui[1]
        # # coefficients de l'équation 3.10 après la discretisation
        # B__1_sui = B_inf_sui(dt, dx, B_0_act, h__1_sui, C_0_act)
        # B_1_sui  = B_sup_sui(dt, dx, B_0_act, h_1_sui, h__1_sui, C_0_act)
        # # Valeurs des éléments (0, 0) et (0, 1) de la matrice
        # K[0, 0] = B_dem_sui(A_0_act, dt, C_0_act, dx)
        # K[0, 1] = B__1_sui + B_1_sui
        
    # Remplir les noeuds internes
    for i in range(1, N - 1):
        # Valeurs dans le noeud i
        h_i_act     = h_all_act[i + k] # : humidité relative
        Sr_i_act    = pol(Sr, h_i_act) # : taux de saturation du béton
        dSrdh_i_act = pol(dSrdh, h_i_act) # : dérivée du taux de saturation du béton
        k_rl_i_act  = k_rl(Sr_i_act, m) # : pérmeabilité relative à l'eau liquide
        Rdg_i_act   = Rdg(phi, Sr_i_act, A, B) # : résistance du milieu poreux ...
        # Coefficients de la équation 1.4 de [1]
        A_i_act = _A(p_l, phi, dSrdh_i_act, Pvsat, M_l, R, temp, Sr_i_act, h_i_act)
        B_i_act = _B(p_l, k_l, n_l, k_rl_i_act, R, temp, M_l, h_i_act)
        C_i_act = _C(p_l, k_l, n_l, k_rl_i_act, R, temp, M_l, h_i_act, Pvsat, 
                     Rdg_i_act, D_0)
        # Humidité relative suivante supposée dans les noeuds `i - 1` et `i + 1`
        h_inf_sui = h_all_sui_0[i - 1 + k]
        h_sup_sui = h_all_sui_0[i + 1 + k]
        # Valeurs des éléments `(i, i - 1)`, `(i, i)` et `(i, i + 1)` de la matrice
        K[i, i - 1] = B_inf_sui(dt, dx, B_i_act, h_inf_sui, C_i_act)
        K[i, i]     = B_dem_sui(A_i_act, dt, C_i_act, dx)
        K[i, i + 1] = B_sup_sui(dt, dx, B_i_act, h_sup_sui, h_inf_sui, C_i_act)
        # Valeur de l'élément `(0, 0)` du vecteur
        F[i, 0] = A_i_act * h_i_act
    
    # Remplir l'extreme droit
    if condLim[1] == 0:
        raw_input("condLim[1] == 0 (True) : Je ne suis pas sûr de cette partie :p")
        raw_input("Faire attention !")
        # # valeurs dans le noeud N - 1
        # h_N__1_act     = h_all_act[N - 1] # humidité relative
        # Sr_N__1_act    = pol(Sr, h_N__1_act) # degré de saturation
        # dSrdh_N__1_act = pol(dSrdh, h_N__1_act) # dérivée humidité relative 
        # k_rl_N_1__act  = k_rl(Sr_N__1_act, m) # pérmeabilité relative
        # # coefficients de la équation 3.10
        # A_N__1_act = A(p_l, phi, dSrdh_N__1_act, Pvsat, M_l, R, temp, Sr_N__1_act, 
                       # h_N__1_act)
        # B_N__1_act = B(p_l, k_l, n_l, k_rl_N__1_act, R, temp, M_l, h_N__1_act)
        # C_N__1_act = C(p_l, k_l, n_l, k_rl_N__1_act, R, temp, M_l, h_N__1_act, R_dg, 
                       # D_0)
        # # valeurs supposées dans les noeuds N - 2 et N
        # h_N__2_sui = h_all_sui[N - 2]
        # h_N_sui    = h_all_sui[N]
        # # valeurs des éléments (N, N - 1) et (N, N) de la matrice
        # K[N - 1, N] = B_inf_sui(dt, dx, B_N__1_act, h_N__2_sui, C_N__1_act)
        # K[N, N]     = B_dem_sui(A_N__1_act, dt, C_N__1_act, dx)
    elif condLim[1] == 1:
        # Valeurs dans le noeud `J`
        h_J_act     = h_all_act[J] # : humidité relative
        Sr_J_act    = pol(Sr, h_J_act) # : taux de saturation du béton
        dSrdh_J_act = pol(dSrdh, h_J_act) # : dérivée du taux de saturation du béton
        k_rl_J_act  = k_rl(Sr_J_act, m) # : pérmeabilité relative à l'eau liquide
        Rdg_J_act   = Rdg(phi, Sr_J_act, A, B) # : résistance du milieu poreux ...
        # Coefficients de la équation 1.4 de [1]
        A_J_act = _A(p_l, phi, dSrdh_J_act, Pvsat, M_l, R, temp, Sr_J_act, h_J_act)
        B_J_act = _B(p_l, k_l, n_l, k_rl_J_act, R, temp, M_l, h_J_act)
        C_J_act = _C(p_l, k_l, n_l, k_rl_J_act, R, temp, M_l, h_J_act, Pvsat, 
                     Rdg_J_act, D_0)
        # Humidité relative suivante supposée dans les noeuds `J - 1` et `J + 1`
        h_J__1_sui = h_all_sui_0[J - 1]
        h_J_1_sui  = h_J__1_sui + 2*dx*condLim_L_sui
        # Coefficients de l'équation 1.7 de [1]
        B_J__1_sui = B_inf_sui(dt, dx, B_J_act, h_J__1_sui, C_J_act)
        B_J_sui    = B_dem_sui(A_J_act, dt, C_J_act, dx)
        B_J_1_sui  = B_sup_sui(dt, dx, B_J_act, h_J_1_sui, h_J__1_sui, C_J_act)
        # Valeurs des éléments `N, N - 1` et `N, N` de la matrice
        K[N - 1, N - 2] = B_J__1_sui + B_J_1_sui
        K[N - 1, N - 1] = B_J_sui
        # valeur de l'élément`N - 1, 0` du vecteur colonne
        F[N - 1, 0] = A_J_act * h_J_act - 2*B_J_1_sui*dx*condLim_L_sui
    
    return K, F
    
def resSec(L, J, T, N, erreur_max, h_all_0, h_all_1_0, condLim, condLim_0_all,
           condLim_L_all, temp, p_l, n_l, Pvsat, M_l, R, phi, Sr, k_l, m, D_0, 
           A, B):
    """
    Résoudre le problème de séchage-humidité.
    
    Le problème se résoudre avec l'équation 3.10 de [1] en supposant que les 
    variations temporelles de `p_l`, `k_l`, 'n_l`, `k_rl`, `R`, `temp`, `M_l`, 
    `Pvsat`, `R_dg` et `D_0`sont négligeables. C'est-à-dire, avec l'équation 1.3
    de [2], ou :
    
    \left[ \rho_{l} \phi \frac{\partial S_r}{\partial h} 
          + \frac{P_{vsat} M_v}{R T} \phi 
          \left( 1 - S_r - h \frac{\partial S_r}{\partial h} \right)
    \right] \frac{\partial h}{\partial t}
    = 
    -\rho_l^2 \frac{k_l}{\eta_l} k_{rl} \frac{R T}{M_v} \frac{1}{h^2} 
    \left( \frac{\partial h}{\partial x} \right)^2 
	+ \left( \rho_l^2 \frac{k_l}{\eta_l} k_{rl} \frac{R T}{M_v} \frac{1}{h} 
            + P_{vsat} \frac{M_v}{R T} R_{dg} D_0 
      \right) \frac{\partial ^2 h}{\partial x^2}
    
    Pour résoudre l'équation on suppose, pour les dérivées de premier ordre, 
    différence finies vers en avant et, pour les dérivées de deuxième ordre, 
    différence finies centrées.    
    
    Parameters
    ----------
    L : float
        Longueur [m].
    J : int
        Nombre de divisions [-].
    T : ndarray
        Temps d'analyse [s].
    N : int
        Nombre de temps d'analyse [-].
    erreur_max : float
        Erreur maximum absolu relative [-].
    h_all_0 : ndarray
        Humidité relative initiale [-].
    h_all_1_0 : ndarray
        Humidité relative supposée pour le premier pas [-].
    condLim : ndarray
        Type de conditions limites [-].
        a, b = 0 : Condition de frontière type Dirichlet.
        a, b = 1 : Condition de frontière type Neumann.
            (a, b) :
            a : int
                Condition limite en x = 0.
            b : int
                Condition limite en x = L.
    condLim_0_all : ndarray
        Condition limite à gauche [-] ou [m⁻¹].
    condLim_L_all : ndarray
        Condition limite à droite [-] ou [m⁻¹].
    temp : ndarray
        Température [K].
    p_l : ndarray
        Masse volumique de l'eau [kg/m³].
    n_l : ndarray
        Viscosité dynamique de l'eau [Pa·s].
    Pvsat : ndarray
        Pression de vapeur saturante [Pa].
    M_l : float
        Masse molaire de l'eau [kg/mol].
    R : float
        Constante des gaz parfaits [J/mol·K].
    phi : float
        Porosité [-].
    Sr : ndarray
        Coefficients du polynôme du taux de la saturation du béton [-].
    k_l : float
        Perméabilité à l'eau liquide [m²].
    m : float
        Paramètre de Van Genuchten [-].
    D_0 : float
        Coefficient de diffusion de la vapeur d'eau hors milieu poreux [m²/s].
    A : float
        Coefficient A [-].
    B : float
        Coefficient B [-].
    
    Returns
    ------
    h_all_all : ndarray
        Solution h(x, t).
    
    Notes
    -----
    --> Cambiar la entrada Sr para poder recibir una función. Permitiría tener en 
        cuenta más que sólo polinomios.
    --> ¿Es que acaso k_l, m. D_0, A y B cambian en el tiempo?
    --> Posibilité de prendre en compte la variation de la porosité dans l'espace ?
    
    References
    ----------
    .. [1] Jonathan MAI-NHU. Corrosion des armatures du béton: couplage carbonatation/chlorures
		   en présence de cycles hydriques, 2013.
    .. [2] Cristian Ramirez. Stage de fin d'études Modélisation numérique du 
           couplage chloration/carbonatation sous conditions atmosphériques en 
           milieu poreux non saturé, 2016.
    
    """
    
    # Des autres variables
    dt = array(T[1:] - T[0:-1]) # : pas du temps
    h_all_all = zeros((J + 1, N + 1)) # : solution h(x, t)
    h_all_all[:, 0] = h_all_0 # enregistrer la solution h(x, 0)
    h_all_all[:, 1] = h_all_1_0 # enregistrer l'hypothèse h(x, 1)
    
    # Résoudre pour chaque pas du temps

    # Itérer dans chaque pas du temps
    for t in range(N): 
        
        print("----- Pas du temps numéro : {} -----".format(t + 1))
        j = 1 # Nombre d'itération
        print("********** Itération : {} **********".format(j))
        
        # Assambler la matrice de coefficients et le vecteur colonne
        K, F = assSys(L, J, dt[t], h_all_all[:, t], h_all_all[:, t + 1], condLim, 
                      condLim_0_all[t + 1], condLim_L_all[t + 1], temp[t], p_l[t], 
                      n_l[t], Pvsat[t], M_l, R, phi, Sr, k_l, m, D_0, A, B)

        # Résoudre l'humidite
        h_all_sui_i = solMatTri(K, F)
        if condLim[0] == 0:
            h_all_sui_i = append(condLim_0_all[t + 1], h_all_sui_i)
        elif condLim[1] == 0:
            h_all_sui_i = append(h_all_sui_i, condLim_L_all[t + 1])

        # Itérer jusqu'à trouver la bonne solution

        erreur = True # Inicialment force l'erreur à etrê vrai
        # Itérer jusqu'à trouver erreur faux
        while erreur:
            i = 0 # : noeud
            # Verifier l'erreur de chaque noeud
            for erreur_i in (absolute((h_all_sui_i - h_all_all[:, t + 1]) / 
                                      h_all_all[:, t + 1])):
                print("Erreur noeud {:4d} : {:10.10e}".format(i, erreur_i))
                # Verifier erreur entre les solutions début / final < erreur_max
                if erreur_i > erreur_max:
                    erreur = True
                    print("--- Nouvelle iteration ---")
                    j = j + 1
                    print("********** Iteration : {} **********".format(j))
                    # Itérer avec le nouveau vecteur
                    h_all_all[:, t + 1] = h_all_sui_i
                    # Assambler la matrice de coefficients et le vecteur colonne
                    K, F = assSys(L, J, dt[t], h_all_all[:, t], h_all_all[:, t + 1], 
                                  condLim, condLim_0_all[t + 1], 
                                  condLim_L_all[t + 1], temp[t], p_l[t], n_l[t], 
                                  Pvsat[t], M_l, R, phi, Sr, k_l, m, D_0, A, B)
                    # Résoudre l'humidite
                    h_all_sui_i = solMatTri(K, F)
                    if condLim[0] == 0:
                        h_all_sui_i = append(condLim_0_all[t + 1], h_all_sui_i)
                    elif condLim[1] == 0:
                        h_all_sui_i = append(h_all_sui_i, condLim_L_all[t + 1])
                    break
                else:
                    i += 1
                    erreur = False
        
        # Enregistrer la solution d'humidité relative
        h_all_all[:, t + 1] = h_all_sui_i
        
    return h_all_all