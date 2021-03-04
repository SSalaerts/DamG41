import numpy as np

"""
Analyse cinématique et dynamique en vue de dimensionner la section d'une bielle.

Le moteur étudié est un moteur Diesel à 4 temps
largement basé sur le moteur trouvé au lien https://powersuite.cummins.com/PS5/PS5Content/SiteContent/en/Binary_Asset/pdf/KentData/DataSheets/FR30233.pdf

La longueur de la bielle n'étant pas référencée, nous avons décidé d'un ratio beta arbitraire de 3.
"""

# Grandeurs géométriques

D = 0.095          # diamètre du piston en mètres
C = 0.115          # course du piston en mètres
Vc = (3.3E-3)/4    # cylindrée d'un piston
R = C/2            # longueur de la manivelle obtenue grâce à la longueur de la course
tau = 18           # taux de compression
beta = 3           # ratio entre la longueur de la bielle et la manivelle
L = R*beta         # longueur de la bielle que l'on cherche à dimensionner
gamma = 1.3

m_piston = 0.25 # masse du piston en kg
m_bielle = 0.35 # masse de la bielle

Qtot = 165000*((2.5E5)*Vc)/(287.1*303.15) #énergie par kg pour le Diesel * la masse trouvée grâce au premier principe : pV = m*R#*T
                                          # (pression d'admission de 2.5 bar car presence de turbo, temperature du gaz de 30°C = 303.15 K
                                          # et R# = 287.1 pour de l'air (néglige le carburant) TODO vérifier que les approx sont légitimes

# Fonctions calculant le volume, la chaleur ainsi que leur dérivée par rapport à theta

V_theta = lambda theta : (Vc/2)*(1 - np.cos(theta) + beta - np.sqrt(theta*theta - np.sin(theta)**2)) + Vc/(tau - 1)
dVdtheta = lambda theta : (Vc/2)*(np.sin(theta) + (np.sin(theta)*np.cos(theta))/np.sqrt(beta*beta - np.sin(theta)**2))

Q_theta = lambda theta, theta_d, theta_comb : (Qtot/2)*(1 - np.cos(np.pi*((theta - theta_d)/theta_comb)))
dQdtheta = lambda theta, theta_d, theta_comb : (Qtot*np.pi/2*theta_comb)*np.sin(np.pi*(theta-theta_d)/theta_comb)

# Fonction calculant la dérivée par rapport à theta de la pression

dPdtheta = lambda p, theta, theta_d, theta_comb : -gamma*p*dVdtheta(theta)/V_theta(theta) + (gamma - 1)*dQdtheta(theta, theta_d, theta_comb)/V_theta(theta)

# Euler explicite en fait
# A vérifier si c'est bien comme implémentation

def eulerExplicite(P0, theta_d, theta_comb):

    P = np.zeros(360)
    P[0] = P0

    for theta in range(360):
        P[theta+1] = P[theta] + dPdtheta(P[theta], theta, theta_d, theta_comb)
        if(theta_d == theta_comb): # TODO les conditions sont pas bonnes encore fin ça dépends de la définition de la fonction donc maybe ya pas besoin
            pass
        else:
            pass

    return P