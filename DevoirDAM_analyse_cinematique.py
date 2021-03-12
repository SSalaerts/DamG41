import numpy as np

"""
Analyse cinématique et dynamique en vue de dimensionner la section d'une bielle.

Le moteur étudié est un moteur Diesel à 4 temps
largement basé sur le moteur trouvé au lien https://powersuite.cummins.com/PS5/PS5Content/SiteContent/en/Binary_Asset/pdf/KentData/DataSheets/FR30233.pdf

La longueur de la bielle n'étant pas référencée, nous avons décidé d'un ratio beta arbitraire de 3.
"""

# Grandeurs géométriques

tau = 18           # taux de compression
D = 0.095          # diamètre du piston en mètres
C = 0.115          # course du piston en mètres
R = C/2            # longueur de la manivelle obtenue grâce à la longueur de la course
beta = 3           # ratio entre la longueur de la bielle et la manivelle
L = R*beta         # longueur de la bielle que l'on cherche à dimensionner
mpiston = 0.25     # masse du piston en kg
mbielle = 0.35     # masse de la bielle
Q = 165000         # valeur chaleur emise par fuel par kg de melange admis (Diesel) en J

Vc = (3.3e-3)/4    # cylindrée d'un piston

 #énergie par kg pour le Diesel * la masse trouvée grâce au premier principe : pV = m*R#*T
                                          # (pression d'admission de 2.5 bar car presence de turbo, temperature du gaz de 30°C = 303.15 K
                                          # et R# = 287.1 pour de l'air (néglige le carburant) TODO vérifier que les approx sont légitimes

gamma = 1.3

# Fonctions calculant le volume, la chaleur ainsi que leur dérivée par rapport à theta

V_theta = lambda theta : (Vc/2)*(1 - np.cos(theta) + beta - np.sqrt(theta*theta - np.sin(theta)**2)) + Vc/(tau - 1)
dVdtheta = lambda theta : (Vc/2)*(np.sin(theta) + (np.sin(theta)*np.cos(theta))/np.sqrt(beta*beta - np.sin(theta)**2))

Q_theta = lambda Qtot, theta, thetaC, deltaThetaC : (Qtot / 2) * (1 - np.cos(np.pi * ((theta - thetaC) / deltaThetaC)))
dQdtheta = lambda Qtot, theta, thetaC, deltaThetaC : (Qtot*np.pi/2*deltaThetaC)*np.sin(np.pi*(theta-thetaC)/deltaThetaC)

# Fonction calculant la dérivée par rapport à theta de la pression

dPdtheta = lambda p, theta, thetaC, deltaThetaC : -gamma*p*dVdtheta(theta)/V_theta(theta) + (gamma - 1)*dQdtheta(theta, thetaC, deltaThetaC)/V_theta(theta)

# Euler explicite en fait
# A vérifier si c'est bien comme implémentation

def eulerExplicite(thetaC, deltaThetaC):

    P = np.zeros(360)
    P[0] = 2.5e5

    for theta in range(360):
        P[theta+1] = P[theta] + dPdtheta(P[theta], theta, thetaC, deltaThetaC)
        if(thetaC == deltaThetaC): # TODO les conditions sont pas bonnes / maybe yen a pas besoin dépendant de ma fonction
            pass
        else:
            pass

    return P


def myfunc(rpm, s, theta, thetaC, deltaThetaC):

    Qtot = Q*(s*(1e5)*Vc)/(287.1*303.15)

    omega = rpm*2*np.pi/60

    V_output = V_theta(theta)
    Q_output = Q_theta(Qtot, theta, thetaC, deltaThetaC)

    dV = dVdtheta(theta)
    dQ = dQdtheta(Qtot, theta, thetaC, deltaThetaC)

    p_output = 0
    F_pied_output = 0
    F_tete_output = 0
    t = 0

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)

