import numpy as np
import matplotlib.pyplot as plt

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
mbielle = 0.35     # masse de la bielle en kg
Q = 1650e3         # valeur chaleur emise par fuel par kg de melange admis (Diesel) en J
Vc = 3.3e-3/4      # cylindrée d'un piston


def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    """
    Fonction calculant pour un moteur Diesel 4 cylindres l'épaisseur critique t de la bielle
    afin qu'elle résiste au flambage

        Parameters:
        -----------
            rpm: nombre de tours par minute du moteur

            s: taux de suralimentation

            theta: angle de rotation sous forme de numpy array de -180 degrés à 180 degrés

            thetaC: angle d'allumage en degré d'angle de vilebrequin avant le PMH

            deltaThetaC: durée de combustion en degrés d'angle de vilebrequin

        Returns:
        -----------
            V_output: l'évolution du volume en [m^3] du cylindre en fonction de l'angle de rotation theta

            Q_output: l'évolution de l'apport de chaleur en [J] en fonction de l'angle de rotation theta

            F_pied_output: l'évolution de la force en [N] s'appliquant sur le pied de bielle

            F_tete_output: l'évolution de la force en [N] s'appliquant sur la tete de bielle

            p_output: l'évolution de la pression en [Pa] dans le cylindre en fonction de l'angle de rotation theta

            t: l'épaisseur critique en [m] de la bielle en forme de I
    """

    """Fonctions calculant le volume et sa dérivée par rapport à theta pour un angle theta donné"""
    V_theta = lambda theta: (Vc/2) * (1 - np.cos(theta) + beta - np.sqrt(beta*beta - np.sin(theta)**2)) + Vc/(tau-1)
    dVdtheta = lambda theta: (Vc/2) * (np.sin(theta) + (np.sin(theta)*np.cos(theta))/np.sqrt(beta*beta - np.sin(theta)**2))

    """Passage de degré en radian pour tous les paramètres le nécessitant"""
    DegtoRad = 2*np.pi/360
    thetaRadian = theta*DegtoRad
    thetaCRadian = thetaC*DegtoRad
    deltaThetaCRadian = deltaThetaC*DegtoRad

    size = theta.size

    """Détermination des constantes"""
    omega = rpm*2*np.pi/60
    p_admission = s*1e5
    Qtot = Q * (p_admission * Vc) / (287.1 * 303.15)
    gamma = 1.3

    p_output = np.zeros(size)
    F_pied_output = np.zeros(size)
    F_tete_output = np.zeros(size)

    V_output = V_theta(thetaRadian)
    dV = dVdtheta(thetaRadian)

    """Fonction calculant l'apport de chaleur par rapport à theta"""
    dQdtheta = lambda theta, thetaC, deltaThetaC: (Qtot * np.pi) * np.sin(
        np.pi * (theta - thetaC) / deltaThetaC) / (2 * deltaThetaC)

    Q_output = dQdtheta(thetaRadian, -thetaCRadian, deltaThetaCRadian)
    Q_output[:180 - thetaC:] = 0                          # Apport de chaleur uniquement entre thetaC et thetaC + deltaThetaC
    Q_output[180 - thetaC + deltaThetaC:] = 0

    dPdtheta = lambda i: (-gamma * p_output[i] * dV[i] + (gamma - 1) * Q_output[i])/V_output[i]

    """Calcul de p par Euler explicite"""
    p_output[0] = p_admission
    for i in range(size - 1):
        p_output[i + 1] = p_output[i] + DegtoRad*dPdtheta(i)

    """Calcul de F_pied_output et F_tete_output"""
    pression = (np.pi*D*D/4)*p_output
    acceleration = R*omega*omega*np.cos(thetaRadian)
    F_pied_output = pression - mpiston*acceleration
    F_tete_output = -pression + (mpiston + mbielle)*acceleration

    """Détermination de la force critique"""
    Fcompression = np.minimum(-F_tete_output, F_pied_output)
    Fcrit = np.max(Fcompression)

    """Calcul de t"""
    sigma = 450e6   # résistance de compression 450 MPa
    E = 200e9       # module d'élasticité 200 GPa
    Kx = 1          # facteur de correction dans le plan du mouvement (axe x)
    Ky = 0.5        # facteur de correction dans le plan perpendiculaire au mouvemement (axe y)
    Ixx = 419/12    # moment d'inertie dans l'axe x
    Iyy = 131/12    # moment d'inertie dans l'axe y

    coeffEuler = (np.pi*np.pi*E)/(L*L)

    ax = coeffEuler*Ixx/(Fcrit*Kx*Kx)
    bx = -coeffEuler*Ixx/(11*sigma*Kx*Kx)
    tx = np.roots([ax, 0, bx, -1])

    ay = coeffEuler*Iyy/(Fcrit*Ky*Ky)
    by = -coeffEuler*Iyy/(11*sigma*Ky*Ky)
    ty = np.roots([ay, 0, by, 0, -1])

    t = max(max(np.real(tx)), max(np.real(ty)))

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)


rpm = 1500
s = 1.8
theta = np.arange(-180, 181)
thetaC = 35
deltaThetaC = 41

V_output, Q_output, F_pied_output, F_tete_output, p_output, t = myfunc(rpm, s, theta, thetaC, deltaThetaC)


def beauxPlots():

    plt.figure()
    plt.plot(theta, V_output)
    plt.title("Volume par rapport a theta en [m^3]")

    plt.figure()
    plt.plot(theta, Q_output)
    plt.title("Ajout de chaleur par rapport a theta en [J]")

    plt.figure()
    plt.plot(theta, F_pied_output, label="F_pied")
    plt.plot(theta, F_tete_output, label="F_tete")
    plt.title("Force sur le pied et la tete de bielle en [N]")
    plt.legend()

    plt.figure()
    plt.plot(theta, p_output/1e5)
    plt.title("Pression en bar par rapport a theta en [bar]")

    print("la section t de la bielle vaut: {} [m]".format(t))

    plt.show()


beauxPlots()
