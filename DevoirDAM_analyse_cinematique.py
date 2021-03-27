import numpy as np
import matplotlib.pyplot as plt
import time
"""
Analyse cinématique et dynamique en vue de dimensionner la section d'une bielle.

Le moteur étudié est un moteur Diesel à 4 temps
largement basé sur le moteur trouvé au lien https://powersuite.cummins.com/PS5/PS5Content/SiteContent/en/Binary_Asset/pdf/KentData/DataSheets/FR30233.pdf

La longueur de la bielle n'étant pas référencée, nous avons décidé d'un ratio beta arbitraire de 3.
"""

# Grandeurs géométriques

tau = 18  # [-]
D = 0.095  # [m]
C = 0.115  # [m]
R = C/2  # [m] longueur de la manivelle
beta = 3  # [-] ratio entre la longueur de la bielle et celle de la manivelle
L = beta*R  # [m] longueur de la bielle
mpiston = 0.25  # [kg]
mbielle = 0.35  # [kg]
Q = 1650e3  # [J/kg_inlet gas]
Vc = 3.3e-3/4  # [m^3] cylindrée d'un piston
print(Vc)


# D = 0.1 #diamètre piston [m]
# L = 0.15 #longueur bielle [m]
# tau = 20#taux de compression [-]
# C = 0.1 #longueur course [m]
# R = C/2
# mpiston = 0.25#valeur masse piston [kg]
# mbielle = 0.35 #masse bielle [kg]
# Q = 1650000 #[J/kg_inlet gas]
# Vc = (np.pi*D**2)/4*C
# print(Vc)


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

    #=== Détermination des constantes ===#
    size = theta.size
    PI = np.pi
    omega = rpm * 2 * PI / 60
    p_admission = s * 1e5
    Qtot = Q * (p_admission * Vc) / (287.1 * 303.15)
    gamma = 1.3
    beta = 3

    #=== Fonctions calculant le volume et sa dérivée par rapport à theta pour un angle theta donné ===#
    V_theta = lambda theta: (Vc/2) * (1 - np.cos(theta) + beta - np.sqrt(beta*beta - np.sin(theta)**2)) + Vc/(tau-1)
    dVdtheta = lambda theta: (Vc/2) * (np.sin(theta) + (np.sin(theta)*np.cos(theta))/np.sqrt(beta*beta - np.sin(theta)**2))

    #=== Passage de degré en radian pour tous les paramètres le nécessitant ===#
    DegtoRad = PI/180
    thetaRadian = theta*DegtoRad
    thetaCRadian = thetaC*DegtoRad
    deltaThetaCRadian = deltaThetaC*DegtoRad

    #=== Calcul de V_output et de la variation de volume par rapport à theta ===#
    V_output = V_theta(thetaRadian)
    dV = dVdtheta(thetaRadian)

    #=== Fonction calculant l'apport de chaleur par rapport à theta ===#
    dQdtheta = lambda theta, thetaC, deltaThetaC: (Qtot * PI) * np.sin(PI * (theta - thetaC) / deltaThetaC) / (2 * deltaThetaC)
    Q_theta = lambda theta, thetaC, deltaThetaC: (Qtot/2)*(1 - np.cos(PI*(theta - thetaC)/deltaThetaC))
    dQ = dQdtheta(thetaRadian, -thetaCRadian, deltaThetaCRadian)
    Q_output = Q_theta(thetaRadian, -thetaCRadian, deltaThetaCRadian)

    count = 0
    for i in range(size):
        if theta[i] < -thetaC or theta[i] > -thetaC + deltaThetaC:
            dQ[i] = 0
            Q_output[i] = 0

        if theta[i] > -thetaC + deltaThetaC:
            count += 1
            dQ[i] = 0
            Q_output[i] = Q_output[i - count]

    #=== Calcul de p par Euler explicite ===#
    p_output = np.zeros(size)
    dPdtheta = lambda i: (-gamma * p_output[i] * dV[i] + (gamma - 1) * dQ[i])/V_output[i]
    p_output[0] = p_admission
    h = (theta[-1] - theta[0])/(size-1)

    for i in range(size - 1):
        p_output[i + 1] = p_output[i] + h*DegtoRad*dPdtheta(i)

    #=== Calcul de F_pied_output et F_tete_output ===#
    F_pression = (PI*D*D/4)*p_output
    acceleration = R*omega*omega*np.cos(thetaRadian)
    F_pied_output = F_pression - mpiston*acceleration
    F_tete_output = -F_pression + (mpiston + mbielle)*acceleration

    #=== Détermination de la force critique ===#
    Fcompression = np.minimum(-F_tete_output, F_pied_output)
    Fcrit = np.max(Fcompression)

    #=== Calcul de t ===#
    sigma = 450e6   # résistance de compression à 450 MPa
    E = 200e9       # module d'élasticité à 200 GPa
    Kx = 1          # facteur de correction dans le plan du mouvement (axe x)
    Ky = 0.5        # facteur de correction dans le plan perpendiculaire au mouvemement (axe y)
    Ixx = 419/12    # coefficient du moment d'inertie dans l'axe x
    Iyy = 131/12    # coefficient du moment d'inertie dans l'axe y

    coeffEuler = (PI*PI*E)/(L*L)

    ax = coeffEuler*Ixx/(Fcrit*Kx*Kx)
    bx = -coeffEuler*Ixx/(11*sigma*Kx*Kx)
    tx = np.roots([ax, 0, bx, -1])

    ay = coeffEuler*Iyy/(Fcrit*Ky*Ky)
    by = -coeffEuler*Iyy/(11*sigma*Ky*Ky)
    ty = np.roots([ay, 0, by, 0, -1])

    t = max(max(np.real(tx)), max(np.real(ty)))

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)


rpm = 4467
s = 0.8
theta = np.linspace(-180, 180, 1001)
thetaC = 40.0
deltaThetaC = 58.0


t1 = time.perf_counter()
V_output, Q_output, F_pied_output, F_tete_output, p_output, t = myfunc(rpm, s, theta, thetaC, deltaThetaC)
t2 = time.perf_counter()
print("time taken =", t2-t1, "[s]")


def beauxPlots():

    # print(V_output)
    #
    # print(Q_output)
    #
    # print(F_pied_output)
    #
    # print(F_tete_output)
    #
    # print(p_output)
    #
    # print(t)

    plt.figure()
    plt.plot(theta, V_output)
    plt.title("Volume par rapport a theta en [m^3]")

    plt.figure()
    plt.plot(theta, Q_output)
    plt.title("Apport de chaleur par rapport a theta en [J]")

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

