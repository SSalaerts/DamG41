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
Vc = 3.3e-3/4    # cylindrée d'un piston

# Fonctions calculant le volume, la chaleur ainsi que leur dérivée par rapport à theta

V_theta = lambda theta : (Vc/2)*(1 - np.cos(theta) + beta - np.sqrt(beta*beta - np.sin(theta)**2)) + Vc/(tau - 1)
dVdtheta = lambda theta : (Vc/2)*(np.sin(theta) + (np.sin(theta)*np.cos(theta))/np.sqrt(beta*beta - np.sin(theta)**2))

Q_theta = lambda Qtot, theta, thetaC, deltaThetaC : (Qtot / 2) * (1 - np.cos(np.pi * ((theta - thetaC) / deltaThetaC)))
dQdtheta = lambda Qtot, theta, thetaC, deltaThetaC : (Qtot*np.pi/2*deltaThetaC)*np.sin(np.pi*(theta-thetaC)/deltaThetaC)


# Fonction calculant la dérivée par rapport à theta de la pression
def myfunc(rpm, s, theta, thetaC, deltaThetaC):

    DegtoRad = 2*np.pi/360
    thetaRadian = theta*DegtoRad
    thetaCRadian = thetaC*DegtoRad
    deltaThetaCRadian = deltaThetaC*DegtoRad

    h = np.abs(thetaRadian[0] - thetaRadian[1])
    omega = rpm*2*np.pi/60
    p_admission = s*1e5
    Qtot = Q * (p_admission * Vc) / (287.1 * 303.15)
    gamma = 1.3

    size = theta.size
    p_output = np.zeros(size)
    F_pied_output = np.zeros(size)
    F_tete_output = np.zeros(size)
    Fcrit = 0

    V_output = V_theta(thetaRadian)
    Q_output = Q_theta(Qtot, thetaRadian, thetaCRadian, deltaThetaCRadian)

    dV = dVdtheta(thetaRadian)
    dQ = dQdtheta(Qtot, thetaRadian, thetaCRadian, deltaThetaCRadian)

    dPdtheta = lambda p, i: -gamma * p * dV[i] / V_output[i] + (
                gamma - 1) * dQ[i] / V_output[i]

    """Calcul de p par Euler explicite"""
    p_output[0] = p_admission
    for i in range(size - 1):
        p_output[i + 1] = p_output[i] + h*dPdtheta(p_output[i], i)

    """Calcul de F_pied_output et F_tete_output"""
    for i in range(size):
        pression = (np.pi*D*D/4)*p_output[i]
        acceleration = R*omega*omega*np.cos(thetaRadian[i])
        F_pied_output[i] = pression - mpiston*acceleration
        F_tete_output[i] = -pression + (mpiston + mbielle)*acceleration

    """Calcul de la force critique"""

    for i in range(size):
        F_pied = F_pied_output[i]
        F_tete = F_tete_output[i]

        if(F_pied >= 0 and F_tete <= 0):
            F_compression = min(F_pied, -F_tete)
            if(Fcrit < F_compression):
                Fcrit = F_compression


    """Calcul de t """ #TODO ça pue encore du cul, on a 2.67e17 m de section, c'est beaucoup trop, aussi j'ai pas fait dans le sens perpendiculaire au mouvemement
    sigma = 450e6   # résistance de compression 450 MPa
    E = 200e9       # module d'élasticité 200 GPa
    Kx = 1          # facteur de correction dans le plan du mouvement

    coeffEuler = (419*np.pi*np.pi*E)/(12*Kx*Kx*L*L)
    a = coeffEuler/Fcrit
    b = coeffEuler/11*sigma

    racineDelta = np.sqrt(b*b + 4*a)

    t1 = np.sqrt((b - racineDelta)/2*a)
    t2 = np.sqrt((b + racineDelta)/2*a)

    t = max(t1, t2)

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)


rpm = 3000
s = 1
theta = np.arange(-180, 181)
thetaC = 40
deltaThetaC = 70

V_output, Q_output, F_pied_output, F_tete_output, p_output, t = myfunc(rpm, s, theta, thetaC, deltaThetaC)


def beauPlot():
    plt.figure()
    plt.plot(theta, V_output)
    plt.title("Volume par rapport a theta en [m^3]")

    plt.figure()
    plt.plot(theta, Q_output)
    plt.title("Chaleur par rapport a theta en [J]")

    plt.figure()
    plt.plot(theta, F_pied_output, label="F_pied")
    plt.plot(theta, F_tete_output, label="F_tete")
    plt.title("Force sur le pied et la tete de bielle en [J]")
    plt.legend()

    plt.figure()
    plt.plot(theta, p_output)
    plt.title("Pression en bar par rapport a theta en [Pa]")

    print("la section t de la bielle vaut: {} en [m]".format(t))

    plt.show()


beauPlot()
