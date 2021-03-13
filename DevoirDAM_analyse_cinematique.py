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

# Fonctions calculant le volume, la chaleur ainsi que leur dérivée par rapport à theta

V_theta = lambda theta : (Vc/2)*(1 - np.cos(theta) + beta - np.sqrt(beta*beta - np.sin(theta)**2)) + Vc/(tau - 1)
dVdtheta = lambda theta : (Vc/2)*(np.sin(theta) + (np.sin(theta)*np.cos(theta))/np.sqrt(beta*beta - np.sin(theta)**2))

dQdtheta = lambda Qtot, theta, thetaC, deltaThetaC : (Qtot*np.pi)*np.sin(np.pi*(theta - thetaC)/deltaThetaC)/(2*deltaThetaC)

def myfunc(rpm, s, theta, thetaC, deltaThetaC):

    DegtoRad = 2*np.pi/360                    # qui est aussi l'écart entre deux éléments du tableau thetaRadian
    thetaRadian = theta*DegtoRad
    thetaCRadian = thetaC*DegtoRad
    deltaThetaCRadian = deltaThetaC*DegtoRad

    size = theta.size

    omega = rpm*2*np.pi/60
    p_admission = s*1e5
    Qtot = Q * (p_admission * Vc) / (287.1 * 303.15)
    gamma = 1.3

    p_output = np.zeros(size)
    F_pied_output = np.zeros(size)
    F_tete_output = np.zeros(size)
    Fcrit = 0

    V_output = V_theta(thetaRadian)

    dV = dVdtheta(thetaRadian)
    dQ = dQdtheta(Qtot, thetaRadian, -thetaCRadian, deltaThetaCRadian)              # -thetaCRadian car thetaC est défini comme de la bite
    dQ[:180 - thetaC:] = 0                                                          # Apport de chaleur uniquement entre thetaC et thetaC + deltaThetaC, donc on mets à 0 les autres valeurs
    dQ[180 - thetaC + deltaThetaC:] = 0

    dPdtheta = lambda p, i: (-gamma * p * dV[i] + (gamma - 1) * dQ[i])/V_output[i]

    """Calcul de p par Euler explicite"""
    p_output[0] = p_admission
    for i in range(size - 1):
        p_output[i + 1] = p_output[i] + DegtoRad*dPdtheta(p_output[i], i)

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

        if(F_pied >= 0 >= F_tete):
            F_compression = min(F_pied, -F_tete)
            if(Fcrit < F_compression):
                Fcrit = F_compression

    print("Force critique {} [N], ps je suis dans myfunc()".format(Fcrit))

    """Calcul de t """ # TODO ça pue encore du cul, on a 6e-7 m de section, c'est beaucoup trop peu, aussi j'ai pas fait dans le sens perpendiculaire au mouvemement
    sigma = 450e6   # résistance de compression 450 MPa
    E = 200e9       # module d'élasticité 200 GPa
    Kx = 1          # facteur de correction dans le plan du mouvement
    Ky = 0.5        # facteur de correction dans le plan perpendiculaire au mouvemement

    coeffEulerX = (419*np.pi*np.pi*E)/(12*Kx*Kx*L*L)
    ax = coeffEulerX/Fcrit
    bx = coeffEulerX/11*sigma

    rootsX = np.roots([-1, 0, -bx, 0, ax])

    tx = np.amax(np.real(rootsX))

    coeffEulerY = (131 * np.pi * np.pi * E) / (12 * Ky * Ky * L * L)
    ay = coeffEulerY / Fcrit
    by = coeffEulerY / 11 * sigma

    rootsY = np.roots([-1, 0, -by, 0, ay])

    ty = np.amax(np.real(rootsY))

    t = 0

    return (V_output, dQ, F_pied_output, F_tete_output, p_output, t)


rpm = 3000
s = 1.8
theta = np.arange(-180, 181)
thetaC = 35
deltaThetaC = 70

V_output, Q_output, F_pied_output, F_tete_output, p_output, t = myfunc(rpm, s, theta, thetaC, deltaThetaC)


def beauxPlots():
    """
    beaux plots :)
    :return:
    """
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


def beauPlot():
    """
    fonction débile car les échelles ne sont pas compatibles
    :return:
    """
    plt.figure()
    plt.plot(theta, V_output, label="Volume par rapport a theta en [m^3]")

    plt.plot(theta, Q_output, label="Ajout de chaleur par rapport a theta en [J]")

    plt.plot(theta, F_pied_output, label="F_pied en [N]")
    plt.plot(theta, F_tete_output, label="F_tete en [N]")

    plt.plot(theta, p_output/1e5, label="Pression en bar par rapport a theta en [bar]")

    plt.title("Tous les beaux graphes en 1")
    plt.legend()

    print("la section t de la bielle vaut: {} [m]".format(t))

    plt.show()


# beauPlot()
beauxPlots()
