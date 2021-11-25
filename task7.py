import numpy as np
import cmath as cm

# Part B - task 7
# symmetrical fault studies
# all results used from part A task 5 are written in this code, previous tasks not imported or included here

# function to convert a vector of phase current vector in pu into a vector of its actual value [A]
def fromputoactualcurrent(puCurrentVector):
    sbase = 100 * 10**6
    vbase = 132000
    ibase = sbase / vbase
    actualCurrentVector = [[], [], []]
    actualCurrentVector[0] = puCurrentVector[0] * ibase
    actualCurrentVector[1] = puCurrentVector[1] * ibase
    actualCurrentVector[2] = puCurrentVector[2] * ibase
    return actualCurrentVector

# function to convert from pu to actual value [kV]
def fromputoactualvoltage(puVoltageVector):
    vbase = 132000
    c = complex(0, 0)
    actualVoltageVector = np.array([[c, c, c], [c, c, c], [c, c, c], [c, c, c], [c, c, c]])
    for i in range(5):
        for l in range(3):
            actualVoltageVector[i][l] = puVoltageVector[i][l] * vbase / 1000
    return actualVoltageVector


# question a-------------------------------------------------------------------
# obtain z_bus

def positivesequencebus():
    # off-diagonal elements
    Y12 = complex(-0.769, 3.846)  # -1 / z12_1
    Y21 = Y12
    Y13 = complex(0, 0)
    Y31 = Y13
    Y14 = complex(-0.495, 4.95)
    Y41 = Y14
    Y15 = 0
    Y51 = Y15
    Y23 = complex(-1.538, 7.692)
    Y32 = Y23
    Y24 = complex(-0.495, 4.950)
    Y42 = Y24
    Y25 = complex(0, 0)
    Y52 = Y25
    Y34 = complex(0, 0)
    Y43 = Y34
    Y35 = complex(0, 0)
    Y53 = Y35
    Y45 = complex(0.999, 9.999)
    Y54 = Y45

    # Second: diagonal elements
    Y11 = complex(10.495, -18.95)  # (1 / zg) + (1 / (0.05+j*0.25)) + (1 / (0.02+j*0.2))
    Y22 = complex(2.802, -16.489)
    Y33 = complex(1.538, -11.692)
    Y44 = complex(1.98, -18.802)
    Y55 = complex(0.99, -9.901)

    ybus = np.array([[Y11, Y12, Y13, Y14, Y15], [Y21, Y22, Y23, Y24, Y25], [Y31, Y32, Y33, Y34, Y35], [Y41, Y42, Y43, Y44,Y45], [Y51, Y52, Y53, Y54, Y55]])
    zbus = np.linalg.inv(ybus)
    return zbus

#  question b-------------------------------------------------------------------
#  find the fault currents at bus 5


def symmetricalfaultcurrent(zbus):
    v5 = cm.rect(0.93408, -0.263)  # from results in part A task 5
    Ifa = v5 / zbus[4][4]
    Ifb = Ifa * (cm.rect(1, (2 * np.pi) / 3))**2
    Ifc = Ifa * (cm.rect(1, (2 * np.pi) / 3))
    print("Symmetrical currents in [PU]:")
    print("IFa: ", Ifa, "IFb: ", Ifb, "IFc: ", Ifc)
    print("\nSymmetrical currents in actual values [A]:")
    actualValues = fromputoactualcurrent([Ifa, Ifb, Ifc])
    print("IFa: ", actualValues[0], "IFb: ", actualValues[1], "IFc: ", actualValues[2])
    return Ifa, Ifb, Ifc


#  question c-------------------------------------------------------------------
#  voltages at all buses
def busvoltages(zbus):
    # pre-fault values in pu
    v1pre = cm.rect(1, -0.0631)
    v2pre = cm.rect(0.95707, -0.115)
    v3pre = cm.rect(1, 0)
    v4pre = cm.rect(0.98342, -0.212)
    v5pre = cm.rect(0.93408, -0.263)
    preVoltages = [v1pre, v2pre, v3pre, v4pre, v5pre]

    # obtaining delta voltages
    a = cm.rect(1, (2 * np.pi / 3))
    matrixA = np.array([[1, 1, 1], [1, a**2, a], [1, a, a**2]])  # [A]
    c = complex(0, 0)  # 0 + 0*j
    deltaVoltages = np.array([[c, c, c], [c, c, c], [c, c, c], [c, c, c], [c, c, c]])
    postVoltages = np.array([[c, c, c], [c, c, c], [c, c, c], [c, c, c], [c, c, c]])
    postPhaseVoltages = np.array([[c, c, c], [c, c, c], [c, c, c], [c, c, c], [c, c, c]])
    for i in range(5):
        b = (-zbus[i][4] / zbus[4][4]) * v5pre
        deltaVoltages[i][1] = b
        postVoltages[i][1] = preVoltages[i] + deltaVoltages[i][1]  # V_post = V_pre + deltaV
        postPhaseVoltages[i] = np.dot(matrixA, postVoltages[i])  # [Vabc] = [A] * [V012]
    print("\nSymmetrical post-fault phase voltages in [PU]:\n", postPhaseVoltages)
    actualVoltages = fromputoactualvoltage(postPhaseVoltages)
    print("\nSymmetrical post-fault phase voltages in actual values [kV]:\n", actualVoltages)
    return postPhaseVoltages, postVoltages

#  question d -------------------------------------------------------------------
#  find the current flowing from bus 4 to bus 5

def currentfrom4to5(postVoltages):
    c = complex(0, 0)  # 0 + 0*j
    z45 = complex(0.01, 0.1)
    a = cm.rect(1, (2 * np.pi / 3))
    matrixA = np.array([[1, 1, 1], [1, a ** 2, a], [1, a, a ** 2]])  # [A]
    currentVector = [c, c, c]
    currentVector[1] = (postVoltages[3][1]-postVoltages[4][1]) / z45
    currentPhaseVector = np.dot(matrixA, currentVector)
    print("\nSymmetrical phase currents from bus 4 to bus 5 in [PU]:")
    print("I45a: ", currentPhaseVector[0], "I45b: ", currentPhaseVector[1], "I45c: ", currentPhaseVector[2])
    actualCurrents = fromputoactualcurrent(currentPhaseVector)
    print("\nSymmetrical phase currents from bus 4 to bus 5 in actual value [A]: ", currentPhaseVector)
    print("I45a: ", actualCurrents[0], "I45b: ", actualCurrents[1], "I45c: ", actualCurrents[2])
    return currentPhaseVector

