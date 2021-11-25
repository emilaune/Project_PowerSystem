import numpy as np
import cmath as cm

# Part B - task 7
# symmetrical fault studies
# all results used from part A task 5 are written in this code, previous tasks not imported or included here

def fromputoactualcurrent(puCurrentVector):
    sbase = 100 * 10**6
    vbase = 132000
    ibase = sbase / vbase
    actualCurrentVector = [[], [], []]
    actualCurrentVector[0] = puCurrentVector[0] * ibase
    actualCurrentVector[1] = puCurrentVector[1] * ibase
    actualCurrentVector[2] = puCurrentVector[2] * ibase

    return actualCurrentVector

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

def negativesequencebus():
    # off-diagonal elements
    Y12 = complex(-0.769, 3.846)
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
    Y11 = complex(1.264, -15.463)
    Y22 = complex(2.802, -16.489)
    Y33 = complex(1.538, -14.359)
    Y44 = complex(1.98, -18.802)
    Y55 = complex(0.99, -9.901)

    ybus = np.array([[Y11, Y12, Y13, Y14, Y15], [Y21, Y22, Y23, Y24, Y25], [Y31, Y32, Y33, Y34, Y35], [Y41, Y42, Y43, Y44,Y45], [Y51, Y52, Y53, Y54, Y55]])
    zbus = np.linalg.inv(ybus)
    return zbus

def zerosequencebus():
    # off-diagonal elements
    Y12 = complex(-0.102, 0.898)
    Y21 = Y12
    Y13 = complex(0, 0)
    Y31 = Y13
    Y14 = complex(-0.064, 1.133)
    Y41 = Y14
    Y15 = 0
    Y51 = Y15
    Y23 = complex(-0.204, 1.795)
    Y32 = Y23
    Y24 = complex(-0.064, 1.133)
    Y42 = Y24
    Y25 = complex(0, 0)
    Y52 = Y25
    Y34 = complex(0, 0)
    Y43 = Y34
    Y35 = complex(0, 0)
    Y53 = Y35
    Y45 = complex(-0.129, 2.265)
    Y54 = Y45

    # Second: diagonal elements
    Y11 = complex(0.166, -3.312)
    Y22 = complex(0.370, -3.825)
    Y33 = complex(0.204, -15.128)
    Y44 = complex(0.257, -4.531)
    Y55 = complex(0.129, -2.265)

    ybus = np.array([[Y11, Y12, Y13, Y14, Y15], [Y21, Y22, Y23, Y24, Y25], [Y31, Y32, Y33, Y34, Y35], [Y41, Y42, Y43, Y44,Y45], [Y51, Y52, Y53, Y54, Y55]])
    zbus = np.linalg.inv(ybus)
    return zbus

#  question b-------------------------------------------------------------------
#  find the unsymmetrical fault currents at bus 5 for all three phases


def unsymmetricalfaultcurrent(zeroZbus, posZbus, negZbus):
    v5 = cm.rect(0.93408, -0.263)  # from results in part A task 5
    # declaring matrix A
    a = cm.rect(1, (2 * np.pi / 3))
    matrixA = np.array([[1, 1, 1], [1, a ** 2, a], [1, a, a ** 2]])  # [A]
    # obtaining zero-, positive- and negative-sequence fault currents
    If0 = v5 / (zeroZbus[4][4] + + posZbus[4][4] + negZbus[4][4])
    If1 = If0
    If2 = If1
    sequenceCurrentVector = [[If0], [If1], [If2]]
    # obtaining the fault current in the three phases a, b and c
    phaseCurrents = np.dot(matrixA, sequenceCurrentVector)
    Ifa = phaseCurrents[0]
    Ifb = phaseCurrents[1]
    Ifc = phaseCurrents[2]
    print("Unsymmetrical fault currents in [PU]:")
    print("IFa: ", Ifa, "IFb: ", Ifb, "IFc: ", Ifc)
    print("\nUnsymmetrical fault currents in actual values [A]:")
    actualValues = fromputoactualcurrent([Ifa, Ifb, Ifc])
    print("IFa: ", actualValues[0], "IFb: ", actualValues[1], "IFc: ", actualValues[2])
    return If0


#  question c-------------------------------------------------------------------
#  obtain the phase voltages at all buses
def unsymmetricalbusvoltages(zeroZbus, posZbus, negZbus, unsSeqCurrent):
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
    # have to declare empty matrices with complex type, could not find any other way to do it, hard coded...
    # deltaVs = np.zeros((5, 3)) would not work
    deltaVoltages = np.array([[c, c, c], [c, c, c], [c, c, c], [c, c, c], [c, c, c]])  # all sequence voltages
    postVoltages = np.array([[c, c, c], [c, c, c], [c, c, c], [c, c, c], [c, c, c]])  # all sequence voltages
    postPhaseVoltages = np.array([[c, c, c], [c, c, c], [c, c, c], [c, c, c], [c, c, c]])  # all phase voltages
    for i in range(5):
        dV0 = (-zeroZbus[i][4] * unsSeqCurrent) # dV0 = -Zik0 * If0
        dV1 = (-posZbus[i][4] * unsSeqCurrent) # dV1 = -Zik1 * If1
        dV2 = (-negZbus[i][4] * unsSeqCurrent) # dV2 = -Zik2 * If2
        # inserting delta voltages to the delta vector - this makes it a matrix with all the sequence delta voltages
        deltaVoltages[i][0] = dV0
        deltaVoltages[i][1] = dV1
        deltaVoltages[i][2] = dV2
        # postV = preV + deltaV
        postVoltages[i][0] = deltaVoltages[i][0]  # postV0 = 0 + deltaV0
        postVoltages[i][1] = preVoltages[i] + deltaVoltages[i][1]  # postV1 = preV + deltaV1
        postVoltages[i][2] = deltaVoltages[i][2]  # postV2 = 0 + deltaV2
        # [Vabc] = [A] * [V012]
        postPhaseVoltages[i] = np.dot(matrixA, postVoltages[i])
    print("\nUnsymmetrical post-fault phase voltages in [PU]:\n", postPhaseVoltages)
    actualVoltages = fromputoactualvoltage(postPhaseVoltages)
    print("\nUnsymmetrical post-fault phase voltages in actual values [kV]:\n", actualVoltages)
    return postPhaseVoltages, postVoltages

#  question d -------------------------------------------------------------------
#  find the current flowing from bus 4 to bus 5 in all three phases

def unsymmetricalcurrentfrom4to5(postPhaseVoltages):
    c = complex(0, 0)  # 0 + 0*j
    z450 = complex(0.025, 0.44)
    z451 = complex(0.01, 0.1)
    z452 = z451
    a = cm.rect(1, (2 * np.pi / 3))
    matrixA = np.array([[1, 1, 1], [1, a ** 2, a], [1, a, a ** 2]])  # [A]
    currentVector = [c, c, c]  # declaring current vector with complex type

    currentVector[0] = (postPhaseVoltages[3][0] - postPhaseVoltages[4][0]) / z450
    currentVector[1] = (postPhaseVoltages[3][1]-postPhaseVoltages[4][1]) / z451
    currentVector[2] = (postPhaseVoltages[3][2] - postPhaseVoltages[4][2]) / z452

    currentPhaseVector = np.dot(matrixA, currentVector)
    print("\nUnsymmetrical phase currents from bus 4 to bus 5 in [PU]:")
    print("I45a: ", currentPhaseVector[0], "I45b: ", currentPhaseVector[1], "I45c: ", currentPhaseVector[2])
    actualCurrents = fromputoactualcurrent(currentPhaseVector)
    print("\nUnsymmetrical phase currents from bus 4 to bus 5 in actual value [A]: ", currentPhaseVector)
    print("I45a: ", actualCurrents[0], "I45b: ", actualCurrents[1], "I45c: ", actualCurrents[2])
    return currentPhaseVector