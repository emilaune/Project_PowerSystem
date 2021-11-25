from numpy import *
import numpy as np
from cmath import *
import cmath as cm


'''class Bus:
    def __init__(self,data): # data has to include following in this exact order [bus_nr, P_spec, Q_spec, v, d]
        self.number = data[0]
        self.P_specified = data[1]
        self.Q_specified = data[2]
        self.v = data[3]
        self.d = data[4]
        self.P = None
        self.Q = None

    def printBus(self):
        print(self.number)
        return None

    def loadflowsolution(self):
        print('Load flow solution in the following order: P, Q, v, d'/n)
        return self.P, self.Q, self.v, self.d'''

def obtain_Jacobian(Ybus, Vvector, Dvector):
    # writing this function because i could not figure how to create functions with variables and derive them
    # creating a zero matrix
    jacobianMatrix = np.zeros((3, 3))

    # declaring every element in the first row (P2)
    jacobianMatrix[0][0] = Vvector[1]*Vvector[0]*Ybus[1][0][0]*np.sin(Ybus[1][0][1]-Dvector[1]+Dvector[0])
    jacobianMatrix[0][1] = 0
    jacobianMatrix[0][2] = 0

    # declaring every element in second row (P3)
    jacobianMatrix[1][0] = 0
    jacobianMatrix[1][1] = Vvector[2]*Vvector[0]*Ybus[2][0][0]*np.sin(Ybus[2][0][1]-Dvector[2]+Dvector[0])
    jacobianMatrix[1][2] = Vvector[0]*Ybus[2][0][0]*np.cos(Ybus[2][0][1]-Dvector[2]+Dvector[0]) + 2*Vvector[2] \
                           * Ybus[2][2][0]*np.cos(Ybus[2][2][1])

    # declaring every element on 3rd row (Q3)
    jacobianMatrix[2][0] = 0
    jacobianMatrix[2][1] = Vvector[2]*Vvector[0]*Ybus[2][0][0]*np.cos(Ybus[2][0][1]-Dvector[2]+Dvector[0])
    jacobianMatrix[2][2] = -Vvector[0]*Ybus[2][0][0]*np.sin(Ybus[2][0][1]-Dvector[2]+Dvector[0]) - 2*Vvector[2] \
                           * Ybus[2][2][0]*np.sin(Ybus[2][2][1])

    return jacobianMatrix


def calc_powers(Ybus, V, D):
    P = np.zeros(len(V))
    Q = np.zeros(len(V))
    U = np.zeros(len(V))
    for j in range(3):
        for k in range(3):
            rad = Ybus[j][k][1]
            P[j] += V[j]*V[k]*Ybus[j][k][0]*np.cos(rad - D[j] + D[k])
            Q[j] -= V[j]*V[k]*Ybus[j][k][0]*np.sin(rad - D[j] + D[k])
    U[0], U[1], U[2] = P[1], P[2], Q[2]
    return U

def obtain_ybus(r12, x12, r13, x13):
    # First: off-diagonal elements
    z12 = complex(r12, x12)
    Y12_rect = -1 / z12
    Y12 = cm.polar(Y12_rect)
    Y21 = Y12

    z13 = complex(r13, x13)
    Y13_rect = -1 / z13
    Y13 = cm.polar(Y13_rect)
    Y31 = Y13

    Y23_rect = 0
    Y23 = cm.polar(Y23_rect)
    Y32 = Y23

    # Second: diagonal elements
    Y11_rect = 1 / z12 + 1/z13
    Y11 = cm.polar(Y11_rect)

    z22 = z12
    Y22_rect = 1 / z22
    Y22 = cm.polar(Y22_rect)

    z33 = z13
    Y33_rect = 1 / z33
    Y33 = cm.polar(Y33_rect)

    # constructing Y_bus matrix
    ybus = array([[Y11, Y12, Y13], [Y21, Y22, Y23], [Y31, Y32, Y33]])

    return ybus


def my_threebus_NR_10076(r12, x12, r13, x13, p2, v2, p3, q3):
    myangles = [0, 0, 0]
    myVs = [1, v2, 1]  # replaced mvVs[2] = 1 with v2
    # vector with p2, p3 and q3
    vectorU = np.zeros(3)
    vectorU[0], vectorU[1], vectorU[2] = p2, p3, q3
    # vector with d2, d3 and v3
    vectorX = np.zeros(3)
    vectorX[0], vectorX[1], vectorX[2] = myangles[1], myangles[2], myVs[2]

    # obtain Y-bus
    Y_BUS = obtain_ybus(r12, x12, r13, x13)
    # Load flow equations: obtaining P2_cal, P3_cal and Q3_cal

    U_cal = calc_powers(Y_BUS, myVs, myangles)
    # obtaining jacobian matrix
    jacobian = obtain_Jacobian(Y_BUS, myVs, myangles)

    # declaring delta vectors
    deltaX = np.zeros(3)
    deltaU = np.zeros(3)

    # iterating until error is small enough
    iteration = 0
    e = 0.01
    while e >= 0.01:
        deltaU = vectorU - U_cal
        deltaX = np.dot(np.linalg.inv(jacobian), deltaU)
        myangles[1] += deltaX[0]
        myangles[2] += deltaX[1]
        myVs[2] += deltaX[2]
        iteration += 1

        #print("Iteration ", iteration, ": deltaU:", deltaU, "deltaX:", deltaX)
        if abs(deltaU[0]) < e and abs(deltaU[1]) < e and abs(deltaU[2]) < e:
            break
        else:
            jacobian = obtain_Jacobian(Y_BUS, myVs, myangles)
            U_cal = calc_powers(Y_BUS, myVs, myangles)

    return myangles, myVs


if __name__ == '__main__':
    import random
    from newton_raphson import check_3bus

    np.set_printoptions(precision=3)
    np.set_printoptions(suppress=True)

    print('Testing 3-bus NR loadflow solver for candidate: 10076')

    for i in range(0, 5):
        # All values in pu
        R12 = random.randint(3, 6) / 100  # R between 0.03 and 0.06
        X12 = random.randint(20, 30) / 100  # X between 0.2 and 0.3

        R13 = random.randint(2, 4) / 100  # R between 0.02 and 0.04
        X13 = random.randint(10, 20) / 100  # X between 0.1 and 0.2

        P2 = random.randint(1, 10) / 10  # P between 0.1 and 1
        V2 = 1 + random.randint(-5, 5) / 100  # V between 0.95 and 1.05

        P3 = -random.randint(1, 10) / 10  # P between -0.1 and -1
        Q3 = -random.randint(1, 10) / 10  # Q between -0.1 and -1

        myangles, myVs = my_threebus_NR_10076(R12, X12, R13, X13, P2, V2, P3, Q3)

        print('Test {}: '.format(i + 1), end='')
        check_3bus(R12, X12, R13, X13, P2, V2, P3, Q3,myangles, myVs,set_print= False)