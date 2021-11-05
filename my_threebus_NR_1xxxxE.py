# Newton-Raphson Algorithm for solving a n-bus Load Flow Problem
# TET4115 Power System Analysis; semester: Fall 2021
# Author: Eivind Falk

import numpy as np
import cmath

def my_threebus_NR_1xxxx(r12, x12, r13, x13, p2, v2, p3, q3):

    linedata = [[1,2,r12,x12],[1,3,r13,x13]]                              # [bus from, to bus, resistance, reactance]
    busdata = [[0,0,0,0,1,0],[1,1,p2,0,v2,0],[2,2,p3,q3,1,0]]             # [bus no., bus type, P, Q, |V|, delta]
    YBus = findYBus(linedata)

    myangles = [0, 0, 0]            # all unknown bus voltage phase angles are initially set to 0 rad
    myVs = [1, v2, 1]               # all unknown bus voltage magnitudes are initially set to 1.0 pu

    [pMismatch,qMismatch,nGen,nLoad,b] = isConverged(myangles, myVs, YBus,busdata)
    J = findJacobian(YBus, busdata,nGen,nLoad)
    inverse = np.linalg.inv(J)
    [myangles,myVs] = updateStateVariables(inverse,myangles,myVs,pMismatch,qMismatch,busdata,nGen,nLoad)

    return myangles, myVs

# isConverged calculates power injections at each bus by Load Flow Equations,
# finds the mismatches vector and returns the mismatches of active- and reactive powers
# of iteration "i" to be used in the iterations. Also, a boolean value is returned to indicate convergence
# or not convergence (i.e. perform another iteration)
def isConverged(myangles, myVs, YBus,busdata):
    nGen = 0                                # number of generator buses in system
    nLoad = 0                               # number of load buses in system
    for i in range(0, len(busdata)):
        if (busdata[i][1] == 1):
            nGen += 1
        elif (busdata[i][1] == 2):
            nLoad += 1
    V = myVs
    d = myangles
    P = np.zeros([nGen+nLoad+1])        # active power injections vector
    Q = np.zeros([nGen+nLoad+1])        # reactive power injections vector
    pMismatch = []                      # active power mismatch vector
    qMismatch = []                      # reactive power mismatch vector
    for i in range(0, nGen+nLoad+1):
        for j in range(0,nGen+nLoad+1):
            [A,angle] = cmath.polar(YBus[i][j])
            if (A != 0):
                P[i] = P[i]+ V[i]*V[j]*A*np.cos(d[i] - d[j] - angle)
                Q[i] = Q[i] + V[i]*V[j]*A*np.sin(d[i] - d[j] - angle)
        if (busdata[i][1]==1)or(busdata[i][1]==2):
            pMismatch.append(busdata[i][2] - P[i])
            if (busdata[i][1]==2):
                qMismatch.append(busdata[i][3] - Q[i])
    b = True
    epsilon = 0.1
    powersMismatch = np.vstack(np.hstack((pMismatch, qMismatch)))
    for i in range(0,len(powersMismatch)):
        if(np.absolute(powersMismatch[i])>epsilon):
            b = False
            break
    return pMismatch,qMismatch,nGen,nLoad,b




# findYBus is a function that takes "linedata" as input and
# returns the Admittance Bus Matrix of the system.
def findYBus(linedata):
    nbus = max(linedata[0])             # no. of network nodes
    branches = len(linedata)            # no. of network branches
    for i in range(0, len(linedata)):
        if (linedata[i][0] > nbus):
            nbus = linedata[i][0]
        elif (linedata[i][1] > nbus):
            nbus = linedata[i][1]
    y = []
    for i in range(0, branches):
        y.append(1 / complex(linedata[i][2], linedata[i][3]))

    YBus = np.zeros([nbus, nbus], dtype=complex)
    # Off-diagonal elements of YBus
    for i in range(0, branches):
        YBus[linedata[i][0] - 1][linedata[i][1] - 1] = -y[i]
        YBus[linedata[i][1] - 1][linedata[i][0] - 1] = YBus[linedata[i][0] - 1][linedata[i][1] - 1]

    # Diagonal elements of YBus
    for i in range(1, nbus + 1):
        for j in range(0, branches):
            if (linedata[j][0] == (i)) or (linedata[j][1] == i):
                YBus[i - 1][i - 1] = YBus[i - 1][i - 1] + y[j]
    #print("YBus:")
    #print("--------------------------------------------------------------")
    #for i in range(0, nbus):
    #     if (i > 0):
    #        print()
    #     for j in range(0, nbus):
    #        admittanceprint = "{:.3f} {:.3f}i".format(YBus[i][j].real, YBus[i][j].imag)
    #        print(admittanceprint, end="\t\t"),
    #print("\n--------------------------------------------------------------")
    return YBus


# findJacobian is the longest and most complex function used in the main program.
# Complexity comes from the calculation and inserting of entry values in the submatrices of the jacobian,
# when generator - and load buses are given at arbitrary positions in "linedata"
def findJacobian(YBus,busdata,nGen,nLoad):

    J1 = np.zeros([nGen+nLoad,nGen+nLoad])
    J2 = np.zeros([nGen+nLoad,nLoad])
    J3 = np.zeros([nLoad,nGen+nLoad])
    J4 = np.zeros([nLoad,nLoad])

    j1count = [0,0]             # (#P x #delta)
    j2count = [0,0]             # (#P x #V)
    j3count = [0,0]             # (#Q x #delta)
    j4count = [0,0]             # (#Q x #V)

    # Calculate the (P-delta) - sensitivities
    for i in range(1,nGen+nLoad+1):
        j1count[1] = 0
        for j in range(1, nGen+nLoad+1):
            [A,angle] = cmath.polar(YBus[i][j])
            if(i!=j)and(A==0):
                J1[j1count[0],j1count[1]] = 0
                j1count[1] += 1
            elif(i!=j)and(A!=0):
                J1[j2count[0],j1count[1]] = busdata[i][4]*busdata[j][4]*A*np.sin(busdata[i][5]-busdata[j][5]-angle)
                j1count[1] += 1
            elif(i==j):
                productP1 = 0
                for k in range(0,nGen+nLoad+1):
                    [A,angle] = cmath.polar(YBus[i][k])
                    if(i!=k)and(A!=0):
                        productP1 = productP1 + A*busdata[k][4]*np.sin(busdata[i][5]-busdata[k][5]-angle)
                J1[j1count[0],j1count[1]] = -1*busdata[i][4]*productP1
                j1count[1] += 1
        j1count[0] += 1

    # Calculates the (P-|V|) - sensitivities
    for i in range(1,nGen+nLoad+1):
        [Ai,anglei] = cmath.polar(YBus[i][i])
        j2count[1] = 0
        for j in range(1,nGen+nLoad+1):
            if(busdata[j][1]==2):
                [A,angle] = cmath.polar(YBus[i][j])
                if (i!=j)and(A==0):
                    J2[j2count[0],j2count[1]] = 0
                    j2count[1] += 1
                elif (i!=j) and (A!=0):
                    J2[j2count[0],j2count[1]] = busdata[i][4]*A*np.cos(busdata[i][5]-busdata[j][5]-angle)
                    j2count[1] += 1
                elif (i==j):
                    productP2 = 0
                    for k in range(0,nGen+nLoad+1):
                        [Ak,anglek] = cmath.polar(YBus[i][k])
                        if(Ak!=0):
                            productP2 = productP2 + Ak*busdata[k][4]*np.cos(busdata[i][5]-busdata[k][5]-anglek)
                    J2[j2count[0],j2count[1]] = busdata[i][4]*Ai*np.cos(anglei) + productP2
                    j2count[1] += 1
        j2count[0] += 1

    # Calculates the (Q-delta) - sensitivities
    for i in range(1,nGen+nLoad+1):
        if (busdata[i][1]==2):
            j3count[1] = 0
            for j in range(1,nGen+nLoad+1):
                [A,angle] = cmath.polar(YBus[i][j])
                if(i!=j)and(A==0):
                    J3[j3count[0],j3count[1]] = 0
                    j3count[1] += 1
                elif(i!=j)and(A!=0):
                    J3[j3count[0],j3count[1]] = -busdata[i][4]*A*np.cos(busdata[i][5]-busdata[j][5]-angle)
                    j3count[1] += 1
                elif(i==j):
                    productQ1 = 0
                    for k in range(0,nGen+nLoad+1):
                        [Ak,anglek] = cmath.polar(YBus[i][k])
                        if(i!=k)and(Ak!=0):
                            productQ1 += productQ1 + Ak*busdata[i][4]*np.cos(busdata[i][5]-busdata[k][5]-anglek)
                    J3[j3count[0],j3count[1]] = busdata[i][4]*productQ1
                    j3count[1] += 1
            j3count[0] += 1

    # Calculates the (Q-|V|) - sensitivities
    for i in range(1,nGen+nLoad+1):
        if (busdata[i][1]==2):
            j4count[1] = 0
            for j in range(1,nGen+nLoad+1):
                if(busdata[j][1]==2):
                    [A,angle] = cmath.polar(YBus[i][j])
                    if(i!=j)and(A==0):
                        J4[j4ccount[0],j4count[1]] = 0
                        j4count[1] += 1
                    elif(i!=j)and(A!=0):
                        J4[j4count[0],j4count[1]] = busdata[i][4]*A*np.sin(busdata[i][5]-busdata[j][5]-angle)
                        j4count[1] += 1
                    elif(i==j):
                        productQ2 = 0
                        for k in range(0,nGen+nLoad+1):
                            [Ak,anglek] = cmath.polar(YBus[i][k])
                            if(Ak!=0):
                                productQ2 = productQ2 + Ak*busdata[k][4]*np.sin(busdata[i][5]-busdata[k][5]-anglek)
                        J4[j4count[0],j4count[1]] = -busdata[i][4]*np.sqrt(YBus[i][i].real**2 + YBus[i][i].imag**2)*np.sin(np.arctan(YBus[i][i].imag/YBus[i][i].real)) + productQ2
                        j4count[1] += 1
            j4count[0] += 1
    return(np.vstack((np.hstack((J1,J2)),np.hstack((J3,J4)))))


# updateStateVariables is a function that finds the updated values of the state variables defined.
def updateStateVariables(inverse,myangles,myVs,pMismatch,qMismatch,busdata,nGen,nLoad):
    powersMismatch = np.vstack(np.hstack((pMismatch, qMismatch)))
    for i in range(1,nGen+nLoad+1):
        for j in range(0,nGen+nLoad+1):
            myangles[i] = myangles[i] + inverse[i-1][j]*powersMismatch[j]
        busdata[i][5] = myangles[i]
    for i in range(1,nGen+nLoad+1):
        currentLoad = nGen+nLoad-1
        if(busdata[i][1]==2):
            currentLoad += 1
            for j in range(0,nGen+nLoad+1):
                myVs[i] = myVs[i] + inverse[currentLoad][j]*powersMismatch[j]
            busdata[i][4] = myVs[i]
    return myangles, myVs


if __name__ == '__main__':
    import random
    from newton_raphson import check_3bus       #imports the function "check_3bus"

    np.set_printoptions(precision=3)            #numbers include 3 decimals
    np.set_printoptions(suppress=True)

    print('Testing 3-bus NR loadflow solver for candidate: XXXXX')

    for i in range(0, 5):                   #conduct 5 tests
        # All values in pu
        R12 = random.randint(3, 6) / 100  # R between 0.03 and 0.06
        X12 = random.randint(20, 30) / 100  # X between 0.2 and 0.3

        R13 = random.randint(2, 4) / 100  # R between 0.02 and 0.04
        X13 = random.randint(10, 20) / 100  # X between 0.1 and 0.2

        P2 = random.randint(1, 10) / 10  # P between 0.1 and 1
        V2 = 1 + random.randint(-5, 5) / 100  # V between 0.95 and 1.05

        P3 = -random.randint(1, 10) / 10  # P between -0.1 and -1
        Q3 = -random.randint(1, 10) / 10  # Q between -0.1 and -1

        myangles, myVs = my_threebus_NR_1xxxx(R12, X12, R13, X13, P2, V2, P3, Q3)

        print('Test {}: '.format(i + 1), end='')
        check_3bus(R12, X12, R13, X13, P2, V2, P3, Q3,myangles, myVs,set_print= False)