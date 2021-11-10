# Newton-Raphson Algorithm for solving a n-bus Load Flow Problem
# TET4115 Power System Analysis; semester: Fall 2021
# Author: Eivind Falk


import numpy as np      # np.sin, np.cos
import cmath            # cmath.polar
import math             # math.nan

def partAquestion4():

    #           [i,     j,  r_ij,   x_ij,   y(OC)_ij]
    linedata = [[1,     2,  0.05,   0.25,   0.5],
                [2,     3,  0.05,   0.25,   0.05],
                [1,     4,  0.02,   0.2,    0.033],
                [2,     4,  0.02,   0.2,    0.033],
                [4,     5,  0.01,   0.1,    0.02]]

    #          [no, type,         P,          Q,   |V|, delta,     P_min,      P_max,      Q_min,     Q_max,     x_shunt]
    busdata = [[ 0,    1,       1.0,   math.nan,   1.0,     0,  math.nan,   math.nan,   math.nan,       0.5,    math.nan],
               [ 1,    2,      -0.6,       -0.3,   1.0,     0,  math.nan,   math.nan,   math.nan,  math.nan,    math.nan],
               [ 2,    0,  math.nan,   math.nan,   1.0,     0,  math.nan,       2.0,    math.nan,       1.5,    math.nan],
               [ 3,    2,      -0.6,       -0.2,   1.0,     0,  math.nan,   math.nan,   math.nan,  math.nan,        -1.0],
               [ 4,    2,      -0.5,       -0.4,   1.0,     0,  math.nan,   math.nan,   math.nan,  math.nan,    math.nan]]

    initVs = [1,1,1,1,1]
    typeswitch = [0,0,0,0,0]

    YBus = findYBus(linedata)
    #YBus = addTransformer(YBus0)           # transformers are not added
    myangles = [0, 0, 0, 0, 0]              # all unknown bus voltage phase angles are initially set to 0 rad
    myVs = [1, 1, 1, 1, 1]                  # all unknown bus voltage magnitudes are initially set to 1.0 pu

    iteration = 0
    print("iteration:",iteration)

    [pMismatch, qMismatch,b,P,Q] = isConverged(myangles, myVs, YBus, busdata)
    printIteration(busdata, pMismatch, qMismatch, P, Q)
    [busdata,typeswitch] = checkReactiveLimits(busdata,Q,typeswitch,initVs)
    [nGen,nLoad] = findSize(busdata)

    while(b == False):                  # while not converged
        iteration += 1
        print("Iteration: ",iteration)
        print(typeswitch)
        # Find new Jacobian
        J = findJacobian(YBus, busdata,nGen,nLoad)
        printJacobian(J)
        inverse = np.linalg.inv(J)

        # Update the state variables
        [myangles,myVs,busdata] = updateStateVariables(inverse,myangles,myVs,pMismatch,qMismatch,busdata,nGen,nLoad)

        # Check convergence
        [pMismatch, qMismatch, b, P, Q] = isConverged(myangles, myVs, YBus, busdata)

        # display the iteration results
        printIteration(busdata, pMismatch, qMismatch,P,Q)
        [busdata,typeswitch] = checkReactiveLimits(busdata,Q,typeswitch,initVs)
        [nGen,nLoad] = findSize(busdata)
    powerTransmissions(linedata, busdata, YBus)


def printJacobian(J):
    dim = len(J)
    gaplength = 5
    gap = ' '*gaplength
    header = 'JACOBIAN'
    space = int((dim)*(14/2)-len(header)/2)
    print("_"*int((space)),header,"_"*(int(space)))
    for i in range(0,dim):
        for j in range(0,dim):
            line = f"{J[i][j]:10f}{gap}"
            print(line,end = "")
        print(' ')


def checkReactiveLimits(busdata,Q,typeswitch,initVs):
    # Q_iL of both generators are 0 --> Q_i = Q_iG
    dim = len(busdata)
    for i in range(0,dim):

        if(busdata[i][1]==1): #maybe redundant4
            if(typeswitch[i]==0)and(Q[i]>busdata[i][9]):    #not typeswitched, violating limit

                typeswitch[i] = 1               # indicate that the bus is type switched
                busdata[i][1] = 2               # switch from P-|V| type to P-Q type
                busdata[i][3] = busdata[i][9]   # set Q_spec to Q_max
                busdata[i][4] = initVs[i]       # set current voltage equal to initial value

                print("-" * 150)
                print("Upper Reactive Power Limit on bus:",busdata[i][0]+1," is violated!!!")
                print("-" * 150)

            elif (typeswitch[i] == 0) and (Q[i] < busdata[i][8]):

                typeswitch[i] = 1                   # indicate that it is transformed
                busdata[i][1] = 2                   # switch from P-|V| to P-Q
                busdata[i][3] = busdata[i][8]       # set Q_spec to Q_min
                busdata[i][4] = initVs[i]           # set current voltage equal to initial value

                print("-" * 150)
                print("Lower Reactive Power Limit on bus:", busdata[i][0] + 1, " is violated!!!")
                print("-" * 150)

            elif(typeswitch[i]==1):#and(Q[i]<busdata[i][9]):  #typeswitched, not violating limit anymore
                print("should typeswitch back here")

                typeswitch[i] = 0               # indicate that it is back-transformed
                busdata[i][1] = 1               # switch from P-Q to P-|V| bus
        print(typeswitch[i], "bus:", i)

    return busdata,typeswitch


def powerTransmissions(linedata,busdata,YBus):

    Sbase = 100

    # header
    gap = ' '*2
    header = f"{'Line':15s}{gap}{'S_ij[MVA]':10s}{gap}{'S_ij(OC)[MVA]':10s}{gap}{'P_ij[MW]':10s}{gap}{'P_ij(OC)[MW]':10s}{gap}{'Q_ij[MVAr]':10s}{gap}{'Q_ij(OC)[MVAr]':10s}{gap}{'P_loss[MW]':10s}{gap}{'Q_loss[MVAr]'}"
    print("-"*150)
    print(header)
    print("-"*150)
    for i in range(0,len(linedata)):
        [A,theta] = cmath.polar(YBus[linedata[i][0]-1][linedata[i][1]-1])

        # transmitted power neglecting operating capacitance
        delta_i = busdata[linedata[i][0]-1][5]
        delta_j = busdata[linedata[i][1]-1][5]
        delta_ij = delta_i+(-delta_j)

        vSending = busdata[linedata[i][0]-1][4]
        vReceiving = busdata[linedata[i][1]-1][4]

        Y_ij = A*(complex(np.cos(theta),(-1)*np.sin(theta)))
        S_ij = (vSending**2-vSending*vReceiving*(complex(np.cos(delta_ij),np.sin(delta_ij))))*(-1)*Y_ij*Sbase

        # transmitted power / operating capacitance
        [yOCMagnitude,yOC_theta] = cmath.polar(complex(0,linedata[i][4]/2))
        S_ijOC = ((busdata[linedata[i][0]-1][4])**2)*complex(yOCMagnitude*np.cos(yOC_theta),-yOCMagnitude*np.sin(yOC_theta))*Sbase

        # transmission losses
        vPhasorSending = busdata[linedata[i][0] - 1][4] * (complex(np.cos(busdata[linedata[i][0] - 1][5]), np.sin(busdata[linedata[i][0] - 1][5])))
        vPhasorReceiving = busdata[linedata[i][1] - 1][4] * (complex(np.cos(busdata[linedata[i][1] - 1][5]), np.sin(busdata[linedata[i][1] - 1][5])))
        S_ij_loss = (np.abs(vPhasorSending-vPhasorReceiving))**2*(-1)*A*(complex(np.cos(theta),np.sin(theta)))*Sbase

        line = f"{linedata[i][0]}{' - '}{linedata[i][1]}{gap}{gap}{S_ij:15.3}{gap}{S_ij+S_ijOC:15.3}{gap}{np.real(S_ij):10.5}{gap}{gap}{np.real(S_ij+S_ijOC):10.5}{gap}{np.imag(S_ij):10.5}{gap}{np.imag(S_ij+S_ijOC):10.5}{gap}{np.abs(np.real(S_ij_loss)):10.5}{gap}{np.abs(np.imag(S_ij_loss)):10.5}"
        print(line)


def printIteration(busdata,pMismatch,qMismatch,P,Q):
    gap = ' '*4
    header = f"{'Bus no':6s}{gap}{'Bus Type':8s}{gap}{'P_spec':5s}{gap}{'P_min':5s}{gap}{'P_max':5s}{gap}{'Q_spec':5s}{gap}{'Q_min':5s}{gap}{'Q_max':8s}{gap}{'|V|':8s}{gap}{'theta':8s}{gap}{'P_calc':8s}{gap}{'Q_calc':4s}{gap}{'P_mismatch'}{gap}{'Q_mismatch'}"
    print("-"*153)
    print(header)
    print("-"*153)

    iterP = 0
    iterQ = 0
    for i in range(0,len(busdata)):

        if(busdata[i][1]==2):
            lines = f"{(busdata[i][0]+1):6d}{gap}{'Load':9}{gap}{busdata[i][2]:5.2f}{gap}{busdata[i][6]:5.2f}{gap}{busdata[i][7]:5.2f}{gap}{busdata[i][3]:5.2f}{gap}{busdata[i][8]:5.2f}{gap}{busdata[i][9]:5.2f}{gap}{busdata[i][4]:8.5}{gap}{float(busdata[i][5]):10.3}{gap}{P[i]:8.2}{gap}{Q[i]:8.2}{gap}{np.abs(pMismatch[iterP]):8.2}{gap}{np.abs(qMismatch[iterQ]):6.2}"

            iterQ += 1
            iterP += 1

        elif(busdata[i][1]==1):
            lines = f"{(busdata[i][0] + 1):6d}{gap}{'Generator':9}{gap}{busdata[i][2]:5.2f}{gap}{busdata[i][6]:5.2f}{gap}{busdata[i][7]:5.2f}{gap}{busdata[i][3]:5.2f}{gap}{busdata[i][8]:5.2f}{gap}{busdata[i][9]:5.2f}{gap}{busdata[i][4]:8.5}{gap}{float(busdata[i][5]):10.3}{gap}{P[i]:8.2}{gap}{Q[i]:8.2}{gap}{np.abs(pMismatch[iterP]):8.2}{gap}{math.nan:6}"

            iterP += 1

        else:
            lines = f"{(busdata[i][0] + 1):6d}{gap}{'Slack bus':9}{gap}{busdata[i][2]:5.2f}{gap}{busdata[i][6]:5.2f}{gap}{busdata[i][7]:5.2f}{gap}{busdata[i][3]:5.2f}{gap}{busdata[i][8]:5.2f}{gap}{busdata[i][9]:5.2f}{gap}{busdata[i][4]:8.5}{gap}{float(busdata[i][5]):10.3}{gap}{P[i]:8.2}{gap}{Q[i]:8.2}{gap}{math.nan:8}{gap}{math.nan:6}"

        print(lines)
    print('\n')

def findSize(busdata):
    nGen = 0  # number of generator buses in system
    nLoad = 0  # number of load buses in system
    n = len(busdata)
    for i in range(0, n):
        if (busdata[i][1] == 1):
            nGen += 1
        elif (busdata[i][1] == 2):
            nLoad += 1
    return nGen,nLoad



# isConverged calculates power injections at each bus by Load Flow Equations,
# finds the mismatches vector and returns the mismatches of active- and reactive powers
# of iteration "i" to be used in the iterations. Also, a boolean value is returned to indicate convergence
# or not convergence (i.e. perform another iteration)
def isConverged(myangles, myVs, YBus,busdata):

    n = len(busdata)
    V = myVs
    d = myangles
    P = np.zeros([n])                   # active power injections vector
    Q = np.zeros([n])                   # reactive power injections vector
    pMismatch = []                      # active power mismatch vector
    qMismatch = []                      # reactive power mismatch vector
    for i in range(0, n):
        for j in range(0,n):
            [A,angle] = cmath.polar(YBus[i][j])
            if (A != 0):
                P[i] = P[i] + V[i]*V[j]*A*np.cos(d[i] - d[j] - angle)
                Q[i] = Q[i] + V[i]*V[j]*A*np.sin(d[i] - d[j] - angle)
        if (busdata[i][1]==1)or(busdata[i][1]==2):
            pMismatch.append(busdata[i][2] - P[i])
            if (busdata[i][1]==2):
                qMismatch.append(busdata[i][3] - Q[i])

    b = True      #assume that the solution is converged
    tolerance = 0.1
    powersMismatch = np.vstack(np.hstack((pMismatch, qMismatch)))       #just stack all mismatches to run through all at once...
    for i in range(0,len(powersMismatch)):
        if(np.absolute(powersMismatch[i])>tolerance):                   #check tolerance
            b = False
            break
    return pMismatch,qMismatch,b,P,Q





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
                YBus[i - 1][i - 1] = YBus[i - 1][i - 1] + y[j] + linedata[j][4]/2       #adding half-admittance on line on last term

    dim = len(YBus)
    gaplength = 5
    gap = ' '*gaplength
    header = 'YBUS'
    space = int((dim)*(18/2)+int(len(header)/2))
    print("-"*space,header,"-"*space)
    for i in range(0,dim):
        for j in range(0,dim):
            [A,angle] = cmath.polar(YBus[i][j])
            line = f"{'|'}{A:6.3f}{'|'}{','}{angle:7.3f}{gap}"
            print(line,end = "")
        print('')
    print('\n')
    return YBus


def addTransformer(YBus):
    z_t1 = complex(0,0.2)
    z_t2 = complex(0,0.1)

    y_t1 = 1/z_t1
    y_t2 = 1/z_t2

    a_t1 = 1.00
    a_t2 = 0.98

    theta_t1 = (4/360)*2*3.14

    # new data on YBus

    # Y_11 + y_pq                       # Y_44 + (-y_pq)
    # Y_14 + (-y_pq/(as+jbs)            # Y_45 + (-y_pq/a)
    # Y_41 + (-y_pq/(as-jbs)            # Y_55 + y_pq/a**2
    # Y_44 + y_pq/(as**2+bs**2)         # Y_54 + (-y_pq/a)

    Y_11t1 = y_t1
    Y_14t1 = -y_t1/complex(a_t1*np.cos(theta_t1),a_t1*np.sin(theta_t1))
    Y_41t1 = -y_t1/complex(a_t1*np.cos(theta_t1),-a_t1*np.sin(theta_t1))
    Y_44t1 = y_t1/((a_t1*np.cos(theta_t1))**2+(a_t1*np.sin(theta_t1))**2)

    Y_44t2 = y_t2
    Y_45t2 = -y_t2/a_t2
    Y_54t2 = -y_t2/a_t2
    Y_55t2 = y_t2/a_t2**2

    YBus[0][0] = YBus[0][0] + Y_14t1
    YBus[0][3] = YBus[0][3] + Y_11t1
    YBus[3][0] = YBus[3][0] + Y_41t1
    YBus[3][3] = YBus[3][3] + Y_44t1

    YBus[3][3] = YBus[3][3] + Y_44t2
    YBus[3][4] = YBus[3][4] + Y_45t2
    YBus[4][3] = YBus[4][3] + Y_54t2
    YBus[4][4] = YBus[4][4] + Y_44t2

    return YBus





# findJacobian is the longest and most complex function used in the main program.
# Complexity comes from the calculation and inserting of entry values in the submatrices of the jacobian,
# when generator - and load buses are given at arbitrary positions in "linedata"
def findJacobian(YBus,busdata,nGen,nLoad):
    dim = len(busdata)
    J1 = np.zeros([nGen+nLoad,nGen+nLoad])
    J2 = np.zeros([nGen+nLoad,nLoad])
    J3 = np.zeros([nLoad,nGen+nLoad])
    J4 = np.zeros([nLoad,nLoad])

    j1count = [0,0]             # (#P x #delta)
    j2count = [0,0]             # (#P x #V)
    j3count = [0,0]             # (#Q x #delta)
    j4count = [0,0]             # (#Q x #V)

    # Calculate the (P-delta) - sensitivities
    for i in range(0,dim):
        if (busdata[i][1]==1)or(busdata[i][1]==2):      #find possible generator and load buses
            j1count[1] = 0
            for j in range(0,dim):
                if(busdata[j][1]==1)or(busdata[j][1]==2):   #find load or generator bus
                    [A,angle] = cmath.polar(YBus[i][j])
                    if(i!=j)and(A==0):
                        J1[j1count[0],j1count[1]] = 0
                        j1count[1] += 1
                    elif(i!=j)and(A!=0):
                        J1[j1count[0],j1count[1]] = busdata[i][4]*busdata[j][4]*A*np.sin(busdata[i][5]-busdata[j][5]-angle)
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
    for i in range(0,dim):
        if(busdata[i][1]==1)or(busdata[i][1]==2):
            [Ai,anglei] = cmath.polar(YBus[i][i])
            j2count[1] = 0
            for j in range(0,dim):
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
    for i in range(0,dim):
        if (busdata[i][1]==2):
            j3count[1] = 0
            for j in range(0,dim):
                if(busdata[j][1]==1)or(busdata[j][1]==2):
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
    for i in range(0,dim):
        if (busdata[i][1]==2):
            j4count[1] = 0
            for j in range(0,dim):
                if(busdata[j][1]==2):
                    [A,angle] = cmath.polar(YBus[i][j])
                    if(i!=j)and(A==0):
                        J4[j4count[0],j4count[1]] = 0
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
    dim = len(busdata)
    iter = 0

    for i in range(0,dim):
        if (busdata[i][1] == 1) or (busdata[i][1] == 2):
            for j in range(0,len(powersMismatch)):
                myangles[i] = myangles[i] + inverse[iter][j]*powersMismatch[j]
            busdata[i][5] = float(myangles[i])
            iter += 1

    for i in range(0,dim):
        if(busdata[i][1]==2):
            for j in range(0,len(powersMismatch)):
                myVs[i] = myVs[i] + inverse[iter][j]*powersMismatch[j]
            busdata[i][4] = float(myVs[i])
            iter += 1
    return myangles, myVs,busdata


partAquestion4()