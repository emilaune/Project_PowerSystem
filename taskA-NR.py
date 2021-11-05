# Newton-Raphson Algorithm for solving a n-bus Load Flow Problem
# TET4115 Power System Analysis; semester: Fall 2021
# Author: Eivind Falk

import numpy as np
import cmath
import xlrd

class Bus:
    def __init__(self, data):  # data given as [number, type, P_specified, Q_specified, delta, V,
        self.number = data[0]
        self.type = data[1]
        self.P_specified = data[2]
        self.Q_specified = data[3]
        self.delta = data[4]
        self.V = data[5]
        self.P_lim = data[6]
        self.Q_lim = data[7]

class Line:
    def __init__(self, linedata):
        self.num = linedata[0]
        self.R = linedata[1]
        self.X = linedata[2]
        self.y = 1/(self.R+1j*self.X)


def getbusdata(filepath):
    loc = filepath
    op = xlrd.open_workbook(loc)
    a = op.sheet_by_index(0)
    busmatrix = [[], [], [], [], []]  # hard?
    for i in range(5):
        busdata = [0, '', 0, 0, 0, 0, 0, 0]  # had to define this in every loop, unless - hard?
        for l in range(8):  # the busmatrix values is updated as busdata gets updated and contains only values from bus5
            busdata[l] = a.cell_value(i+1, l)
        busmatrix[i] = busdata
    return busmatrix


def getbuses(filepath):  # takes the filepath as input parameter and uses getbusdata to obtain data for all the buses
    busmatrix = getbusdata(filepath)
    # initializing all the buses with data from busmatrix - to hard - softer?
    bus1 = Bus(busmatrix[0])
    bus2 = Bus(busmatrix[1])
    bus3 = Bus(busmatrix[2])
    bus4 = Bus(busmatrix[3])
    bus5 = Bus(busmatrix[4])
    buses = [bus1, bus2, bus3, bus4, bus5]
    return buses  # returns a list with all 5 buses

def getlinedata(filepath):
    loc = filepath
    op = xlrd.open_workbook(loc)
    a = op.sheet_by_index(0)
    linematrix = [[], [], [], [], []]
    for i in range (5):
        linedata = [0,0,0,0]
        for l in range(4):
            linedata[l] = a.cell_value(i+1,l)
        linematrix[i] = linedata
    return linematrix

def getlines(filepath):  
    linematrix = getlinedata(filepath)
    bus1 = Bus(linematrix[0])
    bus2 = Bus(linematrix[1])
    bus3 = Bus(linematrix[2])
    bus4 = Bus(linematrix[3])
    bus5 = Bus(linematrix[4])
    buses = [bus1, bus2, bus3, bus4, bus5]
    return buses 

    
def my_threebus_NR_1xxxx():
    

    linedata = linematrix        # [from bus to bus, resistance, reactance]
    busdata = busmatrix       # [bus no., bus type, P, Q, |V|, delta]
    YBus = findYBus(linedata)

    myangles = [0, 0, 0, 0, 0]     # all unknown bus voltage phase angles are initially set to 0 rad
    myVs = [1, 1, 1, 1, 1]         # all unknown bus voltage magnitudes are initially set to 1.0 pu

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
    busmatrix = getbusdata("C:/Users/jonas/Downloads/project_busdata1.xlsx")
    linematrix = getlinedata("C:/Users/jonas/Downloads/project_busdata2.xlsx")
    buses = getbuses("C:/Users/jonas/Downloads/project_busdata1.xlsx")
    print(busmatrix[0])
    print(linematrix[0])
