from numpy import *
import numpy as np
from cmath import *
import cmath as cm
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


if __name__ == '__main__':
<<<<<<<
<<<<<<< HEAD
    busmatrix = getbusdata("C:/Users/jonas/Downloads/project_busdata1.xlsx")
    linematrix = getlinedata("C:/Users/jonas/Downloads/project_busdata2.xlsx")
    buses = getbuses("C:/Users/jonas/Downloads/project_busdata1.xlsx")
    #print(busmatrix[0])
    print(linematrix)
   
=======
    busmatrix = getbusdata("/Users/emiljakobsen/Documents/NTNU 4.klasse/7. semester/Power System Analysis/Project/Project_PowerSystem/project_busdata.xlsx")
=======
    busmatrix = getbusdata(https://github.com/emilaune/Project_PowerSystem/blob/main/project_busdata.xlsx)
>>>>>>>
    buses = getbuses("/Users/emiljakobsen/Documents/NTNU 4.klasse/7. semester/Power System Analysis/Project/project_busdata.xlsx")
    print(busmatrix)
    print("heipåeskil")
>>>>>>> readme-edits
print("hei")