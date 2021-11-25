from task7 import *
from task8 import *
if __name__ == '__main__':
    print("\nTask 7 \nQuestion a)")
    posZbus = positivesequencebus()
    print("\n", posZbus)
    print("\nquestion b) \n")
    Ifa, Ifb, Ifc = symmetricalfaultcurrent(posZbus)
    print("\nquestion c) \n")
    postPhaseVoltages, postVoltages = busvoltages(posZbus)
    print("\nquestion d) \n")
    currentfrom4to5(postVoltages)

    print("\n\nTask 8:\nQuestion a)\n")
    zeroZbus = zerosequencebus()
    print("\nZero sequence Zbus:\n", zeroZbus)
    negZbus = negativesequencebus()
    print("\nNegative sequence Zbus:\n", negZbus)
    print("\nQuestion b)\n")
    unsCurrent = unsymmetricalfaultcurrent(zeroZbus, posZbus, negZbus)
    print("\nQuestion c)\n")
    unsymmetricalPostPhaseVs, unsymmetricalPostVs = unsymmetricalbusvoltages(zeroZbus, posZbus, negZbus, unsCurrent)
    print("\nQuestion d)\n")
    current45 = unsymmetricalcurrentfrom4to5(unsymmetricalPostVs)



