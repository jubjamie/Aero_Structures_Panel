from scipy.optimize import minimize
import math
import matplotlib.pyplot as plt
import numpy as np

mref = 93.566  # kNm
ibchord = 3.333  # m
Esk = 72  # GPa
Est = 88  # GPa
kk = 2.3  # Skin Buckling Coefficient
kt = 0.18  # Stiffener Buckling Coefficient
maxstrain = 0.0036  # m/m
density = 1600  # kg/m3 (of CFRP)
ribSpace = 1200  # mm
bsk = 200  # mm
flangewidth = 35

Nx = 1.5*mref/(0.1*0.4*ibchord*ibchord*1000)  # Safe Nx with 1.5 Safety Factor

# Storage for iteration results
bstRecord = []
tstRecord = []
tskRecord = []
massRecord = []


def mass(dim):
    """
    Objective Function
    :param dim: Length 3 dimension array in order bst, tst, tsk
    :return: The unit mass of the panel
    """
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]

    panel_area = (bsk*tsk) + (bst*tst) + (flangewidth*tst)  # Panel Surface Area
    unitvol = panel_area*1000  # mm3
    unitvolm = unitvol/1000000000  # m3
    panelmass = unitvolm * density  # kg
    unitmass = panelmass / (bsk/1000)  # kg
    return unitmass


def skinBuckle(dim):
    """
    Calculate the reserve factor for Skin Buckling by calculating the critical load and comparing to the applied load.
    Return the reserve factor adjusted for inequalities based around ==0
    :param dim: Length 3 dimension array in order bst, tst, tsk
    :return: Inequality adjusted reserve factor
    """
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]

    epsilonk = kk*((tsk/bsk))**2
    Et = (Esk*tsk)+(Est*((bst*tst)/bsk))
    Nsk = Et*epsilonk  # Critical Load
    rsf = Nsk/Nx
    return rsf - 1  # Using a target Reserve Factor of 1


def stiffenerBuckle(dim):
    """
    Calculate the reserve factor for Stiffener Buckling by calculating the critical load and
    comparing to the applied load.
    Return the reserve factor adjusted for inequalities based around ==0
    :param dim: Length 3 dimension array in order bst, tst, tsk
    :return: Inequality adjusted reserve factor
    """
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]

    epsilont = kt * ((tst / bst)) ** 2
    Et = (Esk * tsk) + (Est * ((bst * tst) / bsk))
    Nst = Et*epsilont  # Critical Load
    rsf = Nst/Nx
    return rsf - 1  # Using a target Reserve Factor of 1


def matFail(dim):
    """
    Calculate the reserve factor for Material Failure by calculating the critical load and
    comparing to the applied load.
    Return the reserve factor adjusted for inequalities based around =>0
    :param dim: Length 3 dimension array in order bst, tst, tsk
    :return: Inequality adjusted reserve factor
    """
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]

    Et = (Esk * tsk) + (Est * ((bst * tst) / bsk))
    Nmat = Et*maxstrain  # Critical Load
    rsf = Nmat/Nx
    return rsf - 1.1  # Using a target Reserve Factor of >=1.1


def eulerBuckle(dim):
    """
    Calculate the reserve factor for Euler Buckling Failure by calculating the critical load and
    comparing to the applied load.
    Return the reserve factor adjusted for inequalities based around =>0
    :param dim: Length 3 dimension array in order bst, tst, tsk
    :return: Inequality adjusted reserve factor
    """
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]

    ZEAZ = (Est*bst*tst*((tsk/2)+(bst/2)))
    ZEA = (Est*bst*tst)+(Esk*tsk*bsk)
    zbar = ZEAZ/ZEA  # Neutral Axis

    EIbar = ((Esk*bsk*(tsk**3))/12)+(Esk*bsk*tsk*(zbar**2))+((Est*tst*bst**3)/12)+\
            (Est*bst*tst*(((bst/2)+(tsk/2)-zbar)**2))  # Using Parallel Axis Theorm
    NxEuler = ((math.pi**2)*EIbar)/(ribSpace**2*bsk)  # Critical Load
    rsf = NxEuler/Nx
    return rsf - 1.1  # Using a target Reserve Factor of >=1.1


def tester(dim):
    """
    This function allows for testing of a set of dimensions, printing key info as a result.
    :param dim: Length 3 dimension array in order bst, tst, tsk
    :return: Printed output to stdout
    """
    print("Unit Mass: " + str(mass(dim)))
    print("Skin Buckle RSF: " + str(skinBuckle(dim)+1))
    print("Stiffener Buckle RSF: " + str(stiffenerBuckle(dim)+1))
    print("Mat Fail RSF: " + str(matFail(dim)+1.1))
    print("Euler Fail RSF: " + str(eulerBuckle(dim)+1.1))


def callbackMonitor(dim):
    """
    The callback function for the optimize process. Manages recording of iteration values and monitoring current
    optimiser solutions.
    :param dim: Length 3 dimension array in order bst, tst, tsk
    :return: Printed Iteration solutions to stdout.
    """
    bstRecord.append(dim[0])
    tstRecord.append(dim[1])
    tskRecord.append(dim[2])
    massRecord.append(mass(dim))
    print(dim, mass(dim))


# Define a dictionary of constraints for the optimiser.
conds = ({'type': 'eq', 'fun': skinBuckle},  # Define Equality for Skin Buckling
         {'type': 'eq', 'fun': stiffenerBuckle},  # Define Equality for Skin Buckling
         {'type': 'ineq', 'fun': matFail},  # Define Inequality for Material Failure
         {'type': 'ineq', 'fun': eulerBuckle})  # Define Inequality for Euler Buckling Failure

bnds = ((0, None), (0, None), (0, None))  # Define all solutions to be non-negative

# Main Optimiser Function
res = minimize(mass, [90, 15, 15], method='SLSQP', bounds=bnds, constraints=conds, callback=callbackMonitor)

print(tester(res['x']))

bstRecordAdjust = [float(k)/10 for k in bstRecord]  # Convert bst to cm to fit on graph nicely.

# Define tick for graphing
major_ticks1 = np.arange(3, 10, 1)
minor_ticks1 = np.arange(3, 10, 0.25)
major_ticks2 = np.arange(900, 2000, 100)
minor_ticks2 = np.arange(500, 2000, 50)

# Plot graphs
fig, ax1 = plt.subplots()
ax1.plot(bstRecordAdjust, 'b-', label='bst cm')
ax1.plot(tstRecord, 'g--', label='tst mm')
ax1.plot(tskRecord, 'c-.', label='tsk mm')
ax1.set_xlabel('Iteration')
ax1.set_ylabel('Value Attempt', color='b')
ax1.set_yticks(major_ticks1)
ax1.set_yticks(minor_ticks1, minor=True)
ax1.tick_params('y', colors='b')
ax1.legend(loc='upper right', shadow=True)

# Plot 2nd y-axis graph
ax2 = ax1.twinx()
ax2.plot(massRecord, 'r-')
ax2.set_ylabel('Unit Mass', color='r')
ax2.tick_params('y', colors='r')

# Create matplotlib windows with GUI
plt.show()
