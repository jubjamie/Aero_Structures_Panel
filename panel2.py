from scipy.optimize import minimize
import math
import matplotlib.pyplot as plt

mref = 93.566
ibchord = 3.333
Esk = 72
Est = 88
kk = 2.3
kt = 0.18
maxstrain = 0.0036
ribSpace = 1200
bsk = 200

Nx = 1.5*mref/(0.1*0.4*ibchord*ibchord*1000)

dimRecord = []

def area(dim):
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]

    panel_area = (bsk*tsk) + (bst * tst)
    return panel_area

def skinBuckle(dim):
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]
    epsilonk = kk*((tsk/bsk))**2
    Et = (Esk*tsk)+(Est*((bst*tst)/bsk))
    Nsk = Et*epsilonk
    rsf = Nsk/Nx
    return rsf - 1


def stiffenerBuckle(dim):
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]
    epsilont = kt * ((tst / bst)) ** 2
    Et = (Esk * tsk) + (Est * ((bst * tst) / bsk))
    Nst = Et*epsilont
    rsf = Nst/Nx
    return rsf - 1


def matFail(dim):
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]
    Et = (Esk * tsk) + (Est * ((bst * tst) / bsk))
    Nmat = Et*maxstrain
    rsf = Nmat/Nx
    return rsf - 1.1


def eulerBuckle(dim):
    bst = dim[0]
    tst = dim[1]
    tsk = dim[2]
    ZEAZ = (Est*bst*tst*((tsk/2)+(bst/2)))
    ZEA = (Est*bst*tst)+(Esk*tsk*bsk)
    zbar = ZEAZ/ZEA
    EIbar = ((Esk*bsk*((tsk)**3))/12)+(Esk*bsk*tsk*((zbar)**2))+((Est*tst*(bst)**3)/12)+(Est*bst*tst*(((bst/2)+(tsk/2)-zbar)**2))
    NxEuler = ((math.pi**2)*EIbar)/(ribSpace**2*bsk)
    rsf = NxEuler/Nx
    return rsf - 1.1


def tester(dim):
    print("Area: " + str(area(dim)))
    print("Skin Buckle RSF: " + str(skinBuckle(dim)+1))
    print("Stiffener Buckle RSF: " + str(stiffenerBuckle(dim)+1))
    print("Mat Fail RSF: " + str(matFail(dim)+1.1))
    print("Euler Fail RSF: " + str(eulerBuckle(dim)+1.1))


def callbackMonitor(dim):
    dimRecord.append(dim)
    print(dim, area(dim))

# print(tester([47.5, 3.3468, 3.94]))

conds = ({'type': 'eq', 'fun': skinBuckle},
         {'type': 'eq', 'fun': stiffenerBuckle},
         {'type': 'ineq', 'fun': matFail},
         {'type': 'ineq', 'fun': eulerBuckle})

bnds = ((0, None), (0, None), (0, None))

res = minimize(area, [90, 10, 10], method='SLSQP', bounds=bnds, constraints=conds, callback=callbackMonitor)
print(tester(res['x']))

plt.plot([bst[0] for bst in dimRecord])
plt.ylabel('bst Attempts')
plt.show()