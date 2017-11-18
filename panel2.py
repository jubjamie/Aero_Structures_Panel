from scipy.optimize import minimize
import math

mref = 93.566
ibchord = 3.333
Esk = 72
Est = 88
kk = 2.3
kt = 0.18
maxstrain = 0.0036
ribSpace = 1200
bsk = 200

def area(bst, tst, tsk):
    """
    Calculate the area as the objective function
    :param bst: Stiffener Height
    :param tst: Stiffener Width
    :param tsk: Skin Height
    :return: area of panel
    """
    panel_area = (bsk*tsk) + (bst * tst)
    return panel_area

print(area(47.5, 3.3468, 3.94))