#Master Problem
import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

UB = np.inf
LB = -np.inf
k = 0
O = []

#construct a concrete model
master = pyo.ConcreteModel(name="MP2")

#sets
master.S = pyo.RangeSet(0,2)
master.D = pyo.RangeSet(0,2)
#master.D = pyo.RangeSet(0,2,1)

#data
#f = [400, 414, 326] #fixed cost at site i
#a = [18, 25, 20] #unit capacity cost for site i
#C = [[22, 33, 24],
#     [33, 23, 30],
#     [20, 25, 27], ] #unit transportation cost from site i to site j

#constant
K = 800

#variables
master.y = pyo.Var(master.S, within = pyo.Binary)  #facility location variable
master.z = pyo.Var(master.S, within = pyo.NonNegativeReals)  #capacity variable
master.eta = pyo.Var(within = pyo.NonNegativeReals, initialize = 0)
master.g = pyo.Var(master.D, bounds = (0,1)) #variables in uncertainty set

#demand = [206 + 40 * master.g[0], 274 + 40 * sub.g[1], 220 + 40 * sub.g[2]] #uncertain demand at site i



#objective
def master_obj_rule(model):
    return (400*master.y[0] + 414*master.y[1] + 326*master.y[2]+
               18*master.z[0] + 25*master.z[1] + 24*master.z[2] + master.eta)
master.obj = pyo.Objective(rule = master_obj_rule)

#constraints
def C1(model,s):
    return master.z[s] <= 800*master.y[s]
master.C1 = pyo.Constraint(master.S, rule=C1)

#def C2(model, d):
#    return sum(master.g[d] for d in master.D) <= 1.8
#master.C2 = pyo.Constraint(master.D, rule=C2)

#def C3(model):
#    return master.g[0]+master.g[1] <= 1.2
#master.C3 = pyo.Constraint(rule=C3)

def C4(model):
    return sum(master.z[s] for s in master.S) >= 206 + 274 + 220 + 1.8*40
master.C4 = pyo.Constraint(rule=C4)
#This is to guarantee feasibility

#solve
master_opt = pyo.SolverFactory('glpk')
master_opt.solve(master) 
LowerB = master.obj() + pyo.value(master.eta)
LowerB













