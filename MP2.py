"""
@author: Wzq
"""

#Master Problem
import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

UpperB = np.inf
LowerB = -np.inf
#construct a concrete model
master = pyo.ConcreteModel(name="MP2")

#sets
S = range(3)
D = range(3)

#constant
K = 800

#variables
master.y = pyo.Var(S, within = pyo.Binary)  #facility location variable
master.z = pyo.Var(S, within = pyo.NonNegativeReals)  #capacity variable
master.eta = pyo.Var(within = pyo.NonNegativeReals, initialize = 0)

#objective
def master_obj_rule(model):
    return (400*master.y[0] + 414*master.y[1] + 326*master.y[2]+
               18*master.z[0] + 25*master.z[1] + 24*master.z[2] + master.eta)
master.obj = pyo.Objective(rule = master_obj_rule)

#constraints
def C1(model,s):
    return master.z[s] <= 800*master.y[s]
master.C1 = pyo.Constraint(S, rule=C1)

def C2(model):
    return sum(master.z[s] for s in S) >= 206 + 274 + 220 + 1.8*40
master.C2 = pyo.Constraint(rule=C2) #This is to guarantee feasibility

#solve MP2 and update lower bound
master_opt = pyo.SolverFactory('glpk')
master_opt.solve(master) 
LowerB = master.obj()
LowerB














