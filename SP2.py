#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Wzq
"""
import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

import MP2 

#create model
sub = pyo.ConcreteModel(name="SP2")

#constants
C = [[22, 33, 24],
     [33, 23, 30],
     [20, 25, 27]] #unit costs
max_D = [206 + 40, 274 + 40, 220 + 40]

#sets
S = range(3)
D = range(3)

#variables
sub.x = pyo.Var(S, D, within = pyo.NonNegativeReals) #transportation variables
sub.g = pyo.Var(D, bounds = (0,1)) #variables in uncertainty set
sub.lambd = pyo.Var(D, within=pyo.NonNegativeReals) #dual variable of the second stage inner problem
sub.pi = pyo.Var(S, within=pyo.NonNegativeReals) #dual variable
sub.beta = pyo.Var(S, within=pyo.Binary)
sub.v = pyo.Var(D, within=pyo.Binary)
sub.alpha = pyo.Var(S, D, within=pyo.Binary)

#uncertainty set
demand = [206 + 40 * sub.g[0], 274 + 40 * sub.g[1], 220 + 40 * sub.g[2]] #uncertain demand at site i

#big M
Ma = np.zeros((3, 3))
Mpi = np.zeros((3))
Mlambd = np.zeros((3))
for s in S:
    for d in D:
        Ma[s,d] = max_D[d]
        Mpi[s] = max(C[s][0], C[s][1], C[s][2])
        Mlambd[d] = max(C[0][d], C[1][d], C[2][d])

#objective
def sub_obj_rule(model):
    return sum(C[s][d]*sub.x[s,d] for s in S for d in D)
sub.obj = pyo.Objective(rule = sub_obj_rule, sense=pyo.maximize)
#sub.obj = pyo.Objective(expr = sum(C[s,d]*sub.x[s,d] for s in S for d in D),
#                        sense = pyo.maximize)

# original constraints
def con1(model,s):
    return sum(sub.x[s,d] for d in D) <= pyo.value(MP2.master.z[s]) 
sub.capacityLimit = pyo.Constraint(S,rule=con1)

def con2(model, d):
    return sum(sub.x[s,d] for s in S) >= demand[d] 
sub.meetDemand = pyo.Constraint(D, rule=con2)

#uncertainty set constraints
def con3(model):
    return sub.g[0] + sub.g[1] + sub.g[2] <= 1.8
sub.G1 = pyo.Constraint(rule=con3)

def con4(model):
    return sub.g[0]+sub.g[1] <= 1.2
sub.G2 = pyo.Constraint(rule=con4)

#dual variable constraints
def con5(model, s, d):
    return sub.lambd[d] - sub.pi[s] <= C[s][d]
sub.dualCons = pyo.Constraint(S, D, rule=con5)

#complementary slackness
def con6(model, s):
    return sub.pi[s] <= Mpi[s]*sub.beta[s]
sub.compSlack1 = pyo.Constraint(S, rule = con6)

def con7(model, s):
    return MP2.master.z[s].value - sum(sub.x[s,d] for d in D) <= Mpi[s]*(1-sub.beta[s])
sub.compSlack2 = pyo.Constraint(S, rule=con7)

def con8(model, d):
    return sub.lambd[d] <= Mlambd[d]*sub.v[d]
sub.compSlack3 = pyo.Constraint(D, rule=con8)

def con9(model, d):
    return sum(sub.x[s,d] for s in S) - demand[d] <= Mlambd[d]*(1-sub.v[d])
sub.compSlack4 = pyo.Constraint(D, rule=con9)

def con10(model,s,d):
    return sub.x[s,d] <= Ma[s,d]*sub.alpha[s,d]
sub.compSlack5 = pyo.Constraint(S, D, rule=con10)

def con11(model,s,d):
    return C[s][d] - sub.lambd[d] + sub.pi[s] <= Ma[s,d]*(1-sub.alpha[s,d])
sub.compSlack6 = pyo.Constraint(S, D, rule=con11)

#solve SP2 and update upper bound
sub_opt = pyo.SolverFactory('glpk') #glpk
sub_opt.solve(sub)
demand = [206 + 40 * sub.g[0].value, 274 + 40 * sub.g[1].value, 220 + 40 * sub.g[2].value]
Q = sub.obj()
MP2.UpperB = min(MP2.UpperB, MP2.master.obj() - pyo.value(MP2.master.eta) + Q)
MP2.UpperB





    



