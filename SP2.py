#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 22:13:09 2021

@author: apple
"""
import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

import MP2 

#create model
sub = pyo.ConcreteModel(name="SP2")

#constants
#denom = [206, 274, 220] #nominal demand
#deunc = [40, 40, 40]  #demand uncertainty coefficients
C = np.array([[22, 33, 24],
              [33, 23, 30],
              [20, 25, 27]]) #unit costs

#variables
sub.x = pyo.Var(MP2.master.S, MP2.master.D, within = pyo.NonNegativeReals) #transportation variables
sub.g = pyo.Var(MP2.master.D, bounds = (0,1)) #variables in uncertainty set
sub.lambd = pyo.Var(MP2.master.D, within=pyo.NonNegativeReals) #dual variable of the second stage inner problem
sub.pi = pyo.Var(MP2.master.S, within=pyo.NonNegativeReals) #dual variable
sub.beta = pyo.Var(MP2.master.S, within=pyo.Binary)
sub.v = pyo.Var(MP2.master.D, within=pyo.Binary)
sub.alpha = pyo.Var(MP2.master.S, MP2.master.D, within=pyo.Binary)

#uncertainty set
demand = [206 + 40 * sub.g[0], 274 + 40 * sub.g[1], 220 + 40 * sub.g[2]] #uncertain demand at site i

#big M
#Ma = np.zeros((3, 3))
#Mpi = np.zeros((3))
#Mlambd = np.zeros((3))
#for i in range(3):
#    for j in range(3):
#        Ma[i][j] = min(demand[j], pyo.value(MP2.master.z[i]))
#         Mpi[i] = max(C[i][0], C[i][1], C[i][2])
#        Mlambd[i] = max(C[0][i], C[1][i], C[2][i])
M = 10^9


#objective
#def sub_obj_rule(model):
#    return sum(C[s, d]*sub.x[s,d] for s in MP2.master.S for d in MP2.master.D)
#sub.obj = pyo.Objective(rule = sub_obj_rule, sense=pyo.maximize)
sub.obj = pyo.Objective(expr = sum(C[s,d]*sub.x[s,d] for s in MP2.master.S for d in MP2.master.D),
                        sense = pyo.maximize)

# original constraints
def con1(model, s, d):
    return sum(sub.x[s,d] for d in MP2.master.D) <= pyo.value(MP2.master.z[s]) 
sub.capacityLimit = pyo.Constraint(MP2.master.S, MP2.master.D, rule=con1)

def con2(model, s, d):
    return sum(sub.x[s,d] for s in MP2.master.S) >= demand[d] 
sub.meetDemand = pyo.Constraint(MP2.master.S, MP2.master.D, rule=con2)

#uncertainty set constraints
def con4(model, d):
    return sum(sub.g[d] for d in MP2.master.D) <= 1.8
sub.G1 = pyo.Constraint(MP2.master.D, rule=con4)

def con5(model):
    return sub.g[0]+sub.g[1] <= 1.2
sub.G2 = pyo.Constraint(rule=con5)

#dual variable constraints
def con6(model, s, d):
    return sub.lambd[d] - sub.pi[s] <= C[s,d]
sub.dualCons = pyo.Constraint(MP2.master.S, MP2.master.D, rule=con6)

#complementary slackness
def con7(model, s):
    return sub.pi[s] <= M*sub.beta[s]
sub.compSlack1 = pyo.Constraint(MP2.master.S, rule = con7)

def con8(model, s, d):
    return pyo.value(MP2.master.z[s]) - sum(sub.x[s,d] for d in MP2.master.D) <= M*(1-sub.beta[s])
sub.compSlack2 = pyo.Constraint(MP2.master.S, MP2.master.D, rule=con8)

def con9(model, d):
    return sub.lambd[d] <= M*sub.v[d]
sub.compSlack3 = pyo.Constraint(MP2.master.D, rule=con9)

def con10(model, s, d):
    return sum(sub.x[s,d] for s in MP2.master.S) - demand[d] <= M*(1-sub.v[d])
sub.compSlack4 = pyo.Constraint(MP2.master.S, MP2.master.D, rule=con10)

def con11(model,s,d):
    return sub.x[s,d] <= M*sub.alpha[s,d]
sub.compSlack5 = pyo.Constraint(MP2.master.S, MP2.master.D, rule=con11)

def con12(model,s,d):
    return C[s,d] - sub.lambd[d] + sub.pi[s] <= M*(1-sub.alpha[s,d])
sub.compSlack6 = pyo.Constraint(MP2.master.S, MP2.master.D, rule=con12)

#solve
sub_opt = pyo.SolverFactory('glpk') #glpk
sub_opt.solve(sub)



    



