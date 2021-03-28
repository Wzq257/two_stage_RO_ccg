#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 14:50:31 2021

@author: apple
"""
import pyomo.environ as pyo
from pyomo.opt import SolverFactory 
import MP2
import SP2
import numpy as np

eps = 10**(-4)
S = range(3)
D = range(3)
k = 1
O = [1]

#for o in O:
#    xx = 'add_x' + str(O[k-1])
    #new_x[name] = MP2.master.add_component(name, pyo.Var(S,D,within=pyo.NonNegativeReals))
#    MP2.master.xx = pyo.Var(S,D,within=pyo.NonNegativeReals)
#    new_x[xx] = MP2.master.name
MP2.master.xx1 = pyo.Var(S,D,within=pyo.NonNegativeReals)
MP2.master.xx2 = pyo.Var(S,D,within=pyo.NonNegativeReals)
MP2.master.xx3 = pyo.Var(S,D,within=pyo.NonNegativeReals)
MP2.master.xx4 = pyo.Var(S,D,within=pyo.NonNegativeReals)
MP2.master.xx5 = pyo.Var(S,D,within=pyo.NonNegativeReals)
var_dic = {1:MP2.master.xx1, 2:MP2.master.xx2, 3:MP2.master.xx3, 4:MP2.master.xx4, 5:MP2.master.xx5}

    
while MP2.UpperB - MP2.LowerB > eps:
    print("start")
    #add new constraints
    def add_bound(model):
        return MP2.master.eta >= sum(SP2.C[s][d]* var_dic[k][s,d] for s in S for d in D) 
    MP2.master.addBound = pyo.Constraint(rule=add_bound)
    print("added bound")
        
    def cap_cut(model,s):
        return sum(var_dic[k][s,d] for d in D) <= MP2.master.z[s]
    MP2.master.capCut = pyo.Constraint(S, rule=cap_cut)
    print("added capacity cut")

    def demand_cut(model, d):
        return sum(var_dic[k][s,d] for s in S) >= pyo.value(SP2.demand[d])
    MP2.master.demandCut = pyo.Constraint(D, rule=demand_cut)
    print("added demand cut")
    
    #go to step 2 (solve MP2 and update lower bound)
    master_opt = pyo.SolverFactory('glpk')
    master_opt.solve(MP2.master) 
    MP2.LowerB = MP2.master.obj()
    print("updated lower bound")
    print(MP2.LowerB)
    
    #delete some old constraints in SP2
    SP2.sub.del_component(SP2.sub.capacityLimit)
    SP2.sub.del_component(SP2.sub.meetDemand)
    SP2.sub.del_component(SP2.sub.dualCons)
    SP2.sub.del_component(SP2.sub.compSlack1)
    SP2.sub.del_component(SP2.sub.compSlack2)
    SP2.sub.del_component(SP2.sub.compSlack3)
    SP2.sub.del_component(SP2.sub.compSlack4)
    SP2.sub.del_component(SP2.sub.compSlack5)
    SP2.sub.del_component(SP2.sub.compSlack6)
    
    SP2.sub.del_component(SP2.sub.capacityLimit_index)
    SP2.sub.del_component(SP2.sub.meetDemand_index)
    SP2.sub.del_component(SP2.sub.dualCons_index)
    SP2.sub.del_component(SP2.sub.dualCons_index_0)
    SP2.sub.del_component(SP2.sub.dualCons_index_1)
    SP2.sub.del_component(SP2.sub.compSlack1_index)
    SP2.sub.del_component(SP2.sub.compSlack2_index)
    SP2.sub.del_component(SP2.sub.compSlack3_index)
    SP2.sub.del_component(SP2.sub.compSlack4_index)
    SP2.sub.del_component(SP2.sub.compSlack5_index)
    SP2.sub.del_component(SP2.sub.compSlack5_index_0)
    SP2.sub.del_component(SP2.sub.compSlack5_index_1)
    SP2.sub.del_component(SP2.sub.compSlack6_index)
    SP2.sub.del_component(SP2.sub.compSlack6_index_0)
    SP2.sub.del_component(SP2.sub.compSlack6_index_1)
    print("finished deleting")
    
    #add back constraints with updated values
    for s in S:
        for d in D:
            SP2.Ma[s,d] = SP2.max_D[d]
            SP2.Mpi[s] = max(SP2.C[s][0], SP2.C[s][1], SP2.C[s][2])
            SP2.Mlambd[d] = max(SP2.C[0][d], SP2.C[1][d], SP2.C[2][d])
    print("big Ms set")
    
    def con1(model,s):
        return sum(SP2.sub.x[s,d] for d in D) <= pyo.value(MP2.master.z[s]) 
    SP2.sub.capacityLimit = pyo.Constraint(S, rule=con1)

    def con2(model, s, d):
        return sum(SP2.sub.x[s,d] for s in S) >= SP2.demand[d] 
    SP2.sub.meetDemand = pyo.Constraint(S, D, rule=con2)
    
    #dual variable constraints
    def con5(model, s, d):
        return SP2.sub.lambd[d] - SP2.sub.pi[s] <= SP2.C[s][d]
    SP2.sub.dualCons = pyo.Constraint(S, D, rule=con5)

    #complementary slackness
    def con6(model, s):
        return SP2.sub.pi[s] <= SP2.Mpi[s]*SP2.sub.beta[s]
    SP2.sub.compSlack1 = pyo.Constraint(S, rule = con6)

    def con7(model, s, d):
        return pyo.value(MP2.master.z[s]) - sum(SP2.sub.x[s,d] for d in D) <= SP2.Mpi[s]*(1-SP2.sub.beta[s])
    SP2.sub.compSlack2 = pyo.Constraint(S, D, rule=con7)

    def con8(model, d):
        return SP2.sub.lambd[d] <= SP2.Mlambd[d]*SP2.sub.v[d]
    SP2.sub.compSlack3 = pyo.Constraint(D, rule=con8)

    def con9(model, s, d):
        return sum(SP2.sub.x[s,d] for s in S) - SP2.demand[d] <= SP2.Mlambd[d]*(1-SP2.sub.v[d])
    SP2.sub.compSlack4 = pyo.Constraint(S, D, rule=con9)

    def con10(model,s,d):
        return SP2.sub.x[s,d] <= SP2.Ma[s,d]*SP2.sub.alpha[s,d]
    SP2.sub.compSlack5 = pyo.Constraint(S, D, rule=con10)

    def con11(model,s,d):
        return SP2.C[s][d] - SP2.sub.lambd[d] + SP2.sub.pi[s] <= SP2.Ma[s,d]*(1-SP2.sub.alpha[s,d])
    SP2.sub.compSlack6 = pyo.Constraint(S, D, rule=con11)
    print("constraints redefined")
    
    #solve the subprobelm; update upper bound
    sub_opt = pyo.SolverFactory('glpk') #glpk
    sub_opt.solve(SP2.sub)
    Q = SP2.sub.obj()
    MP2.UpperB = min(MP2.UpperB, MP2.master.obj() - pyo.value(MP2.master.eta) + Q)
    print("updated upper bound")
    print(MP2.UpperB)
    
    k = k+1
    print(k)
    
print(MP2.LowerB)
    
    
        
        
    
