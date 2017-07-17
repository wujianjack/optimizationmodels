from __future__ import division, print_function

from gurobipy import *


rollwidth = 115
size = [25, 40, 50, 55, 70]
amount = [50, 36, 24, 8, 30]
nwidth = 5

MAX_CGTIME = 1000

def reportRMP(model):
    if model.status == GRB.OPTIMAL:
        print("Using ", model.objVal, " rolls\n")
        
        var = model.getVars()
        for i in range(model.numVars):
            print(var[i].varName, " = ", var[i].x)
        print("\n")
        
        con = model.getConstrs()
        for i in range(model.numConstrs):
            print(con[i].constrName, " = ", con[i].pi)
        print("\n")
        

def reportSUB(model):
    if model.status == GRB.OPTIMAL:
        print("Pi: ", model.objVal, "\n")
        
        if model.objVal <= 1e-6:
            var = model.getVars()
            
            for i in range(model.numVars):
                print(var[i].varName, " = ", var[i].x)
            print("\n")

def reportMIP(model):
    if model.status == GRB.OPTIMAL:
        print("Best MIP Solution: ", model.objVal, " rolls\n")
        
        var = model.getVars()
        for i in range(model.numVars):
            print(var[i].varName, " = ", var[i].x)

try:
    rmp = Model("rmp")
    sub = Model("sub")
    
    rmp.setParam("OutputFlag", 0)
    sub.setParam("OutputFlag", 0)
    
    # construct RMP
    rmp_var = []
    for i in range(nwidth):
        rmp_var.append(rmp.addVar(0.0, GRB.INFINITY, 1.0, GRB.CONTINUOUS, "rmp_" + str(i)))
    
    rmp_con = []
    row_coeff = [0.0] * nwidth
    for i in range(nwidth):
        row_coeff[i] = int(rollwidth / size[i])
        rmp_con.append(rmp.addConstr(quicksum(rmp_var[j] * row_coeff[j] for j in range(nwidth)) >= amount[i], "rmpcon_" + str(i)))
        row_coeff[i] = 0.0
    
    rmp.setAttr("ModelSense", GRB.MINIMIZE)
    # end RMP
    
    # construct SUB
    sub_var = []
    for i in range(nwidth):
        sub_var.append(sub.addVar(0.0, GRB.INFINITY, 0.0, GRB.INTEGER, "sub_" + str(i)))
    
    sub.addConstr(quicksum(sub_var[i] * size[i] for i in range(nwidth)) <= rollwidth, "subcon")
    # end SUB
    
    print("               *** Column Generation Loop ***               \n")
    for i in range(MAX_CGTIME):
        print("Iteration: ", i, "\n")
        
        rmp.optimize()
        reportRMP(rmp)
        
        rmp_pi = rmp.getAttr("Pi", rmp.getConstrs())
        
        sub.setObjective(1 - quicksum(sub_var[i] * rmp_pi[i] for i in range(nwidth)), GRB.MINIMIZE)
        sub.optimize()
        reportSUB(sub)
        
        if sub.objVal > -1e-6:
            break
        
        rmp_coeff = sub.getAttr("X", sub.getVars())
        rmp_col = Column(rmp_coeff, rmp_con)
        rmp.addVar(0.0, GRB.INFINITY, 1.0, GRB.CONTINUOUS, "cg_" + str(i), rmp_col)
    
    print("               *** End Loop ***               \n")
    
    mip_var = rmp.getVars()
    for i in range(rmp.numVars):
        mip_var[i].setAttr("VType", GRB.INTEGER)
    
    rmp.optimize()
    reportMIP(rmp)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError:
    print('Encountered an attribute error')