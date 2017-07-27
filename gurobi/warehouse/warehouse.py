from __future__ import division, print_function

import gurobipy as GRBPY


# awkward restriction for 'callback'
def cbwarehouse(model, where):
    if where == GRBPY.GRB.Callback.MIPSOL:
        if model._iter >= 1:
            for i in range(model._nwarehouse):
                model._csupply[i].rhs = model.cbGetSolution(model._vmbuild[i]) * model._supply[i]
            
        model._sub.optimize()

        if model._sub.status == GRBPY.GRB.INFEASIBLE:
            print("Iteration: ", model._iter)
            print("Adding feasibility cut...\n")

            lazycut = GRBPY.quicksum(model._csupply[i].farkasdual * model._supply[i] * model._vmbuild[i] \
                                     for i in range(model._nwarehouse)) + \
                      sum(model._cdemand[i].farkasdual * model._demand[i] for i in range(model._nstore))
                
            model.cbLazy(lazycut >= 0)
            
            model._iter += 1
        elif model._sub.status == GRBPY.GRB.OPTIMAL:
            if model._sub.objval > model.cbGetSolution(model._maxshipcost) + 1e-6:
                print("Iteration: ", model._iter)
                print("Adding optimality cut...\n")

                lazycut = GRBPY.quicksum(model._csupply[i].pi * model._supply[i] * model._vmbuild[i] \
                                         for i in range(model._nwarehouse)) + \
                          sum(model._cdemand[i].pi * model._demand[i] for i in range(model._nstore))
                
                model.cbLazy(model._maxshipcost >= lazycut)
                
                model._iter += 1
        else:
            model.terminate()

class WareHouse:
    def __init__(self):
        # initialize data
        self.nwarehouse = 0
        self.nstore = 0
        self.supply = []
        self.demand = []
        self.fixcost = []
        self.varcost = []

        # initialize variables and constraints
        self.vmbuild = []
        self.vship = []
        self.csupply = []
        self.cdemand = []
    
    def read(self, filename):
        # input data
        with open(filename, "r") as data:
            self.nwarehouse = int(data.readline())
            self.nstore = int(data.readline())

            column = data.readline().split()
            for i in range(self.nwarehouse):
                self.supply.append(float(column[i]))
            
            column = data.readline().split()
            for i in range(self.nstore):
                self.demand.append(float(column[i]))
            
            column = data.readline().split()
            for i in range(self.nwarehouse):
                self.fixcost.append(float(column[i]))
            
            for i in range(self.nwarehouse):
                column = data.readline().split()
                lvarcost = []
                for j in range(self.nstore):
                    lvarcost.append(float(column[j]))
                self.varcost.append(lvarcost)
    
    def build(self):
        try:
            # define models for 'master' and 'sub'
            self.master = GRBPY.Model("master")
            self.sub = GRBPY.Model("sub")

            # disable log information
            self.master.setParam("OutputFlag", 0)
            self.sub.setParam("OutputFlag", 0)

            # use lazy constraints
            self.master.setParam("LazyConstraints", 1)

            # disable presolving in subproblem
            self.sub.setParam("Presolve", 0)

            # required to obtain farkas dual
            self.sub.setParam("InfUnbdInfo", 1)

            # use dual simplex
            self.sub.setParam("Method", 1)

            # construct master problem
            for i in range(self.nwarehouse):
                self.vmbuild.append(self.master.addVar(0.0, 1.0, 0.0, GRBPY.GRB.BINARY))
            
            self.maxshipcost = self.master.addVar(0.0, GRBPY.GRB.INFINITY, 0.0, GRBPY.GRB.CONTINUOUS)

            self.master.setObjective(GRBPY.quicksum(self.fixcost[i] * self.vmbuild[i] \
                                                    for i in range(self.nwarehouse)) + \
                                                    self.maxshipcost, GRBPY.GRB.MINIMIZE)

            # construct subproblem
            for i in range(self.nwarehouse):
                lvship = []
                for j in range(self.nstore):
                    lvship.append(self.sub.addVar(0.0, GRBPY.GRB.INFINITY, 0.0, GRBPY.GRB.CONTINUOUS))
                self.vship.append(lvship)
            
            for i in range(self.nwarehouse):
                self.csupply.append(self.sub.addConstr(GRBPY.quicksum(self.vship[i][j] \
                                                                      for j in range(self.nstore)) \
                                                                      <= self.supply[i] * 1.0))
            
            for j in range(self.nstore):
                self.cdemand.append(self.sub.addConstr(GRBPY.quicksum(self.vship[i][j] \
                                                                      for i in range(self.nwarehouse)) \
                                                                      == self.demand[j]))
            
            self.sub.setObjective(GRBPY.quicksum(self.varcost[i][j] * self.vship[i][j] \
                                                 for i in range(self.nwarehouse) \
                                                 for j in range(self.nstore)), GRBPY.GRB.MINIMIZE)
        except GRBPY.GurobiError as e:
            print('Error code' + str(e.errno) + ': ' + str(e))
        except AttributeError as e:
            print('Encountered an attribute error: ' + str(e))
    
    def solve(self):
        # build 'master' and 'sub'
        self.build()

        # register callback
        self.master._iter = 0
        self.master._nwarehouse = self.nwarehouse
        self.master._nstore = self.nstore
        self.master._supply = self.supply
        self.master._demand = self.demand
        
        self.master._csupply = self.csupply
        self.master._cdemand = self.cdemand
        self.master._vmbuild = self.vmbuild
        self.master._maxshipcost = self.maxshipcost

        self.master._sub = self.sub

        # optimize master problem
        print("               *** Benders Decomposition Loop ***               ")
        self.master.optimize(cbwarehouse)
        print("                        *** End Loop ***                        ")

        # it seems that 64-bit needs this extra work
        for i in range(self.nwarehouse):
            self.csupply[i].rhs = self.vmbuild[i].x * self.supply[i]
        
        self.sub.optimize()

    def report(self):
        print("               *** Summary Report ***               ")
        print("Objective: %.6f" % self.master.objval)
        print("Variables:")
        for i in range(self.nwarehouse):
            if abs(self.vmbuild[i].x) > 1e-6:
                print("  Build[%d] = %.0f" % (i, self.vmbuild[i].x))
        
        for i in range(self.nwarehouse):
            for j in range(self.nstore):
                if abs(self.vship[i][j].x) > 1e-6:
                    print("  Ship[%d][%d] = %.6f" % (i, j, self.vship[i][j].x))

if __name__ == "__main__":
    warehouse = WareHouse()
    warehouse.read("warehouse.dat")
    warehouse.solve()
    warehouse.report()