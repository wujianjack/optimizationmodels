from __future__ import division, print_function

import gurobipy as GRBPY


class Multicommodity:
    def __init__(self):
        # initialize data
        self.norig = 0
        self.ndest = 0
        self.nprod = 0
        self.supply = []
        self.demand = []
        self.limit = []
        self.cost = []

        # initialize variables and constraints
        self.cmulti = []
        self.vtrans = []
        self.vweight = []

        # initialize parameters
        self.iter = 0
        self.priceconvex = 1.0
        self.price = []
        self.propcost = []
    
    def read(self, filename):
        # input data
        with open(filename, "r") as data:
            self.norig = int(data.readline())
            self.ndest = int(data.readline())
            self.nprod = int(data.readline())

            for i in range(self.norig):
                column = data.readline().split()
                lsupply = []
                for k in range(self.nprod):
                    lsupply.append(float(column[k]))
                self.supply.append(lsupply)

            for j in range(self.ndest):
                column = data.readline().split()
                ldemand = []
                for k in range(self.nprod):
                    ldemand.append(float(column[k]))
                self.demand.append(ldemand)

            for i in range(self.norig):
                column = data.readline().split()
                llimit = []
                for j in range(self.ndest):
                    llimit.append(float(column[j]))
                self.limit.append(llimit)

            for i in range(self.norig):
                lcost = []
                for j in range(self.ndest):
                    column = data.readline().split()
                    llcost = []
                    for k in range(self.nprod):
                        llcost.append(float(column[k]))
                    lcost.append(llcost)
                self.cost.append(lcost)
    
    def build(self):
        try:
            # define models for 'master' and 'sub'
            self.master = GRBPY.Model("master")
            self.sub = GRBPY.Model("sub")

            # disable log information
            self.master.setParam("OutputFlag", 0)
            self.sub.setParam("OutputFlag", 0)

            # add variable 'excess' for master problem
            self.vexcess = self.master.addVar(0.0, GRBPY.GRB.INFINITY, 0.0, GRBPY.GRB.CONTINUOUS)

            # add temporary constraints 'multi' for master problem
            for i in range(self.norig):
                lmulti = []
                for j in range(self.ndest):
                    lmulti.append(self.master.addConstr(-self.vexcess <= self.limit[i][j]))
                self.cmulti.append(lmulti)
            
            # tricky approach to add temporary constraint 'convex' for master problem
            self.cconvex = self.master.addConstr(1e-6 * self.vexcess == 1.0)
            
            # add variables 'trans' for subproblem
            for i in range(self.norig):
                ltrans = []
                for j in range(self.ndest):
                    lltrans = []
                    for k in range(self.nprod):
                        lltrans.append(self.sub.addVar(0.0, GRBPY.GRB.INFINITY, 0.0, GRBPY.GRB.CONTINUOUS))
                    ltrans.append(lltrans)
                self.vtrans.append(ltrans)
            
            # add constraints 'supply' for subproblem
            for i in range(self.norig):
                for k in range(self.nprod):
                    self.sub.addConstr(GRBPY.quicksum(self.vtrans[i][j][k] for j in range(self.ndest)) \
                                                      == self.supply[i][k])
            
            # add constraints 'demand' for subproblem
            for j in range(self.ndest):
                for k in range(self.nprod):
                    self.sub.addConstr(GRBPY.quicksum(self.vtrans[i][j][k] for i in range(self.norig)) \
                                                      == self.demand[j][k])
        except GRBPY.GurobiError as e:
            print('Error code' + str(e.errno) + ': ' + str(e))
        except AttributeError as e:
            print('Encountered an attribute error: ' + str(e))
    
    def solvePhaseI(self):
        # initialize parameters
        for i in range(self.norig):
            lprice = [0.0] * self.ndest
            self.price.append(lprice)
        
        obj_masteri = self.vexcess
        obj_subi = -GRBPY.quicksum(self.price[i][j] * self.vtrans[i][j][k] \
                                   for i in range(self.norig) \
                                   for j in range(self.ndest) \
                                   for k in range(self.nprod)) - self.priceconvex

        # set objective for master problem in 'Phase I'
        self.master.setObjective(obj_masteri, GRBPY.GRB.MINIMIZE)

        # set objective for subproblem in 'Phase I'
        self.sub.setObjective(obj_subi, GRBPY.GRB.MINIMIZE)

        # ok, forget this
        self.master.chgCoeff(self.cconvex, self.vexcess, 0.0)

        # 'Phase I' of dantzig-wolfe decomposition
        print("Phase I: ")

        while True:
            print("Iteration: ", self.iter)

            # solve subproblem in 'Phase I'
            self.sub.optimize()

            if self.sub.objval >= -1e-6:
                print("No feasible solution...")
                break
            else:
                self.iter += 1

                # calculate parameters for master problem
                sum_costxtrans = sum(self.cost[i][j][k] * self.vtrans[i][j][k].x \
                                     for i in range(self.norig) \
                                     for j in range(self.ndest) \
                                     for k in range(self.nprod))

                self.propcost.append(sum_costxtrans)

                # update constraints in master problem
                col = GRBPY.Column()
                for i in range(self.norig):
                    for j in range(self.ndest):
                        col.addTerms(sum(self.vtrans[i][j][k].x for k in range(self.nprod)), self.cmulti[i][j])
                
                col.addTerms(1.0, self.cconvex)

                # add variable 'weight'
                self.vweight.append(self.master.addVar(0.0, GRBPY.GRB.INFINITY, 0.0, GRBPY.GRB.CONTINUOUS, "", col))

                # solve master problem in 'Phase I'
                self.master.optimize()

                # update price
                if self.master.objval <= 1e-6:
                    break
                else:
                    for i in range(self.norig):
                        for j in range(self.ndest):
                            self.price[i][j] = self.cmulti[i][j].pi
                    
                    self.priceconvex = self.cconvex.pi
                
                # reset objective for subproblem in 'Phase I'
                obj_subi = -GRBPY.quicksum(self.price[i][j] * self.vtrans[i][j][k] \
                                           for i in range(self.norig) \
                                           for j in range(self.ndest) \
                                           for k in range(self.nprod)) - self.priceconvex

                self.sub.setObjective(obj_subi, GRBPY.GRB.MINIMIZE)
    
    def setupPhaseII(self):
        # setting up for 'Phase II'
        print("Setting up for Phase II...")

        # set objective for master problem in 'Phase II'
        obj_masterii = GRBPY.quicksum(self.propcost[i] * self.vweight[i] for i in range(len(self.propcost)))

        self.master.setObjective(obj_masterii, GRBPY.GRB.MINIMIZE)

        # fix variable 'excess'
        self.vexcess.lb = self.vexcess.x
        self.vexcess.ub = self.vexcess.x

        # solve master problem in 'Phase II'
        self.master.optimize()

        # update price
        for i in range(self.norig):
            for j in range(self.ndest):
                self.price[i][j] = self.cmulti[i][j].pi
        
        self.priceconvex = self.cconvex.pi

        # set objective for subproblem in 'Phase II'
        obj_subii = GRBPY.quicksum((self.cost[i][j][k] - self.price[i][j]) * self.vtrans[i][j][k] \
                                   for i in range(self.norig) \
                                   for j in range(self.ndest) \
                                   for k in range(self.nprod)) - self.priceconvex
        
        self.sub.setObjective(obj_subii, GRBPY.GRB.MINIMIZE)

        # increase iteration count
        self.iter += 1

    def solvePhaseII(self):
        # 'Phase II' of dantzig-wolfe decomposition
        print("Phase II: ")

        while True:
            print("Iteration: ", self.iter)

            # solve subproblem in 'Phase II'
            self.sub.optimize()

            if self.sub.objval >= -1e-6:
                print("Optimal solution...")
                break
            else:
                self.iter += 1

                # calculate parameters for master problem
                sum_costxtrans = sum(self.cost[i][j][k] * self.vtrans[i][j][k].x \
                                     for i in range(self.norig) \
                                     for j in range(self.ndest) \
                                     for k in range(self.nprod))
                
                # update constraints in master problem
                col = GRBPY.Column()
                for i in range(self.norig):
                    for j in range(self.ndest):
                        col.addTerms(sum(self.vtrans[i][j][k].x for k in range(self.nprod)), self.cmulti[i][j])
                
                col.addTerms(1.0, self.cconvex)

                # add variable 'weight'
                self.vweight.append(self.master.addVar(0.0, GRBPY.GRB.INFINITY, sum_costxtrans, GRBPY.GRB.CONTINUOUS, "", col))

                # solve master problem in 'Phase II'
                self.master.optimize()

                # update price
                for i in range(self.norig):
                    for j in range(self.ndest):
                        self.price[i][j] = self.cmulti[i][j].pi
                
                self.priceconvex = self.cconvex.pi

                # reset objective for subproblem in 'Phase II'
                obj_subii = GRBPY.quicksum((self.cost[i][j][k] - self.price[i][j]) * self.vtrans[i][j][k] \
                                            for i in range(self.norig) \
                                            for j in range(self.ndest) \
                                            for k in range(self.nprod)) - self.priceconvex
                
                self.sub.setObjective(obj_subii, GRBPY.GRB.MINIMIZE)
    
    def solvePhaseIII(self):
        optship = []
        for i in range(self.norig):
            loptship = [0.0] * self.ndest
            optship.append(loptship)
        
        # 'Phase III' of dantzig-wolfe decomposition
        print("Phase III:")

        # set objective for master problem in 'Phase III'
        for i in range(self.norig):
            for j in range(self.ndest):
                optship[i][j] = self.limit[i][j] + self.vexcess.x - self.cmulti[i][j].slack
        
        obj_masteriii = GRBPY.quicksum(self.cost[i][j][k] * self.vtrans[i][j][k] \
                                       for i in range(self.norig) \
                                       for j in range(self.ndest) \
                                       for k in range(self.nprod))
        
        self.sub.setObjective(obj_masteriii, GRBPY.GRB.MINIMIZE)

        for i in range(self.norig):
            for j in range(self.ndest):
                self.sub.addConstr(GRBPY.quicksum(self.vtrans[i][j][k] for k in range(self.nprod)) == optship[i][j])
        
        # solve master problem in 'Phase III'
        self.sub.optimize()

    def solve(self):
        # build 'master' and 'sub'
        self.build()

        # dantzig-wolfe decomposition
        print("               *** Dantzig-Wolfe Decomposition ***               ")
        self.solvePhaseI()

        self.setupPhaseII()

        self.solvePhaseII()

        self.solvePhaseIII()
        print("                        *** End Loop ***                        ")
    
    def report(self):
        # report solution
        print("               *** Summary Report ***                           ")
        print("Objective: ", self.sub.objval)
        print("Variables: ")
        for i in range(self.norig):
            for j in range(self.ndest):
                for k in range(self.nprod):
                    if abs(self.vtrans[i][j][k].x) >= 1e-6:
                        print("  trans[%d][%d][%d] = %12.6f" % (i, j, k, self.vtrans[i][j][k].x))

if __name__ == "__main__":
    multicomm = Multicommodity()
    multicomm.read("multicommodity.dat")
    multicomm.solve()
    multicomm.report()