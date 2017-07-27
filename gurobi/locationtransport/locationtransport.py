from __future__ import division, print_function

import gurobipy as GRBPY


class LocationTransport:
    def __init__(self, name=None):
        # initialize data
        self.buildlimit = 0
        self.ncites = 0
        self.supply = []
        self.demand = []
        self.shipcost = []

        # initialize parameters
        self.iterlimit = 100
        self.samelimit = 3
        self.steplog = []
        self.scalelog = []
        self.xLBlog = []
        self.xUBlog = []

        if name is not None:
            self.name = name
        else:
            self.name = "demo"
        
        self.vship = []
        self.vbuild = []
        self.crelax = []
    
    def read(self, filename):
        with open(filename, "r") as data:
            self.buildlimit = int(data.readline())
            self.ncites = int(data.readline())
            
            column = data.readline().split()
            for i in range(self.ncites):
                self.supply.append(float(column[i]))
            
            column = data.readline().split()
            for i in range(self.ncites):
                self.demand.append(float(column[i]))
            
            for i in range(self.ncites):
                column = data.readline().split()
                lshipcost = []
                for j in range(self.ncites):
                    lshipcost.append(float(column[j]))
                self.shipcost.append(lshipcost)

    def build(self):
        try:
            self.mtrans = GRBPY.Model(self.name)

            # discard output information
            self.mtrans.setParam("OutputFlag", 0)

            # construct model
            for i in range(self.ncites):
                shipvar = []
                for j in range(self.ncites):
                    shipvar.append(self.mtrans.addVar(0.0, self.demand[j], 0.0, GRBPY.GRB.INTEGER))
                self.vship.append(shipvar)
            
            for i in range(self.ncites):
                self.vbuild.append(self.mtrans.addVar(0.0, 1.0, 0.0, GRBPY.GRB.BINARY))
            
            for i in range(self.ncites):
                self.mtrans.addConstr(GRBPY.quicksum(self.vship[i][j] for j in range(self.ncites)) \
                                      <= self.supply[i] * self.vbuild[i])
            
            self.mtrans.addConstr(GRBPY.quicksum(self.vbuild[i] for i in range(self.ncites)) \
                                  <= self.buildlimit)

            for j in range(self.ncites):
                self.crelax.append(self.mtrans.addConstr(GRBPY.quicksum(self.vship[i][j] \
                                   for i in range(self.ncites)) \
                                   >= self.demand[j]))
            
            self.mtrans.setObjective(GRBPY.quicksum(self.vship[i][j] * self.shipcost[i][j] \
                                     for i in range(self.ncites) \
                                     for j in range(self.ncites)), GRBPY.GRB.MINIMIZE)

            # update is necessary
            self.mtrans.update()
        except GRBPY.GurobiError as e:
            print('Error code' + str(e.errno) + ': ' + str(e))
        except AttributeError as e:
            print('Encountered an attribute error: ' + str(e))

    def solve(self):
        # 'Lagrangian Relaxation' parameters
        same = 0
        norm = 0.0
        step = 0.0
        scale = 1.0
        xLB = 0.0
        xUB = 0.0
        xlambda = [0.0] * self.ncites
        slack = [0.0] * self.ncites

        # build model
        self.build()

        # initial 'xLB'
        xLB = self.relaxUB(self.mtrans)
        
        # initial 'xUB'
        for i in range(self.ncites):
            xUB += max(self.shipcost[i])

        # temporary linear expression
        obj_shipcost = GRBPY.quicksum(self.vship[i][j] * self.shipcost[i][j] \
                                      for i in range(self.ncites) \
                                      for j in range(self.ncites))

        # sentinel flag
        lbmodel = 0

        # main 'Lagrangian Relaxation' loop
        for iter in range(self.iterlimit):
            # solve lower bound
            if lbmodel == 0:
                lbmodel = 1
                
                lenrelax = len(self.crelax)
                for i in range(lenrelax):
                    self.mtrans.remove(self.crelax[i])
                self.crelax = []
            
            obj_lagrangian = GRBPY.quicksum(xlambda[j] * (self.demand[j] - \
                                            GRBPY.quicksum(self.vship[i][j] for i in range(self.ncites))) \
                                            for j in range(self.ncites))
            self.mtrans.setObjective(obj_shipcost + obj_lagrangian, GRBPY.GRB.MINIMIZE)

            # 'LB' model
            self.mtrans.optimize()

            # calculate 'slack'
            for j in range(self.ncites):
                slack[j] = sum(self.vship[i][j].x for i in range(self.ncites)) - self.demand[j]

            # improve lower bound
            if self.mtrans.objval > xLB + 1e-6:
                xLB = self.mtrans.objval
                same = 0
            else:
                same += 1
            
            # update 'scale' if no improvement in 'samelimit' iteration
            if same == self.samelimit:
                scale /= 2.0
                same = 0
            
            # calculate 'norm'
            norm = sum(slack[i]**2.0 for i in range(self.ncites))

            # update 'step'
            step = scale * (xUB - self.mtrans.objval) / norm

            # update 'lambda'
            for i in range(self.ncites):
                if xlambda[i] > step * slack[i]:
                    xlambda[i] -= step * slack[i]
                else:
                    xlambda[i] = 0.0
            
            # solve upper bound
            sumsbval = sum(self.supply[i] * self.vbuild[i].x for i in range(self.ncites))
            sumdemand = sum(self.demand)

            if sumsbval - sumdemand >= 1e-6:
                lbmodel = 0
                
                for j in range(self.ncites):
                    self.crelax.append(self.mtrans.addConstr(GRBPY.quicksum(self.vship[i][j] \
                                                             for i in range(self.ncites)) >= self.demand[j]))
                
                # retrieve solution from LB model and fix it
                for i in range(self.ncites):
                    self.vbuild[i].lb = self.vbuild[i].x
                    self.vbuild[i].ub = self.vbuild[i].x
                
                self.mtrans.setObjective(obj_shipcost, GRBPY.GRB.MINIMIZE)

                self.mtrans.optimize()

                xUB = min(xUB, self.mtrans.objval)

                # reset to initial bound
                for i in range(self.ncites):
                    self.vbuild[i].lb = 0.0
                    self.vbuild[i].ub = 1.0
            
            # update 'xLBlog', 'xUBlog', 'steplog', 'scalelog'
            self.xLBlog.append(xLB)
            self.xUBlog.append(xUB)
            self.steplog.append(step)
            self.scalelog.append(scale)

    def relaxUB(self, mtrans):
        mrelax = mtrans.relax()

        mrelax.optimize()

        return mrelax.objval

    def report(self):
        print("\n               *** Summary Report ***               \n")
        print("  Iter        LB               UB          scale        step")

        for i in range(len(self.xLBlog)):
            print("  %3d    %12.6f    %12.6f    %8.6f    %8.6f" \
                  % (i, self.xLBlog[i], self.xUBlog[i], self.scalelog[i], self.steplog[i]))

if __name__ == "__main__":
    loctrans = LocationTransport()
    loctrans.read("loctrans.dat")
    loctrans.solve()
    loctrans.report()