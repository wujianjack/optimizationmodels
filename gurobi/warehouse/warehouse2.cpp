#include "gurobi_c++.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace std;

class warehousebenders: public GRBCallback {
public:
    // models
    GRBModel master;
    GRBModel sub;
    
    // variables in master problem
    GRBVar *mbuild;
    GRBVar maxshipcost;
    
    // variables and constraints in subproblem
    GRBVar **ship;
    GRBConstr *consupply;
    GRBConstr *condemand;
    
    // data to build model
    size_t nwarehouse;
    size_t nstore;
    
    double *supply;
    double *demand;
    double *fixcost;
    
    double **varcost;
    
    // lazycut to be added
    GRBLinExpr lazycut;
    
    // iteration count
    size_t iter;
    
    warehousebenders(GRBEnv &env) : iter(0), master(GRBModel(env)), sub(GRBModel(env)) {}
    
    ~warehousebenders() {
        // deallocate heap memory
        for (size_t i = 0; i < nwarehouse; ++i)
            delete [] ship[i];
        delete [] ship;
        
        delete [] mbuild;
        delete [] consupply;
        delete [] condemand;
        
        for (size_t i = 0; i < nwarehouse; ++i)
            delete [] varcost[i];
        delete [] varcost;
        
        delete [] fixcost;
        delete [] demand;
        delete [] supply;
    }
    
    void read(const char *name) {
        // input data
        std::ifstream data;
        
        data.open(name, std::ifstream::in);
        
        data >> nwarehouse;
        data >> nstore;
        
        supply = new double [nwarehouse];
        demand = new double [nstore];
        fixcost = new double [nwarehouse];
        
        varcost = new double *[nwarehouse];
        for (size_t i = 0; i < nwarehouse; ++i)
            varcost[i] = new double [nstore];
        
        for (size_t i = 0; i < nwarehouse; ++i)
            data >> supply[i];
        
        for (size_t i = 0; i < nstore; ++i)
            data >> demand[i];
        
        for (size_t i = 0; i < nwarehouse; ++i)
            data >> fixcost[i];
        
        for (size_t i = 0; i < nwarehouse; ++i)
            for (size_t j = 0; j < nstore; ++j)
                data >> varcost[i][j];
        
        data.close();
        // end data
    }

    void build() {
        // variables in master problem
        mbuild = new GRBVar [nwarehouse];
    
        // variables and constraints in subproblem
        ship = new GRBVar *[nwarehouse];
        for (size_t i = 0; i < nwarehouse; ++i)
            ship[i] = new GRBVar [nstore];
        
        consupply = new GRBConstr [nwarehouse];
        condemand = new GRBConstr [nstore];
        
        // disable log information
        master.set(GRB_IntParam_OutputFlag, 0);
        sub.set(GRB_IntParam_OutputFlag, 0);
        
        // use lazy constraints
        master.set(GRB_IntParam_LazyConstraints, 1);
        
        // disable presolving in subproblem
        sub.set(GRB_IntParam_Presolve, 0);
        
        // required to obtain farkas dual
        sub.set(GRB_IntParam_InfUnbdInfo, 1);
        
        // use dual simplex
        sub.set(GRB_IntParam_Method, 1);
        
        // construct master problem
        for (size_t i = 0; i < nwarehouse; ++i)
            mbuild[i] = master.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mbuild_" + to_string(i));
        
        maxshipcost = master.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "maxshipcost");
        
        GRBLinExpr obj_totalcost = 0.0;
        for (size_t i = 0; i < nwarehouse; ++i)
            obj_totalcost += fixcost[i] * mbuild[i];
        
        master.setObjective(obj_totalcost + maxshipcost, GRB_MINIMIZE);
        
        // construct subproblem
        for (size_t i = 0; i < nwarehouse; ++i)
            for (size_t j = 0; j < nstore; ++j)
                ship[i][j] = sub.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "ship_" + to_string(i) + "_" + to_string(j));
        
        GRBLinExpr con_supply = 0.0;
        for (size_t i = 0; i < nwarehouse; ++i) {
            for (size_t j = 0; j < nstore; ++j)
                con_supply += ship[i][j];
            
            consupply[i] = sub.addConstr(con_supply <= supply[i] * 1.0, "supply_" + to_string(i));
            con_supply = 0.0;
        }
        
        GRBLinExpr con_demand = 0.0;
        for (size_t j = 0; j < nstore; ++j) {
            for (size_t i = 0; i < nwarehouse; ++i)
                con_demand += ship[i][j];
            
            condemand[j] = sub.addConstr(con_demand == demand[j], "demand_" + to_string(j));
            con_demand = 0.0;
        }
        
        GRBLinExpr obj_shipcost = 0.0;
        for (size_t i = 0; i < nwarehouse; ++i)
            for (size_t j = 0; j < nstore; ++j)
                obj_shipcost += varcost[i][j] * ship[i][j];
        
        sub.setObjective(obj_shipcost, GRB_MINIMIZE);   
    }

    void run() {
        // build master and sub
        build();
        
        // register callback
        master.setCallback(this);
    
        // optimize master
        cout << "               *** Benders Decomposition Loop ***               " << endl;
        master.optimize();
        cout << "               *** End Loop ***               " << endl;
        
        // it seems that 64-bit needs this extra work
        for (size_t i = 0; i < nwarehouse; ++i)
            consupply[i].set(GRB_DoubleAttr_RHS, mbuild[i].get(GRB_DoubleAttr_X) * supply[i]);
        
        sub.optimize();
    }
    
    void solution() {
        // display solution
		cout << endl;
        cout << "              *** Summary Report ***               " << endl;
		printf("Objective: %.6f\n", master.get(GRB_DoubleAttr_ObjVal));
		cout << endl;
		cout << "Variables:: " << endl;
        for (size_t i = 0; i < nwarehouse; ++i) {
			if (fabs(mbuild[i].get(GRB_DoubleAttr_X)) > 1e-6)
				printf("  Build[%d] = %.0f\n", i, mbuild[i].get(GRB_DoubleAttr_X));
		}
        cout << endl;
               
        for (size_t i = 0; i < nwarehouse; ++i) {
            for (size_t j = 0; j < nstore; ++j) {
				if (fabs(ship[i][j].get(GRB_DoubleAttr_X)) > 1e-6)
                    printf("  Ship[%d][%d] = %.6f\n", i, j, ship[i][j].get(GRB_DoubleAttr_X));
			}
        }
		cout << endl;
    }

protected:
    void callback() {
        try {
            if (where == GRB_CB_MIPSOL) {
                if (iter >= 1) {
                    for (size_t i = 0; i < nwarehouse; ++i)
                        consupply[i].set(GRB_DoubleAttr_RHS, getSolution(mbuild[i]) * supply[i]);
                }
                
                sub.optimize();
                
                if (sub.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
                    cout << "Iteration: " << iter << endl;
                    cout << "Adding feasibility cut..." << endl;
                    cout << endl;
                    
                    lazycut = 0.0;
                    
                    for (size_t i = 0; i < nwarehouse; ++i)
                        lazycut += consupply[i].get(GRB_DoubleAttr_FarkasDual) * supply[i] * mbuild[i];
                    
                    for (size_t i = 0; i < nstore; ++i)
                        lazycut += condemand[i].get(GRB_DoubleAttr_FarkasDual) * demand[i];
                    
                    addLazy(lazycut >= 0);
                    
                    ++iter;
                }
                else if (sub.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
                    if (sub.get(GRB_DoubleAttr_ObjVal) > getSolution(maxshipcost) + 1e-6) {
                        cout << "Iteration: " << iter << endl;
                        cout << "Adding optimality cut..." << endl;
                        cout << endl;
                        
                        lazycut = 0.0;
                        
                        for (size_t i = 0; i < nwarehouse; ++i)
                            lazycut += consupply[i].get(GRB_DoubleAttr_Pi) * supply[i] * mbuild[i];
                    
                        for (size_t i = 0; i < nstore; ++i)
                            lazycut += condemand[i].get(GRB_DoubleAttr_Pi) * demand[i];
                        
                        addLazy(maxshipcost >= lazycut);
                        
                        ++iter;
                    }
                }
                else
                    abort();
            }
        }
        catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
        catch (...) {
            cout << "Error during callback" << endl;
        }
    }
};

int main(int argc, char *argv[]) {
    GRBEnv env = GRBEnv();
    
    warehousebenders warehouse = warehousebenders(env);
    
    warehouse.read("warehouse.dat");
    
    warehouse.run();
    
    warehouse.solution();
    
    return 0;
}