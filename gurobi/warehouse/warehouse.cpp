#include "gurobi_c++.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
    try {
        // input data
        ifstream data("warehouse.dat");
        
        size_t nwarehouse = 0;
        size_t nstore = 0;
        
        data >> nwarehouse;
        data >> nstore;
        
        double *supply = new double [nwarehouse];
        double *demand = new double [nstore];
        double *fixcost = new double [nwarehouse];
        
        double **varcost = new double *[nwarehouse];
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
        
        // environment and models
        GRBEnv env = GRBEnv();
        GRBModel master = GRBModel(env);
        GRBModel sub = GRBModel(env);
        
        // variables and parameters in master problem
        GRBVar *mbuild = new GRBVar [nwarehouse];
        GRBVar maxshipcost;
        
        // variables and constraints in subproblem
        GRBVar **ship = new GRBVar *[nwarehouse];
        for (size_t i = 0; i < nwarehouse; ++i)
            ship[i] = new GRBVar [nstore];
        
        GRBVar *sbuild = new GRBVar [nwarehouse];
        
        GRBConstr *consupply = new GRBConstr [nwarehouse];
        GRBConstr *condemand = new GRBConstr [nstore];
        
        // disable log information
        master.set(GRB_IntParam_OutputFlag, 0);
        sub.set(GRB_IntParam_OutputFlag, 0);
        
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
        
        for (size_t i = 0; i < nwarehouse; ++i)
            sbuild[i] = sub.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "sbuild_" + to_string(i));
		
        GRBLinExpr con_supply = 0.0;
        for (size_t i = 0; i < nwarehouse; ++i) {
            for (size_t j = 0; j < nstore; ++j)
                con_supply += ship[i][j];
            
            consupply[i] = sub.addConstr(con_supply <= supply[i] * sbuild[i], "supply_" + to_string(i));
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
        
        // run master once
        master.optimize();
        
        // temporary linear expression to be added
        GRBLinExpr con_lazycut = 0.0;
        
        // set iteration limit
        size_t iterlimit = 100;
		
        // set 'sbuild' value for initial subproblem
        for (size_t i = 0; i < nwarehouse; ++i) {
            sbuild[i].set(GRB_DoubleAttr_LB, 1.0);
            sbuild[i].set(GRB_DoubleAttr_UB, 1.0);
        }
        
        // main benders loop
        cout << "               *** Benders Decomposition Loop ***               " << endl;
        for (size_t iter = 0; iter < iterlimit; ++iter) {
            cout << "Iteration: " << iter << endl;
            
            sub.optimize();
            
            if (sub.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
                cout << "Adding feasibility cut..." << endl;
                cout << endl;
				
                con_lazycut = 0.0;
                
                for (size_t i = 0; i < nwarehouse; ++i)
                    con_lazycut += consupply[i].get(GRB_DoubleAttr_FarkasDual) * supply[i] * mbuild[i];
                    
                for (size_t j = 0; j < nstore; ++j)
                    con_lazycut += condemand[j].get(GRB_DoubleAttr_FarkasDual) * demand[j];
                
                master.addConstr(con_lazycut >= 0);
            }
            else if (sub.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
                if (sub.get(GRB_DoubleAttr_ObjVal) > maxshipcost.get(GRB_DoubleAttr_X) + 1e-6) {
                    cout << "Adding optimality cut..." << endl;
                    cout << endl;
					
                    con_lazycut = 0.0;
                    
                    for (size_t i = 0; i < nwarehouse; ++i)
                        con_lazycut += consupply[i].get(GRB_DoubleAttr_Pi) * supply[i] * mbuild[i];
                    
                    for (size_t j = 0; j < nstore; ++j)
                        con_lazycut += condemand[j].get(GRB_DoubleAttr_Pi) * demand[j];
                
                    master.addConstr(maxshipcost >= con_lazycut);
                }
                else
                    break;
            }
            else
                break;
            
            master.optimize();
            
            for (size_t i = 0; i < nwarehouse; ++i) {
                sbuild[i].set(GRB_DoubleAttr_LB, mbuild[i].get(GRB_DoubleAttr_X));
                sbuild[i].set(GRB_DoubleAttr_UB, mbuild[i].get(GRB_DoubleAttr_X));
            }
        }
        cout << "               *** End Loop ***               " << endl;
        
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
        delete [] sbuild;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch(...) {
        cout << "Exception during optimization" << endl;
    }
    
    return 0;
}