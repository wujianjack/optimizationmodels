#include "gurobi_c++.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
    try {
        // input data
        ifstream data("multicommodity.dat");
        
        size_t norig = 0;
        size_t ndest = 0;
        size_t nprod = 0;
        
        data >> norig;
        data >> ndest;
        data >> nprod;
        
        double **supply = new double *[norig];
        for (size_t i = 0; i < norig; ++i)
            supply[i] = new double [nprod];
        
        double **demand = new double *[ndest];
        for (size_t j = 0; j < ndest; ++j)
            demand[j] = new double [nprod];
        
        double **limit = new double *[norig];
        for (size_t i = 0; i < norig; ++i)
            limit[i] = new double [ndest];
        
        double ***cost = new double **[norig];
        for (size_t i = 0; i < norig; ++i) {
            cost[i] = new double *[ndest];
            
            for (size_t j = 0; j < ndest; ++j)
                cost[i][j] = new double [nprod];
        }
        
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < nprod; ++j)
                data >> supply[i][j];
        
        for (size_t j = 0; j < ndest; ++j)
            for (size_t k = 0; k < nprod; ++k)
                data >> demand[j][k];
        
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                data >> limit[i][j];
        
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                for (size_t k = 0; k < nprod; ++k)
                    data >> cost[i][j][k];
        
        data.close();
        
        // environment and models
        GRBEnv env = GRBEnv();
        GRBModel master = GRBModel(env);
        GRBModel sub = GRBModel(env);
        
        // disable log information
        master.set(GRB_IntParam_OutputFlag, 0);
        sub.set(GRB_IntParam_OutputFlag, 0);
        
        // iteration count
        size_t iter = 0;
        
        // variables, constraints and parameters in mater problem
        GRBConstr **multi = new GRBConstr *[norig];
        for (size_t i = 0; i < norig; ++i)
            multi[i] = new GRBConstr [ndest];
        
        GRBConstr *convex = NULL;
        
        vector<GRBVar> weight;
        GRBVar excess;
        
        GRBLinExpr obj_masteri = 0.0;
        GRBLinExpr obj_masterii = 0.0;
        GRBLinExpr obj_masteriii = 0.0;
        
        GRBColumn col;
        
        double **xtrans = NULL;
        double sum_xtrans = 0.0;
        double sum_costxtrans = 0.0;
        vector<double> propcost;
        
        // variables and parameters in sub problem
        GRBVar ***trans = new GRBVar **[norig];
        for (size_t i = 0; i < norig; ++i) {
            trans[i] = new GRBVar *[ndest];
            
            for (size_t j = 0; j < ndest; ++j)
                trans[i][j] = new GRBVar [nprod];
        }
        
        GRBLinExpr obj_subi = 0.0;
        GRBLinExpr obj_subii = 0.0;
        
        GRBLinExpr con_supply = 0.0;
        GRBLinExpr con_demand = 0.0;
        
        GRBLinExpr con_optmulti = 0.0;
        
        double priceconvex = 1.0;
        double **price = new double *[norig];
        for (size_t i = 0; i < norig; ++i)
            price[i] = new double [ndest];
        
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                price[i][j] = 0.0;
        
        double **optship = new double *[norig];
        for (size_t i = 0; i < norig; ++i)
            optship[i] = new double [ndest];
        
        // add variable 'excess' for master problem
        excess = master.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, string("excess"));
        
        // add temporary constraints 'multi' for master problem
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                multi[i][j] = master.addConstr(-excess <= limit[i][j], "multi" + to_string(i) + to_string(j));
        
        // add temporary constraint 'convex' for master problem
        convex = master.addConstrs(1);
        convex[0].set(GRB_CharAttr_Sense, GRB_EQUAL);
        convex[0].set(GRB_DoubleAttr_RHS, 1.0);
        convex[0].set(GRB_StringAttr_ConstrName, string("convex"));
        
        // set objective for master problem in 'Phase I'
        obj_masteri = excess;
        
        // add variables 'trans' for subproblem
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                for (size_t k = 0; k < nprod; ++k)
                    trans[i][j][k] = sub.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "trans" + to_string(i) + to_string(j) + to_string(k));
        
        // add constraints 'supply' for subproblem
        for (size_t i = 0; i < norig; ++i) {
            for (size_t k = 0; k < nprod; ++k) {
                for (size_t j = 0; j < ndest; ++j)
                    con_supply += trans[i][j][k];
                
                sub.addConstr(con_supply == supply[i][k], "supply" + to_string(i) + to_string(k));
                con_supply = 0.0;
            }
        }
        
        // add constraints 'demand' for subproblem
        for (size_t j = 0; j < ndest; ++j) {
            for (size_t k = 0; k < nprod; ++k) {
                for (size_t i = 0; i < norig; ++i)
                    con_demand += trans[i][j][k];
                
                sub.addConstr(con_demand == demand[j][k], "demand" + to_string(j) + to_string(k));
                con_demand = 0.0;
            }
        }
        
        // set objective for subproblem in 'Phase I'
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                for (size_t k = 0; k < nprod; ++k)
                    obj_subi -= price[i][j] * trans[i][j][k];
        
        obj_subi -= priceconvex;
        
        // dantzig-wolfe decomposition
        cout << "               *** Dantzig-Wolfe Decomposition ***               " << endl;
        cout << endl;
        
        // set objective for subproblem in 'Phase I'
        sub.setObjective(obj_subi, GRB_MINIMIZE);
        
        // set objective for master problem in 'Phase I'
        master.setObjective(obj_masteri, GRB_MINIMIZE);
        
        // 'Phase I' of dantzig-wolfe decomposition
        cout << "Phase I: " << endl;
        for (;;) {
            cout << "Iteration: " << iter << endl;
            
            // solve subproblem in 'Phase I'
            sub.optimize();
            
            if (sub.get(GRB_DoubleAttr_ObjVal) >= -1e-6) {
                cout << "No feasible solution..." << endl;
                break;
            }
            else {
                ++iter;
                
                // calculate parameters for master problem
                xtrans = new double *[norig];
                for (size_t i = 0; i < norig; ++i)
                    xtrans[i] = new double [ndest];
                
                sum_xtrans = 0.0;
                for (size_t i = 0; i < norig; ++i) {
                    for (size_t j = 0; j < ndest; ++j) {
                        for (size_t k = 0; k < nprod; ++k)
                            sum_xtrans += trans[i][j][k].get(GRB_DoubleAttr_X);
                        
                        xtrans[i][j] = sum_xtrans;
                        sum_xtrans = 0.0;
                    }
                }
                
                sum_costxtrans = 0.0;
                for (size_t i = 0; i < norig; ++i)
                    for (size_t j = 0; j < ndest; ++j)
                        for (size_t k = 0; k < nprod; ++k)
                            sum_costxtrans += cost[i][j][k] * trans[i][j][k].get(GRB_DoubleAttr_X);
                
                propcost.push_back(sum_costxtrans);
                
                // update constraints in master problem
                col = GRBColumn();
                for (size_t i = 0; i < norig; ++i)
                    for (size_t j = 0; j < ndest; ++j)
                        col.addTerm(xtrans[i][j], multi[i][j]);
                
                col.addTerm(1.0, convex[0]);
                
                // add variable 'weight'
                weight.push_back(master.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, col, "weight" + to_string(iter)));
                
                // deallocate heap memory
                for (size_t i = 0; i < norig; ++i)
                    delete [] xtrans[i];
                delete [] xtrans;
                
                // solve master problem in 'Phase I'
                master.optimize();
                
                // update price
                if (master.get(GRB_DoubleAttr_ObjVal) <= 1e-6)
                    break;
                else {
                    for (size_t i = 0; i < norig; ++i)
                        for (size_t j = 0; j < ndest; ++j)
                            price[i][j] = multi[i][j].get(GRB_DoubleAttr_Pi);
                    
                    priceconvex = convex[0].get(GRB_DoubleAttr_Pi);
                }
                
                // reset objective for subproblem in 'Phase I'
                obj_subi = 0.0;
                for (size_t i = 0; i < norig; ++i)
                    for (size_t j = 0; j < ndest; ++j)
                        for (size_t k = 0; k < nprod; ++k)
                            obj_subi -= price[i][j] * trans[i][j][k];
                
                obj_subi -= priceconvex;
                
                sub.setObjective(obj_subi, GRB_MINIMIZE);
            }
        }
        
        // setting up for 'Phase II'
        cout << "Setting up for Phase II..." << endl;
        
        // set objective for master problem in 'Phase II'
        for (size_t i = 0; i < weight.size(); ++i)
            obj_masterii += propcost[i] * weight[i];
        
        master.setObjective(obj_masterii, GRB_MINIMIZE);
        
        // fix variable 'excess'
        double xexcess = excess.get(GRB_DoubleAttr_X);
        
        excess.set(GRB_DoubleAttr_LB, xexcess);
        excess.set(GRB_DoubleAttr_UB, xexcess);
        
        // solve master problem in 'Phase II'
        master.optimize();
        
        // update price
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                price[i][j] = multi[i][j].get(GRB_DoubleAttr_Pi);
        
        priceconvex = convex[0].get(GRB_DoubleAttr_Pi);
        
        // set objective for subproblem in 'Phase II'
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                for (size_t k = 0; k < nprod; ++k)
                    obj_subii += (cost[i][j][k] - price[i][j]) * trans[i][j][k];
        
        obj_subii -= priceconvex;
        
        sub.setObjective(obj_subii, GRB_MINIMIZE);
        
        // increase iteration count
        ++iter;
        
        // 'Phase II' of dantzig-wolfe decomposition
        cout << "Phase II: " << endl;
        for (;;) {
            cout << "Iteration: " << iter << endl;
            
            // solve subproblem in 'Phase II'
            sub.optimize();
            
            if (sub.get(GRB_DoubleAttr_ObjVal) >= -1e-6) {
                cout << "Optimal solution..." << endl;
                break;
            }
            else {
                ++iter;
                
                // calculate parameters for master problem
                xtrans = new double *[norig];
                for (size_t i = 0; i < norig; ++i)
                    xtrans[i] = new double [ndest];
                
                sum_xtrans = 0.0;
                for (size_t i = 0; i < norig; ++i) {
                    for (size_t j = 0; j < ndest; ++j) {
                        for (size_t k = 0; k < nprod; ++k)
                            sum_xtrans += trans[i][j][k].get(GRB_DoubleAttr_X);
                        
                        xtrans[i][j] = sum_xtrans;
                        sum_xtrans = 0.0;
                    }
                }
                
                sum_costxtrans = 0.0;
                for (size_t i = 0; i < norig; ++i)
                    for (size_t j = 0; j < ndest; ++j)
                        for (size_t k = 0; k < nprod; ++k)
                            sum_costxtrans += cost[i][j][k] * trans[i][j][k].get(GRB_DoubleAttr_X);
                
                // update constraints in master problem
                col = GRBColumn();
                for (size_t i = 0; i < norig; ++i)
                    for (size_t j = 0; j < ndest; ++j)
                        col.addTerm(xtrans[i][j], multi[i][j]);
                
                col.addTerm(1.0, convex[0]);
                
                // add variable 'weight'
                weight.push_back(master.addVar(0.0, GRB_INFINITY, sum_costxtrans, GRB_CONTINUOUS, col, "weight" + to_string(iter)));
                
                // deallocate heap memory
                for (size_t i = 0; i < norig; ++i)
                    delete [] xtrans[i];
                delete [] xtrans;
                
                // solve master problem in 'Phase II'
                master.optimize();
                
                // update price
                for (size_t i = 0; i < norig; ++i)
                    for (size_t j = 0; j < ndest; ++j)
                        price[i][j] = multi[i][j].get(GRB_DoubleAttr_Pi);
                    
                priceconvex = convex[0].get(GRB_DoubleAttr_Pi);
                
                // set objective for subproblem in 'Phase II'
                obj_subii = 0.0;
                for (size_t i = 0; i < norig; ++i)
                    for (size_t j = 0; j < ndest; ++j)
                        for (size_t k = 0; k < nprod; ++k)
                            obj_subii += (cost[i][j][k] - price[i][j]) * trans[i][j][k];
                
                obj_subii -= priceconvex;
                
                sub.setObjective(obj_subii, GRB_MINIMIZE);
            }
        }

        // 'Phase III' of dantzig-wolfe decomposition
        cout << "Phase III: " << endl;
        
        // set objective for master problem in 'Phase III'
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                optship[i][j] = limit[i][j] + excess.get(GRB_DoubleAttr_X) - multi[i][j].get(GRB_DoubleAttr_Slack);
        
        for (size_t i = 0; i < norig; ++i)
            for (size_t j = 0; j < ndest; ++j)
                for (size_t k = 0; k < nprod; ++k)
                    obj_masteriii += cost[i][j][k] * trans[i][j][k];
        
        sub.setObjective(obj_masteriii, GRB_MINIMIZE);
        
        for (size_t i = 0; i < norig; ++i) {
            for (size_t j = 0; j < ndest; ++j) {
                for (size_t k = 0; k < nprod; ++k)
                    con_optmulti += trans[i][j][k];
                
                sub.addConstr(con_optmulti == optship[i][j], "optmulti" + to_string(i) + to_string(j));
                con_optmulti = 0.0;
            }
        }
        
        // solve master problem in 'Phase III'
        sub.optimize();
        
        // end loop
        cout << "               *** End Loop ***               " << endl;
        
        // report solution
        cout << "               *** Summary Report ***               " << endl;
        cout << "Objective: " << sub.get(GRB_DoubleAttr_ObjVal) << endl;
        cout << "Variables: " << endl;
        for (size_t i = 0; i < norig; ++i) {
            for (size_t j = 0; j < ndest; ++j) {
                for (size_t k = 0; k < nprod; ++k) {
                    if (fabs(trans[i][j][k].get(GRB_DoubleAttr_X)) >= 1e-6)
                        printf("  trans[%d][%d][%d] = %12.6f\n", i, j, k, trans[i][j][k].get(GRB_DoubleAttr_X));
                }
            }
        }
        
        // deallocate heap memory
        for (size_t i = 0; i < norig; ++i)
            delete [] supply[i];
        delete [] supply;
        
        for (size_t i = 0; i < ndest; ++i)
            delete [] demand[i];
        delete [] demand;
        
        for (size_t i = 0; i < norig; ++i)
            delete [] limit[i];
        delete [] limit;
        
        for (size_t i = 0; i < norig; ++i) {
            for (size_t j = 0; j < ndest; ++j)
                delete [] cost[i][j];
            delete [] cost[i];
        }
        delete [] cost;
        
        for (size_t i = 0; i < norig; ++i)
            delete [] multi[i];
        delete [] multi;
        
        for (size_t i = 0; i < norig; ++i) {
            for (size_t j = 0; j < ndest; ++j)
                delete [] trans[i][j];
            delete [] trans[i];
        }
        delete [] trans;
        
        for (size_t i = 0; i < norig; ++i)
            delete [] price[i];
        delete [] price;
        
        for (size_t i = 0; i < norig; ++i)
            delete [] optship[i];
        delete [] optship;
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