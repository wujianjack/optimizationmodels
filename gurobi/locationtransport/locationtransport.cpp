#include "gurobi_c++.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <algorithm>

#ifdef WIN32
    #define NOMINMAX 1
#endif

using namespace std;

double relaxUB(GRBModel &model);
double calculateNorm(double *slack, size_t ncites);
void reportLog(double *LBlog, double *UBlog, double *scalelog, double *steplog, size_t count);

int main(int argc, char *argv[]) {
    try {
        // Input data
        ifstream data("loctrans.dat");
        
        size_t buildlimit = 0;
        size_t ncites = 0;
        
        data >> buildlimit;
        data >> ncites;
        
        double *supply = new double [ncites];
        double *demand = new double [ncites];
        
        double **shipcost = new double *[ncites];
        for (size_t i = 0; i < ncites; ++i)
            shipcost[i] = new double [ncites];
        
        for (size_t i = 0; i < ncites; ++i)
            data >> supply[i];
        
        for (size_t i = 0; i < ncites; ++i)
            data >> demand[i];
        
        for (size_t i = 0; i < ncites; ++i)
            for (size_t j = 0; j < ncites; ++j)
                data >> shipcost[i][j];
        data.close();
        // End data
       
        // Define parameters
        size_t iterlimit = 200;
        size_t samelimit = 3;
        size_t same = 0;
        double norm = 0.0;
        double step = 0.0;
        double scale = 1.0;
        double LB = 0.0;
        double UB = 0.0;
        double *lambda = new double [ncites];
        double *slack = new double [ncites];
        double *steplog = new double [iterlimit];
        double *scalelog = new double [iterlimit];
        double *LBlog = new double [iterlimit];
        double *UBlog = new double [iterlimit];
        // End parameters
        
        // Initialize parameters
        for (size_t i = 0; i < iterlimit; ++i) {
            LBlog[i] = 0.0;
            UBlog[i] = 0.0;
        }
        
        for (size_t i = 0; i < ncites; ++i) {
            lambda[i] = 0.0;
            slack[i] = 0.0;
        }
        // End initialization
        
        GRBEnv env = GRBEnv();
        GRBModel trans = GRBModel(env);
        
        GRBVar **ship = new GRBVar *[ncites];
        for (size_t i = 0; i < ncites; ++i)
            ship[i] = new GRBVar [ncites];
        
        GRBVar *build = new GRBVar[ncites];
        
        for (size_t i = 0; i < ncites; ++i) {
            for (size_t j = 0; j < ncites; ++j)
                ship[i][j] = trans.addVar(0.0, demand[j], 0.0, GRB_INTEGER, "ship_" + to_string(i) + "_" + to_string(j));
        }
        
        for (size_t i = 0; i < ncites; ++i)
            build[i] = trans.addVar(0.0, 1.0, 0.0, GRB_BINARY, "build_" + to_string(i));
        
        GRBLinExpr con_supply = 0.0;
        for (size_t i = 0; i < ncites; ++i) {
            for (size_t j = 0; j < ncites; ++j)
                con_supply += ship[i][j];
            
            trans.addConstr(con_supply <= supply[i] * build[i], "supply_" + to_string(i));
            con_supply = 0.0;
        }
        
        GRBLinExpr con_limit = 0.0;
        for (size_t i = 0; i < ncites; ++i)
            con_limit += build[i];
        
        trans.addConstr(con_limit <= buildlimit, string("limit"));
        
        // update is necessary so far (Gurobi v7.5)
        trans.update();
        
        GRBModel transLB = GRBModel(trans);
        GRBModel transUB = GRBModel(trans);
        
        transLB.set(GRB_IntParam_OutputFlag, 0);
        transUB.set(GRB_IntParam_OutputFlag, 0);
        
        GRBLinExpr con_demand = 0.0;
        for (size_t j = 0; j < ncites; ++j) {
            for (size_t i = 0; i < ncites; ++i)
                con_demand += ship[i][j];
            
            transUB.addConstr(con_demand >= demand[j], "demand_" + to_string(j));
            con_demand = 0.0;
        }
        
        GRBLinExpr obj_shipcost = 0.0;
        for (size_t i = 0; i < ncites; ++i) {
            for (size_t j = 0; j < ncites; ++j)
                obj_shipcost += ship[i][j] * shipcost[i][j];
        }
        
        transUB.setObjective(obj_shipcost, GRB_MINIMIZE);
        
        // initial 'LB'
        LB = relaxUB(transUB);
        
        // initial 'UB'
        double *rowmax = new double [ncites];
        for (size_t i = 0; i < ncites; ++i)
            rowmax[i] = *max_element(shipcost[i], shipcost[i] + ncites);
    
        for (size_t i = 0; i < ncites; ++i)
            UB += rowmax[i];
        
        delete [] rowmax;
        
        GRBLinExpr obj_lagrange = 0.0;
        GRBLinExpr sum_ship = 0.0;
        double sumsbval = 0.0;
        double sumdemand = 0.0;
        // main lagrange relaxation loop
        cout << "               *** Larange Relaxation ***               " << endl;
        for (size_t i = 0; i < iterlimit; ++i) {  
            // solve lower bound
            obj_lagrange = 0.0;
            sum_ship = 0.0;
            for (size_t j = 0; j < ncites; ++j) {
                for (size_t i = 0; i < ncites; ++i) 
                    sum_ship += ship[i][j];
                
                obj_lagrange += lambda[j] * (demand[j] - sum_ship);
                sum_ship = 0.0;
            }
                
            transLB.setObjective(obj_shipcost + obj_lagrange, GRB_MINIMIZE);
            
            // LB model
            transLB.optimize();
            
            // calculate 'slack'
            for (size_t j = 0; j < ncites; ++j) {
                double sum_shipval = 0.0;
                GRBVar shipvar;
                
                for (size_t ii = 0; ii < ncites; ++ii) {
                    shipvar = transLB.getVarByName("ship_" + to_string(ii) + "_" + to_string(j));
                    sum_shipval += shipvar.get(GRB_DoubleAttr_X);
                }
                
                slack[j] = sum_shipval - demand[j];
            }
            
            // improve lower bound
            if (transLB.get(GRB_DoubleAttr_ObjVal) > LB + 1e-6) {
                LB = transLB.get(GRB_DoubleAttr_ObjVal);
                same = 0;
            }
            else
                ++same;
            
            // update 'scale' if no improvment in 'samelimit' iteration
            if (same == samelimit) {
                scale /= 2.0;
                same = 0;
            }
            
            // calculate 'norm'
            norm = calculateNorm(slack, ncites);
            
            // update 'step'
            step = scale * (UB - transLB.get(GRB_DoubleAttr_ObjVal)) / norm;
            
            // update 'lambda'
            for (size_t j = 0; j < ncites; ++j) {
                if (lambda[j] > (step * slack[j]))
                    lambda[j] -= step * slack[j];
                else
                    lambda[j] = 0.0;
            }
            
            // solve upper bound
            sumsbval = 0.0;
            for (size_t j = 0; j < ncites; ++j) {
                GRBVar buildvar;
                
                buildvar = transLB.getVarByName("build_" + to_string(j));
                sumsbval += supply[j] * buildvar.get(GRB_DoubleAttr_X);
            }
            
            sumdemand = 0.0;
            for (size_t j = 0; j < ncites; ++j)
                sumdemand += demand[j];
            
            if (sumsbval - sumdemand >= -1e-6) {
                // retrieve solution from LB model and fix it
                for (size_t j = 0; j < ncites; ++j) {
                    GRBVar buildvarl;
                    GRBVar buildvaru;
                    
                    buildvarl = transLB.getVarByName("build_" + to_string(j));
                    buildvaru = transUB.getVarByName("build_" + to_string(j));
                    buildvaru.set(GRB_DoubleAttr_LB, buildvarl.get(GRB_DoubleAttr_X));
                    buildvaru.set(GRB_DoubleAttr_UB, buildvarl.get(GRB_DoubleAttr_X));
                }
                
                transUB.optimize();
                
                UB = min(UB, transUB.get(GRB_DoubleAttr_ObjVal));
                
                // reset to initial bound
                for (size_t j = 0; j < ncites; ++j) {
                    GRBVar buildvaru;
                    
                    buildvaru = transUB.getVarByName("build_" + to_string(j));
                    buildvaru.set(GRB_DoubleAttr_LB, 0.0);
                    buildvaru.set(GRB_DoubleAttr_UB, 1.0);
                }
            }
            
            // update 'LBlog', 'UBlog', 'steplog', 'scalelog'
            LBlog[i] = LB;
            UBlog[i] = UB;
            steplog[i] = step;
            scalelog[i] = scale;
        }
        
        // show the total process
        reportLog(LBlog, UBlog, scalelog, steplog, iterlimit);
        
        delete [] supply;
        delete [] demand;
        
        for (size_t i = 0; i < ncites; ++i)
            delete [] shipcost[i];
        delete [] shipcost;
        
        for (size_t i = 0; i < ncites; ++i)
            delete [] ship[i];
        delete [] ship;
        
        delete [] build;
        delete [] lambda;
        delete [] slack;
        delete [] steplog;
        delete [] scalelog;
        delete [] LBlog;
        delete [] UBlog;
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

double relaxUB(GRBModel &model) {
    double LB = 0.0;
    size_t numvars = model.get(GRB_IntAttr_NumVars);
    char *vtype = new char [numvars];
    GRBVar *vars = model.getVars();
    
    for (size_t i = 0; i < numvars; ++i) {
        vtype[i] = vars[i].get(GRB_CharAttr_VType);
        vars[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }
    
    model.optimize();
    
    LB = model.get(GRB_DoubleAttr_ObjVal);
    
    for (size_t i = 0; i < numvars; ++i)
        vars[i].set(GRB_CharAttr_VType, vtype[i]);
    
    delete [] vtype;
    
    return LB;
}

double calculateNorm(double *slack, size_t ncites) {
    double norm = 0.0;
    
    for (size_t i = 0; i < ncites; ++i)
        norm += pow(slack[i], 2.0);
    
    return norm;
}

void reportLog(double *LBlog, double *UBlog, double *scalelog, double *steplog, size_t count) {
    printf("\n                *** Summary Report ***               \n");
    printf("  Iter        LB              UB          scale        step\n");
    
    for (size_t i = 0; i < count; ++i)
        printf(" %3d    %12.6f    %12.6f    %8.6f    %8.6f\n", i, LBlog[i], UBlog[i], scalelog[i], steplog[i]);
}