#include "gurobi_c++.h"

#include <iostream>

using namespace std;

#define MAX_CGTIME 1000

void reportRMP(GRBModel &model);
void reportSUB(GRBModel &model);
void reportMIP(GRBModel &model);

int main(int argc, char *argv[]) {
    try {
        // Input data
        int rollwidth = 115;
        int size[5] = {25, 40, 50, 55, 70};
        int amount[5] = {50, 36, 24, 8, 30};
        int nwidth = 5;
        // End data
        
        GRBEnv env = GRBEnv();
        GRBModel rmp = GRBModel(env);
        GRBModel sub = GRBModel(env);
        
        GRBConstr *rmp_con = NULL;
        GRBVar *rmp_var = NULL;
        GRBVar *sub_var = new GRBVar[nwidth];
        
        GRBLinExpr lexpr = 0;
        
        double *rmp_pi = new double [nwidth];
        double *rmp_coeff = new double [nwidth];
        
        rmp.set(GRB_IntParam_OutputFlag, 0);
        sub.set(GRB_IntParam_OutputFlag, 0);
        
        // Construct RMP 
        rmp_con = rmp.addConstrs(nwidth);
        
        for (int i = 0; i < nwidth; ++i) {
            rmp_con[i].set(GRB_CharAttr_Sense, GRB_GREATER_EQUAL);
            rmp_con[i].set(GRB_DoubleAttr_RHS, amount[i]);
        }
        
        for (int i = 0; i < nwidth; ++i)
            rmp_coeff[i] = 0.0;
        
        for (int i = 0; i < nwidth; ++i) {
            rmp_coeff[i] = (int) (rollwidth / size[i]);
            rmp.addVar(0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, nwidth, rmp_con, rmp_coeff);
            rmp_coeff[i] = 0.0;
        }
        
        rmp.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        // End RMP
        
        // Construct SUB
        for (int i = 0; i < nwidth; ++i)
            sub_var[i] = sub.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER);
        
        for (int i = 0; i < nwidth; ++i)
            lexpr += size[i] * sub_var[i];
        
        sub.addConstr(lexpr, GRB_LESS_EQUAL, rollwidth);
        // End SUB
        
        cout << "               *** Main Loop ***               " << endl;
        for (int niter = 0; niter < MAX_CGTIME; ++niter) {
            cout << "Iteration: " << niter + 1 << endl;
            
            rmp.optimize();
            reportRMP(rmp);
            
            for (int i = 0; i < nwidth; ++i)
                rmp_pi[i] = rmp_con[i].get(GRB_DoubleAttr_Pi);
            
            lexpr = 1;
            for (int i = 0; i < nwidth; ++i)
                lexpr += -rmp_pi[i] * sub_var[i];
            
            sub.setObjective(lexpr, GRB_MINIMIZE);
            sub.optimize();
            reportSUB(sub);
            
            if (sub.get(GRB_DoubleAttr_ObjVal) > -1e-6)
                break;
                
            for (int i = 0; i < nwidth; ++i)
                rmp_coeff[i] = sub_var[i].get(GRB_DoubleAttr_X);
                
            rmp.addVar(0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, nwidth, rmp_con, rmp_coeff);
        }
        cout << "               *** End Loop ***               " << endl;
        
        rmp_var = rmp.getVars();
        
        for (int i = 0; i < rmp.get(GRB_IntAttr_NumVars); ++i)
            rmp_var[i].set(GRB_CharAttr_VType, GRB_INTEGER);
        
        rmp.optimize();
        reportMIP(rmp);
        
        delete [] rmp_con;
        delete [] sub_var;
        delete [] rmp_pi;
        delete [] rmp_coeff;
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

void reportRMP(GRBModel &model) {
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        cout << "Using " << model.get(GRB_DoubleAttr_ObjVal) << " rolls" << endl;
        cout << endl;
        
        GRBVar *var = model.getVars();
        for (int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i)
            cout << var[i].get(GRB_StringAttr_VarName) << " = " << var[i].get(GRB_DoubleAttr_X) << endl;
        cout << endl;
        
        GRBConstr *con = model.getConstrs();
        for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i)
            cout << con[i].get(GRB_StringAttr_ConstrName) << " = " << con[i].get(GRB_DoubleAttr_Pi) << endl;
        cout << endl;
    }
}

void reportSUB(GRBModel &model) {
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        cout << "Pi: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
        cout << endl;
        
        if (model.get(GRB_DoubleAttr_ObjVal) <= -1e-6) {
            GRBVar *var = model.getVars();
            
            for (int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i)
                cout << var[i].get(GRB_StringAttr_VarName) << " = " << var[i].get(GRB_DoubleAttr_X) << endl;
            cout << endl;
        }
    }
}

void reportMIP(GRBModel &model) {
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        cout << endl;
        cout << "Best MIP Solution: " << model.get(GRB_DoubleAttr_ObjVal) << " rolls" << endl;
        cout << endl;
        
        GRBVar *var = model.getVars();
        for (int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i)
            cout << var[i].get(GRB_StringAttr_VarName) << " = " << var[i].get(GRB_DoubleAttr_X) << endl;
        cout << endl;
    }
}