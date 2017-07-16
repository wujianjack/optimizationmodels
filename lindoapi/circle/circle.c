#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lindo.h"

#define APIERRORCHECK                                            \
  if (nErrorCode) {                                              \
    if (pEnv) {                                                  \
      LSgetErrorMessage(pEnv, nErrorCode, cErrorMessage);        \
      printf("Errorcode = %d: %s\n", nErrorCode, cErrorMessage); \
    } else                                                       \
      printf("Fatal Error\n");                                   \
                                                                 \
    exit(1);                                                     \
  }

int ntotal = 5;
double rtotal = 5.0;

static void LS_CALLTYPE getLogInfo(pLSmodel pModel, char *line, void *userdata) {
  if (line) printf("%s", line);
}

int loadModel(pLSmodel pModel) {
  int nObjs = 1;
  int *panObjSense = (int *) malloc(nObjs * sizeof(int));
  int *paiObj = (int *) malloc(nObjs * sizeof(int));
  int *panObj = (int *) malloc(nObjs * sizeof(int));

  int nCons = ntotal * (ntotal - 1) / 2 + 5 * ntotal;
  char *pacConType = (char *) malloc(nCons * sizeof(char));
  int *paiRows = (int *) malloc(nCons * sizeof(int));
  int *panRows = (int *) malloc(nCons * sizeof(int));

  int nVars = 2 * ntotal + 1;
  char *pacVarType = (char *) malloc(nVars * sizeof(char));
  int *paiVars = NULL;
  double *padX0 = (double *) malloc(nVars * sizeof(double));
  double *padL = (double *) malloc(nVars * sizeof(double));
  double *padU = (double *) malloc(nVars * sizeof(double));

  int nNums = 3;
  double *padVals = (double *) malloc(nNums * sizeof(double));

  int nCode = 10000;
  int *panCode = (int *) malloc(nCode * sizeof(int));

  int i, j, nErrorCode;
  int ikod = 0, iobj = 0, icon = 0;

  for (i = 0; i < nVars; i++) {
    if (i == nVars - 1) {
      padL[i] = 0.0;
      padU[i] = rtotal;
    } else {
      padL[i] = -rtotal;
      padU[i] = rtotal;
    }

    padX0[i] = 0.00000;
    pacVarType[i] = 'C';
  }

  //
  padVals[0] = 2.0;
  padVals[1] = 4.0;
  padVals[2] = rtotal;
  //

  /*********************
   variable name vs index
   * 0 ~ ntotal-1                 X
   * ntotal ~ 2*ntotal - 1        Y
   * 2*ntotal                     R
   **********************/

  // MAX  R
  panObjSense[iobj] = LS_MAX;
  paiObj[iobj] = ikod;
  panCode[ikod++] = EP_PUSH_VAR;
  panCode[ikod++] = nVars - 1;
  panObj[iobj] = ikod - paiObj[iobj];
  iobj++;

  // @FOR(DIM(I): @FOR(DIM(J)| J #GT# I: (X(I) - X(J))^2 + (Y(I) - Y(J))^2 >= 4*R^2));
  for (i = 0; i < ntotal; i++) {
    for (j = 0; j < ntotal; j++) {
      if (j > i) {
        pacConType[icon] = 'G';
        paiRows[icon] = ikod;
        panCode[ikod++] = EP_PUSH_VAR;
        panCode[ikod++] = i;
        panCode[ikod++] = EP_PUSH_VAR;
        panCode[ikod++] = j;
        panCode[ikod++] = EP_MINUS;
        panCode[ikod++] = EP_PUSH_NUM;
        panCode[ikod++] = 0;
        panCode[ikod++] = EP_POWER;
        panCode[ikod++] = EP_PUSH_VAR;
        panCode[ikod++] = i + ntotal;
        panCode[ikod++] = EP_PUSH_VAR;
        panCode[ikod++] = j + ntotal;
        panCode[ikod++] = EP_MINUS;
        panCode[ikod++] = EP_PUSH_NUM;
        panCode[ikod++] = 0;
        panCode[ikod++] = EP_POWER;
        panCode[ikod++] = EP_PLUS;
        panCode[ikod++] = EP_PUSH_NUM;
        panCode[ikod++] = 1;
        panCode[ikod++] = EP_PUSH_VAR;
        panCode[ikod++] = nVars - 1;
        panCode[ikod++] = EP_PUSH_NUM;
        panCode[ikod++] = 0;
        panCode[ikod++] = EP_POWER;
        panCode[ikod++] = EP_MULTIPLY;
        panCode[ikod++] = EP_MINUS;
        panRows[icon] = ikod - paiRows[icon];
        icon++;
      }
    }
  }

  // @FOR(DIM(I): X(I)^2 + Y(I)^2 <= (C - R)^2);
  for (i = 0; i < ntotal; i++) {
    pacConType[icon] = 'L';
    paiRows[icon] = ikod;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = i;
    panCode[ikod++] = EP_PUSH_NUM;
    panCode[ikod++] = 0;
    panCode[ikod++] = EP_POWER;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = i + ntotal;
    panCode[ikod++] = EP_PUSH_NUM;
    panCode[ikod++] = 0;
    panCode[ikod++] = EP_POWER;
    panCode[ikod++] = EP_PLUS;
    panCode[ikod++] = EP_PUSH_NUM;
    panCode[ikod++] = 2;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = nVars - 1;
    panCode[ikod++] = EP_MINUS;
    panCode[ikod++] = EP_PUSH_NUM;
    panCode[ikod++] = 0;
    panCode[ikod++] = EP_POWER;
    panCode[ikod++] = EP_MINUS;
    panRows[icon] = ikod - paiRows[icon];
    icon++;
  }

  // @FOR(DIM(I): X(I) >= R - C; X(I) <= C - R);
  for (i = 0; i < ntotal; i++) {
    pacConType[icon] = 'G';
    paiRows[icon] = ikod;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = i;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = nVars - 1;
    panCode[ikod++] = EP_MINUS;
    panCode[ikod++] = EP_PUSH_NUM;
    panCode[ikod++] = 1;
    panCode[ikod++] = EP_PLUS;
    panRows[icon] = ikod - paiRows[icon];
    icon++;
  }

  for (i = 0; i < ntotal; i++) {
    pacConType[icon] = 'L';
    paiRows[icon] = ikod;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = i;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = nVars - 1;
    panCode[ikod++] = EP_PLUS;
    panCode[ikod++] = EP_PUSH_NUM;
    panCode[ikod++] = 1;
    panCode[ikod++] = EP_MINUS;
    panRows[icon] = ikod - paiRows[icon];
    icon++;
  }

  // @FOR(DIM(I): Y(I) >= R - C; Y(I) <= C - R);
  for (i = 0; i < ntotal; i++) {
    pacConType[icon] = 'G';
    paiRows[icon] = ikod;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = i + ntotal;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = nVars - 1;
    panCode[ikod++] = EP_MINUS;
    panCode[ikod++] = EP_PUSH_NUM;
    panCode[ikod++] = 1;
    panCode[ikod++] = EP_PLUS;
    panRows[icon] = ikod - paiRows[icon];
    icon++;
  }

  for (i = 0; i < ntotal; i++) {
    pacConType[icon] = 'L';
    paiRows[icon] = ikod;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = i + ntotal;
    panCode[ikod++] = EP_PUSH_VAR;
    panCode[ikod++] = nVars - 1;
    panCode[ikod++] = EP_PLUS;
    panCode[ikod++] = EP_PUSH_NUM;
    panCode[ikod++] = 1;
    panCode[ikod++] = EP_MINUS;
    panRows[icon] = ikod - paiRows[icon];
    icon++;
  }

  nCode = ikod;

  nErrorCode = LSloadInstruct(pModel, nCons, nObjs, nVars, nNums, panObjSense,
                              pacConType, pacVarType, panCode, nCode, paiVars, padVals,
                              padX0, paiObj, panObj, paiRows, panRows, padL, padU);

  free(panObjSense);
  free(paiObj);
  free(panObj);

  free(pacConType);
  free(paiRows);
  free(panRows);

  free(pacVarType);
  // free(paiVars);
  free(padX0);
  free(padL);
  free(padU);

  free(padVals);

  free(panCode);

  return nErrorCode;
}

int main(int argc, char *argv[]) {
  int nErrorCode;
  char cErrorMessage[LS_MAX_ERROR_MESSAGE_LENGTH];

  pLSenv pEnv;
  pLSmodel pModel;

  pEnv = LScreateEnv(&nErrorCode, NULL);
  if (nErrorCode == LSERR_NO_VALID_LICENSE) {
    printf("Invalid License Key!\n");
    exit(1);
  }
  APIERRORCHECK;

  pModel = LScreateModel(pEnv, &nErrorCode);
  APIERRORCHECK;

  nErrorCode = loadModel(pModel);
  APIERRORCHECK;

  nErrorCode = LSsetModelLogfunc(pModel, (printModelLOG_t)getLogInfo, NULL);
  APIERRORCHECK;
  nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_NLP_AUTODERIV, 1);
  APIERRORCHECK;
  nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_NLP_PRINTLEVEL, 1);
  APIERRORCHECK;
  nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_GOP_PRINTLEVEL, 1);
  APIERRORCHECK;
  // nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_GOP_TIMLIM, 100);
  // APIERRORCHECK;

  nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_NLP_SOLVER, LS_NMETHOD_MSW_GRG);
  APIERRORCHECK;
  nErrorCode = LSoptimize(pModel, LS_METHOD_FREE, NULL);
  APIERRORCHECK;

  // nErrorCode = LSsolveGOP(pModel, NULL);
  // APIERRORCHECK;

  printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  {
    int i, nVars, nCons, nStatus;
    double pObjVal;
    double *primalSol;

    nErrorCode = LSgetInfo(pModel, LS_IINFO_MODEL_STATUS, &nStatus);
    APIERRORCHECK;

    if (nStatus == LS_STATUS_OPTIMAL || nStatus == LS_STATUS_BASIC_OPTIMAL ||
        nStatus == LS_STATUS_LOCAL_OPTIMAL || nStatus == LS_STATUS_FEASIBLE) {
      nErrorCode = LSgetInfo(pModel, LS_IINFO_NUM_VARS, &nVars);
      APIERRORCHECK;
      nErrorCode = LSgetInfo(pModel, LS_IINFO_NUM_CONS, &nCons);
      APIERRORCHECK;
      nErrorCode = LSgetInfo(pModel, LS_DINFO_POBJ, &pObjVal);
      APIERRORCHECK;

      primalSol = (double *)malloc(nVars * sizeof(double));
      nErrorCode = LSgetPrimalSolution(pModel, primalSol);
      APIERRORCHECK;

      printf("\n\t  *** Solution Report *** \n");
      printf("Objective value: \n    %.12f\n", pObjVal);

      printf("\nVariables: \n");
      for (i = 0; i < nVars; i++)
        printf("    x[%d] = %.12f\n", i, primalSol[i]);

      free(primalSol);
    } else if (nStatus == LS_STATUS_INFEASIBLE)
      printf("\n\nNo feasible solution.\n");
  }

  nErrorCode = LSdeleteModel(&pModel);
  APIERRORCHECK;

  nErrorCode = LSdeleteEnv(&pEnv);
  APIERRORCHECK;

  return 0;
}