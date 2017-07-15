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

  {  // Problem Setting Block
     /* Number of objectives */
    int nobjs = 1;
    int *objsense = (int *) malloc(nobjs * sizeof(int));
    int *objs_beg = (int *) malloc(nobjs * sizeof(int));
    int *objs_length = (int *) malloc(nobjs * sizeof(int));

    /* Number of constraints */
    int ncons = 1;
    char *ctype = (char *) malloc(ncons * sizeof(char));
    int *cons_beg = (int *) malloc(ncons * sizeof(int));
    int *cons_length = (int *) malloc(ncons * sizeof(int));

    /* Number of variables */
    int nvars = 3;
    char *vtype = (char *) malloc(nvars * sizeof(char));
    double *varval = (double *) malloc(nvars * sizeof(double));
    double *lwrbnd = (double *) malloc(nvars * sizeof(double));
    double *uprbnd = (double *) malloc(nvars * sizeof(double));

    /* Number of real number constants */
    int nnums = 3;
    double *numval = (double *) malloc(nnums * sizeof(double));

    /* Number of items in the instruction lists */
    int lsize = 29;
    int *code = (int *) malloc(lsize * sizeof(int));

    /* Count for instruction code */
    int ikod = 0;
    /* Count for objective row */
    int iobj = 0;
    /* Count for constraint row */
    int icon = 0;

    int i;
    int nLinearz, nAutoDeriv, nConvexRelax, nCRAlgReform;

    /* Numerical constants in model */
    numval[0] = 1.0;
    numval[1] = 50.0;
    numval[2] = exp(1.0);

    /* Lower & upper bounds of variables */
    for (i = 0; i < nvars; ++i) {
      if (i == nvars - 1) {
        lwrbnd[i] = numval[2];
        uprbnd[i] = 1e20;
      } else {
        lwrbnd[i] = 0.0;
        uprbnd[i] = 100.0;
      }
    }

    /* Starting point of variables */
    for (i = 0; i < nvars; ++i) varval[i] = 1.234567;

    /* Variable type, C = continuous, B = binary */
    for (i = 0; i < nvars; ++i) vtype[i] = 'C';

    /*
     *  Instruction code of the objective:
     *  max  sin(x[2]) / x[2] + 1;
     */

    /* Direction of optimization */
    objsense[iobj] = LS_MAX;
    /* Beginning position of objective */
    objs_beg[iobj] = ikod;
    /* Instruction list code */
    code[ikod++] = EP_PUSH_VAR;
    code[ikod++] = 2;
    code[ikod++] = EP_SIN;
    code[ikod++] = EP_PUSH_VAR;
    code[ikod++] = 2;
    code[ikod++] = EP_DIVIDE;
    code[ikod++] = EP_PUSH_NUM;
    code[ikod++] = 0;
    code[ikod++] = EP_PLUS;

    /* Length of objective */
    objs_length[iobj] = ikod - objs_beg[iobj];
    /* Increment the objective count */
    iobj++;

    /*
     *  Instruction code of constraint 0:
     *   x[2] - sqrt(sqr(x[0] - 50.0) + sqr(x[1] - 50)) - exp(1) = 0;
     */

    /* Constraint type */
    ctype[icon] = 'E'; /* less or than or equal to */
    /* Beginning position of constraint 0 */
    cons_beg[icon] = ikod;
    /* Instruction list code */
    code[ikod++] = EP_PUSH_VAR;
    code[ikod++] = 2;
    code[ikod++] = EP_PUSH_VAR;
    code[ikod++] = 0;
    code[ikod++] = EP_PUSH_NUM;
    code[ikod++] = 1;
    code[ikod++] = EP_MINUS;
    code[ikod++] = EP_SQR;
    code[ikod++] = EP_PUSH_VAR;
    code[ikod++] = 1;
    code[ikod++] = EP_PUSH_NUM;
    code[ikod++] = 1;
    code[ikod++] = EP_MINUS;
    code[ikod++] = EP_SQR;
    code[ikod++] = EP_PLUS;
    code[ikod++] = EP_SQRT;
    code[ikod++] = EP_MINUS;
    code[ikod++] = EP_PUSH_NUM;
    code[ikod++] = 2;
    code[ikod++] = EP_MINUS;

    /* Length of constraint 0 */
    cons_length[icon] = ikod - cons_beg[icon];
    /* Increment the constraint count */
    icon++;

    /* Total number of items in the instruction list */
    lsize = ikod;

    /* Set linearization level, before a call to LSloadNLPCode.
     * If not specified, the solver will decide */
    nLinearz = 1;
    nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_NLP_LINEARZ, nLinearz);
    APIERRORCHECK;

    /* Select algebraic reformulation level in convex relaxation*/
    nCRAlgReform = 1;
    nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_NLP_CR_ALG_REFORM, nCRAlgReform);
    APIERRORCHECK;

    /* Select convex relax level */
    nConvexRelax = 0;
    nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_NLP_CONVEXRELAX, nConvexRelax);
    APIERRORCHECK;

    /* Set up automatic differentiation, before a call to LSloadNLPCode.
     * If not specified, the numerical derivative will be applied */
    nAutoDeriv = 1;
    nErrorCode = LSsetModelIntParameter(pModel, LS_IPARAM_NLP_AUTODERIV, nAutoDeriv);
    APIERRORCHECK;

    /* Pass the instruction list to problem structure by a call to
     * LSloadNLPCode() */
    nErrorCode = LSloadInstruct(pModel, ncons, nobjs, nvars, nnums, objsense, ctype,
                                vtype, code, lsize, NULL, numval, varval, objs_beg,
                                objs_length, cons_beg, cons_length, lwrbnd, uprbnd);
    APIERRORCHECK;

    free(objsense);
    free(objs_beg);
    free(objs_length);

    free(ctype);
    free(cons_beg);
    free(cons_length);

    free(vtype);
    free(varval);
    free(lwrbnd);
    free(uprbnd);

    free(numval);

    free(code);
  }

  if (argc == 2) {
    if (strcmp(argv[1], "-gop") == 0) {
      printf("Solving with 'Global Optimizer'...\n");
      nErrorCode = LSsolveGOP(pModel, NULL);
    } else {
      printf("usage: peak.exe -gop\n");
      exit(1);
    }
  } else {
    printf("Solving with 'Local Nonlinear Optimizer'...\n");
    nErrorCode = LSoptimize(pModel, LS_METHOD_FREE, NULL);
  }
  APIERRORCHECK;

  {  // Solution Report Block
    int i, status, nvars;
    double pobjval;
    double *primal;
    char *cStatus;

    /* Report the status of solution */
    nErrorCode = LSgetInfo(pModel, LS_IINFO_MODEL_STATUS, &status);
    APIERRORCHECK;

    /* Get the optimization result */
    nErrorCode = LSgetInfo(pModel, LS_IINFO_NUM_VARS, &nvars);
    APIERRORCHECK;

    nErrorCode = LSgetInfo(pModel, LS_DINFO_POBJ, &pobjval);
    APIERRORCHECK;

    primal = (double *) malloc(nvars * sizeof(double));
    nErrorCode = LSgetPrimalSolution(pModel, primal);
    APIERRORCHECK;

    if (status == LS_STATUS_OPTIMAL || status == LS_STATUS_BASIC_OPTIMAL ||
        status == LS_STATUS_LOCAL_OPTIMAL || status == LS_STATUS_FEASIBLE) {
      printf("\n\t  *** Solution Report *** \n");
      printf("Objective value: \n    %.12f\n", pobjval);

      printf("\nVariables: \n");
      for (i = 0; i < nvars; i++)
        printf("    x[%d] = %.12f\n", i, primal[i]);
    } else if (status == LS_STATUS_INFEASIBLE)
      printf("\n\nNo feasible solution. \n");

    free(primal);
  }

  nErrorCode = LSdeleteModel(&pModel);
  APIERRORCHECK;

  LSdeleteEnv(&pEnv);
  APIERRORCHECK;

  return 0;
}