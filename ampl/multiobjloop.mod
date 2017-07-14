# loop style multi-objective model in AMPL
# ref: http://www.gurobi.com/documentation/7.5/examples/multiobj_cpp_cpp.html
# author: wujianw@stu.xjtu.edu.cn

set Groundset = 0..19;
set Subsets = 0..3;

param Budget;
param ObjCoeff {Subsets, Groundset};

var Elem {Groundset} binary;

maximize MultiObject {i in Subsets}:
    sum {j in Groundset} ObjCoeff[i, j] * Elem[j];

subject to conBudget:
    sum {j in Groundset} Elem[j] <= Budget;

# assign parameter
data;
param Budget := 12;
param ObjCoeff (tr):
          0    1    2    3 :=
      0   1    0    0    0
      1   1    0    0    0
      2   1    0    0    0
      3   1    0    1    1
      4   1    0    1    1
      5   1    1    0    1
      6   1    1    1    0
      7   1    1    1    0
      8   1    1    0    0
      9   1    1    0    1
      10  0    0    0    1
      11  0    0    0    1
      12  0    0    0    0
      13  0    0    1    0
      14  0    0    1    0
      15  0    1    0    1
      16  0    1    1    1
      17  0    1    1    1
      18  0    1    0    0
      19  0    1    0    0;

suffix objpriority IN;
suffix objweight IN;
suffix objreltol IN;
suffix objabstol IN;

# assign objective priority
let MultiObject[0].objpriority := 3;
let MultiObject[1].objpriority := 2;
let MultiObject[2].objpriority := 2;
let MultiObject[3].objpriority := 1;

# assign objective weight
let MultiObject[0].objweight := 1.0;
let MultiObject[1].objweight := 0.25;
let MultiObject[2].objweight := 1.25;
let MultiObject[3].objweight := 1.0;

# assign relative tolerance
let {i in Subsets} MultiObject[i].objreltol := 0.01;

# assign absolute tolerance
let {i in Subsets} MultiObject[i].objabstol := 1.0 + i;

# assign options to solver
option solver gurobi;
option gurobi_options 'poolsolutions 100 multiobj 1';

# solve model
solve;

# display solution
display MultiObject;
display Elem;

# exit program
exit;