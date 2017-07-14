# non-loop style multi-objective model in AMPL
# ref: http://www.gurobi.com/documentation/7.5/examples/multiobj_cpp_cpp.html
# author: wujianw@stu.xjtu.edu.cn

set Groundset = 0..19;
set Subsets = 0..3;

param Budget;
param ObjCoeff {Subsets, Groundset};

var Elem {Groundset} binary;

maximize Objective_one: sum {j in Groundset} ObjCoeff[0, j] * Elem[j];
maximize Objective_two: sum {j in Groundset} ObjCoeff[1, j] * Elem[j];
maximize Objective_three: sum {j in Groundset} ObjCoeff[2, j] * Elem[j];
maximize Objective_four: sum {j in Groundset} ObjCoeff[3, j] * Elem[j];

subject to conBudget:
    sum {j in Groundset} Elem[j] <= Budget;

# assign parameter
data;
param Budget := 12;
param ObjCoeff :
        0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19 :=
    0   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0
    1   0   0   0   0   0   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1
    2   0   0   0   1   1   0   1   1   0   0   0   0   0   1   1   0   1   1   0   0
    3   0   0   0   1   1   1   0   0   0   1   1   1   0   0   0   1   1   1   0   0;

suffix objpriority IN, integer, >= 0, <= 999;
suffix objweight IN;
suffix objreltol IN;
suffix objabstol IN;

# assign objective priority
let Objective_one.objpriority := 3;
let Objective_two.objpriority := 2;
let Objective_three.objpriority := 2;
let Objective_four.objpriority := 1;

# assign objective weight
let Objective_one.objweight := 1.0;
let Objective_two.objweight := 0.25;
let Objective_three.objweight := 1.25;
let Objective_four.objweight := 1.0;

# assign relative tolerance
let Objective_one.objreltol := 0.01;
let Objective_two.objreltol := 0.01;
let Objective_three.objreltol := 0.01;
let Objective_four.objreltol := 0.01;

# assign absolute tolerance
let Objective_one.objabstol := 1.0;
let Objective_two.objabstol := 2.0;
let Objective_three.objabstol := 3.0;
let Objective_four.objabstol := 4.0;

# assign options to solver
option solver gurobi;
option gurobi_options 'poolsolutions 100 multiobj 1';

# solve model
solve;

# display solution
display Objective_one;
display Objective_two;
display Objective_three;
display Objective_four;
display Elem;

# exit program
exit;