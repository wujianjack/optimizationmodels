param n integer;
param niter integer;
param dtime;

param cR = 7.0;
param cr = 1.0;

set dim = 1..n;

var x {dim};
var y {dim};
var r >= 0;

maximize robj: r;

subject to distance{i in dim, j in dim: j > i}:
	(x[i] - x[j])^2 + (y[i] - y[j])^2 >= 4*r^2;

subject to each{i in dim}:
	x[i]^2 + y[i]^2 <= (cR - r)^2;

subject to xbound_rhs{i in dim}:
	x[i] <= cR - r;

subject to xbound_lhs{i in dim}:
    x[i] >= r - cR;

subject to ybound_rhs{i in dim}:
    y[i] <= cR - r;

subject to ybound_lhs{i in dim}:
    y[i] >= r - cR;

let n := floor(cR^2 / cr);

option log_file "circle_loop.log";
option solver_msg 0;
option reset_initial_guesses 1; # unset warm start to avoid poor solution

# option solver knitro;
# option knitro_options 'outlev 0';
option solver ipopt;
option ipopt_options 'outlev 0 linear_solver pardiso';

let niter := 1;
let dtime := 0.0;

repeat {
    solve;
    
    if robj - cr < 1e-6 then {
        let n := n - 1;
        let niter := niter + 1;
        let dtime := dtime + _solve_elapsed_time;
    }
    else
        break;
    
    if niter >= 100 then {
        break;
    }
}

printf "\n\n       ********Solution Report********       \n";
printf "\nObjective:\n";
printf "  %.12f\n", robj;
printf "\nVariables:\n";
printf "  Dim             X                Y\n";
for {i in dim} {
    printf "  %-2d      %15.12f  %15.12f\n", i, x[i], y[i];
}
printf "\nIterations: %d\n", niter;
printf "\nElapsed time: %.3f (seconds)\n", dtime;

exit;