# We fix the dimensions to 3

set DIM := 1..3;
set ATOMS; #set of given atoms
set COVDISTS within {ATOMS, ATOMS}; #index set distances given from covalent bonds + angles
set RHOS within {ATOMS, ATOMS, ATOMS}; #index set over rho
set DISTS within {ATOMS, ATOMS};

param CovDists{COVDISTS}; 
param rho{RHOS}; #rhos/peak volume from NOE-experiement

var sm{RHOS} >=0; #slack variable
var sp{RHOS} >= 0; #slack variable
var x{ATOMS, DIM};
var d2{DISTS} >= 0; #squared inter-hydrogen distances to be determined from NOE-data

minimize ERROR: 
    sum{(a1, a2) in COVDISTS} ((sum{k in DIM} (x[a1, k] - x[a2, k]))^2 - CovDists[a1, a2]^2)^2
    #correct assignment for distances given covalent bonds + angles
    + sum{(a1, a2) in DISTS} ((sum{k in DIM}(x[a1, k] - x[a2, k]))^2 - d2[a1, a2])^2
    # + sum of slack varibles squared??
    + sum{(p, q, r) in RHOS} (sm[p, q, r] + sp[p, q, r]);
    #correct assignment from NOE-data

subject to NOE_lb {(p, q, r) in RHOS}: 
        rho[p,q,r] - sm[p,q,r] <= 1 / d2[p,q]^3;
 subject to NOE_ub {(p, q, r) in RHOS}: 
        1 / d2[p,q]^3 <= rho[p,q,r] + sp[p,q,r];