# We fix the dimensions to 3
param dim = 3;
set DIM := 1..dim;

set ATOMS; #set of given atoms
set COVDISTS within {ATOMS, ATOMS}; #index set distances given from covalent bonds
set ANGDISTS within {ATOMS, ATOMS}; #index set distances given from angles

param numRHOS;
set RHOS := 1..numRHOS;

set NOE_A1{RHOS}; #first set of atoms related to specific rho
set NOE_A2{RHOS}; #second set of atoms realted to specific rho
#set of NOE distance edges
set F{r in RHOS} := {a1 in NOE_A1[r], a2 in NOE_A2[r]};
set FA := union{r in RHOS} F[r];

param CovDists{COVDISTS}; 
param AngDists{ANGDISTS};

param rho{RHOS}; #rhos/peak volume from NOE-experiement

## big-M for bounding variables
param maxCovDist := max{(u,v) in COVDISTS} CovDists[u,v];
param M := card(ATOMS) * maxCovDist / 2;

var x{ATOMS, DIM} >= -M, <= M, default Uniform(-M/10, M/10);

var sm{RHOS} >=0, default 0; #slack variable
var sp{RHOS} >= 0, default 0; #slack variable

#squared inter-hydrogen distances to be determined from NOE-data
#lower bound: vdW radius, upper bound: property of NOE experiment
var d2{FA} >= 1, <= 36, default Uniform(1, 5); 

## slack variables on equations for angle bond distances
param scovCoeff := 0.05; # obj fun coeff for covalent error term
var scov{ANGDISTS} default 0;

minimize ERROR: 
        sum{r in RHOS} (sm[r] + sp[r])
        + scovCoeff*sum{(u,v)in ANGDISTS} scov[u, v]^2;

#relate NOE distances abd rhos
subject to NOE_lb { r in RHOS}: 
        rho[r] - sm[r] <= sum{(u, v) in F[r]}1 / d2[u,v]^3;
subject to NOE_ub {r in RHOS}: 
        sum{(u, v) in F[r]} 1 / d2[u,v]^3 <= rho[r] + sp[r];

#
subject to NOE_dists{(u,v) in FA}:
        sum{k in DIM} (x[u,k] - x[v,k])^2 = d2[u,v];
subject to cov_dists {(u,v) in COVDISTS}:
        sum{k in DIM} (x[u,k] - x[v,k])^2 = CovDists[u,v]^2;
subject to ang_dists {(u,v) in ANGDISTS}:
        sum{k in DIM} (x[u,k] - x[v,k])^2 = AngDists[u,v]^2 + scov[u,v];

# zero centroid
subject to Centroid {k in DIM}: 
        sum{a in ATOMS} x[a, k] = 0;