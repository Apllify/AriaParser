# We fix the dimensions to 3
param dim = 3;
set DIM := 1..dim;

#Scale factor
param scale;

set ATOMS; #set of given atoms
set COVDISTS within {ATOMS, ATOMS}; #index set distances given from covalent bonds
set ANGDISTS within {ATOMS, ATOMS}; #index set distances given from angles
set QDISTS within {ATOMS, ATOMS}; #index set over distances given from Qs
set ALLDISTS := QDISTS union ANGDISTS union COVDISTS;

param numRHOS;
set RHOS := 1..numRHOS;

set NOE_A1{RHOS}; #first set of atoms related to specific rho
set NOE_A2{RHOS}; #second set of atoms realted to specific rho
#set of NOE distance edges
set F{r in RHOS} := {a1 in NOE_A1[r], a2 in NOE_A2[r]: a1 != a2};
set FA := union{r in RHOS} F[r];

<<<<<<< Updated upstream
param CovDists{COVDISTS}; 

#angle dists need more info to properly handle slack
param AngDists{ANGDISTS};
param d1d2{ANGDISTS};
param SinAlpha{ANGDISTS};
=======
set NON_DIST := {ATOMS, ATOMS} diff ALLDISTS diff FA;
>>>>>>> Stashed changes

param AllDists{ALLDISTS};
param rho{RHOS}; #rhos/peak volume from NOE-experiement

## big-M for bounding variables
param maxDist := max{(u,v) in ALLDISTS} AllDists[u,v];
param M := card(ATOMS) * maxDist / 2;

param init_x{ATOMS, DIM};
var x{ATOMS, DIM} >= -M, <= M;

var sm{RHOS} >=0, default 0; #slack variable
var sp{RHOS} >= 0, default 0; #slack variable

#squared inter-hydrogen distances to be determined from NOE-data
#lower bound: vdW radius, upper bound: property of NOE experiment
var d2{FA} >= 1, <= 36, default Uniform(1, 5); 

## slack variables on equations for angle bond distances
<<<<<<< Updated upstream
param AngSCoeff := 10^4; 
var angSlack{ANGDISTS} default 0;

minimize ERROR: 
        sum{r in RHOS} (sm[r] + sp[r])
        + AngSCoeff*sum{(u,v)in ANGDISTS} angSlack[u, v]^2;
=======
param sdistCoeff := 0.05; # obj fun coeff for covalent error term

var sdist{ALLDISTS} default 0;

minimize ERROR: 
        sum{r in RHOS} (sm[r] + sp[r])
        + sdistCoeff*sum{(u,v)in ALLDISTS} sdist[u, v]^2
        #- 0.00000000001*sum{(u,v) in NON_DIST} sum{k in DIM} (x[u,k] - x[v,k])^2
        ;
>>>>>>> Stashed changes

#relate NOE distances abd rhos
subject to NOE_lb { r in RHOS}: 
        rho[r] - sm[r] <= sum{(u, v) in F[r]} scale / d2[u,v]^3;
subject to NOE_ub {r in RHOS}: 
        sum{(u, v) in F[r]} scale / d2[u,v]^3 <= rho[r] + sp[r];

subject to NOE_dists{(u,v) in FA}:
        sum{k in DIM} (x[u,k] - x[v,k])^2 = d2[u,v];
<<<<<<< Updated upstream
subject to cov_dists {(u,v) in COVDISTS}:
        sum{k in DIM} (x[u,k] - x[v,k])^2 = CovDists[u,v]^2;

#only uncomment if using taylor exp or trig function for angSlack
subject to ANG_SLACK{(u, v) in ANGDISTS} : 
        -4 * atan(1)<=angSlack[u,v]<=4 * atan(1);


### NORMAL NAIVE LINEAR ANGLE SLACK VERSION
# subject to ANG_DISTS {(u,v) in ANGDISTS}:
#         sum{k in DIM} (x[u,k] - x[v,k])^2 = AngDists[u,v]^2 + angSlack[u, v];

### TAYLOR EXPANSION VERSION
subject to ANG_DISTS {(u,v) in ANGDISTS}:
        sum{k in DIM} (x[u,k] - x[v,k])^2 = AngDists[u,v]^2 + 
                                            d1d2[u, v] * (SinAlpha[u, v] * (angSlack[u, v] - angSlack[u, v]^3/6) + 
                                                          angSlack[u, v]^2/2 - angSlack[u,v]^4/24); 

### FULL TRIGONOMETRIC VERSION
# subject to ANG_DISTS {(u,v) in ANGDISTS}:
#         sum{k in DIM} (x[u,k] - x[v,k])^2 = AngDists[u,v]^2 + 
#                                             d1d2[u, v] * (SinAlpha[u, v] * sin(angSlack[u, v]) - cos(angSlack[u, v]) + 1);
=======

#distance slacks
subject to dists_constraint {(u,v) in ALLDISTS}:
        sum{k in DIM} (x[u,k] - x[v,k])^2 <= AllDists[u,v]^2+sdist[u,v];
>>>>>>> Stashed changes

# zero centroid
subject to Centroid {k in DIM}: 
        sum{a in ATOMS} x[a, k] = 0;