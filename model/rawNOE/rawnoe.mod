##### rawnoe.mod - rawNOE formulation 

### sets and parameters
param Kmax = 3;
set K := 1..Kmax;

## set of chemical shift indices for C or N atoms
set csCN;
# set of C/N atoms for each chemical shift
set CN{csCN};
# set of C/N atoms
set ACN := union{cs in csCN} CN[cs];

## set of chemical shift indices for H atoms (first column)
set csH1;
# set of H atoms for each chemical shift (first column)
set H1{csH1};
# set of H atoms (first column)
set AH1 := union{cs in csH1} H1[cs];

## set of chemical shift indices for H atoms (second column)
set csH2;
# set of H atoms for each chemical shift (second column)
set H2{csH2};
# set of H atoms (second column)
set AH2 := union{cs in csH2} H2[cs];

## set of NOE distance edges
set F{q in csH1, r in csH2} := {h1 in H1[q], h2 in H2[r]};
set FA := union{q in csH1, r in csH2} F[q,r];

## peak intensity at chemical shift triplet
param rho{csCN,csH1,csH2} >= 0, default 0;

## covalent distances
set V;
set A := ACN union AH1 union AH2 union V; # all atoms
set E within {A,A}; 
param covd{E} >= 0;


## big-M for bounding variables
param maxCovDist := max{(u,v) in E} covd[u,v];
param M := card(A) * maxCovDist / 2;

## atom properties (for output purposes)
param atomId{A} integer, >= 0, default 0;
param atomType{A} symbolic, default "none";
param residue{A} integer, >= 0, default 0;

### decision variables

## realization
var x{A,K} >= -M, <= M, default Uniform(-M/10,M/10);

## square distances for NOE values: at least 1AA, at most 6AA
var d2{FA} >= 1, <= 36, default Uniform(1,5);

## slack variables on inequalities
var sm{csCN,csH1,csH2} >= 0, default 0;
var sp{csCN,csH1,csH2} >= 0, default 0;

## slack variables on equations for covalent bond distances
param scovCoeff := 0.05; # obj fun coeff for covalent error term
var scov{E} default 0;

### objective
minimize error :
  # slack for inequalities
  sum{p in csCN, q in csH1, r in csH2} (sm[p,q,r] + sp[p,q,r])
  # slack for realization equations for covalent bonds
  + scovCoeff * sum{(u,v) in E} scov[u,v]^2;

### constraints

## attribution constraints 
subject to attribution1{p in csCN, q in csH1, r in csH2 : rho[p,q,r]>0}:
  rho[p,q,r] - sm[p,q,r] <= sum{(u,v) in F[q,r]} 1 / d2[u,v]^3;
subject to attribution2{p in csCN, q in csH1, r in csH2 : rho[p,q,r]>0}:
  sum{(u,v) in F[q,r]} 1 / d2[u,v]^3 <= rho[p,q,r] + sp[p,q,r];

## realization constraints
subject to realization_noe{(u,v) in FA}:
  sum{k in K} (x[u,k] - x[v,k])^2 = d2[u,v];
subject to realization_cov{(u,v) in E}:
  sum{k in K} (x[u,k] - x[v,k])^2 = covd[u,v]^2 + scov[u,v];

## zero centroid
subject to zeroCentroid{k in K}: sum{v in V} x[v,k] = 0;
