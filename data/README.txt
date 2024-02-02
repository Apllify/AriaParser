Output from NOE experiments based on chemical shifts

*.peaks: peak information (chemical shift dimensions, intensity, error)
  first 2-4 adjacent numerical columns:
    chemical shifts from atom1, atom2, atom3, atom4
  next two numerical columns: intensity, error
    ignore the error, use slack variables
    ignore negative intensities
  
*.prot: proton information (chemical shift, error, atom name)
  first numerical column: chemical shift
  second numerical column: error (ignore)
  next textual column: atom type

The two files are correlated by the value of chemical shifts:
  Open a .peaks file. Consider the line
  Consider the line
    2 29.646 2.283 4.318 1 U 8.214e+04 0.00e+00 m 0 150 151 149 0
  coming from a 3-dimensional NOE experiment;
  the chemical shifts refer to
    Carbon (col 1), Hydrogen (col 2), Hydrogen (col 3)
      (this information is not written in the file, perhaps in filename?)
  We now look at chemical shift values in the corresponding .prot file:
     150  29.646 0.000  CB  9
     170  29.646 0.000  CG 10
    1262  29.646 0.000  CB 73
  Since neither 2.283 nor 4.318 appear in the .prot file, we discard the peak.
  Consider the line:
    12  52.572   4.548   4.548 5 U  -2.234e+06  0.00e+00 m  0  219  220  220 0
  We look for chemical shift value 52.572 (col1) in the corresponding .prot:
    we find "219 52.572 0.000 CA 13", meaning:
      carbon-alpha on residue 13 has chemical shift 52.572
  next we look for chemical shift 4.548 (col2 and col3) in the .prot:
    we find "220 4.548 0.000 HA 13", meaning
      hydrogen-alpha on residue 13 has chemical shift 4.548
  => p=52.572, q=4.548, r=4.548, I_p={CA13}, I_q={HA13}, I_r={HA13}
     rho_pqr=-2.234e+06
  Since rho_pqr<0, we discard the peak.

*.top: topology information for every protein
  look at records "atom" and extract names

*.par: bond lengths and angles
  with atom names from .top find "bond" records and read distance in AA
  angles are not precise, as they depend on proteins 

  Therese: 
  Les fichiers .par et .top sont uniques pour toute les proteines, 
    donc, il n'y a pas de correspondance individuelle avec .peaks et .prot.
  La correspondance se fait de la facon suivante: 
    Prenons un residu de la sequence de la proteine, par exemple ALA:
    dans .top, on peut trouver la topologie: residue ALA
    avec les lignes atom qui definissent les noms d'atomes et leur type, 
    ainsi que les lignes bond qui definissent les liaisons covalentes. 
    par exemple: atom N  a comme  type NH1  
  Ensuite, pour avoir les liaisons covalentes, 
    il faut aller dans ibp-protein.par et cherche bond. Par exemple:
      bond CH1E NH1   1.000 1.458
    signifie que la liaison covalente entre les types CH1E et NH1 
    a comme longueur 1.458

------
Therese 240130

Concernant les problemes de format, il a ete cree dans les annees 2000 un data-model python CCPN pour encapsuler toute l'information RMN.

Par ailleurs, ARIA produit des fichiers xml qui contiennent les donnees de structures et les donnees de peaks, et qui sont sans  doute plus propres 
du point de vue de la qualite informatique. 

Mais, ensuite ARIA fait appel a CNS pour calculer les structures par recuit simule. CNS est un programme en Fortran ecrit dans les annees 90, et donc
c'est sans doute formatte de maniere moins stricte. Si on regarde l'ordre des elements dans les lignes bond et pas les numeros de colonne, on voit que
les ordres correspondent:
bond C5   CH2E  1.000 1.497
 BOND  C    CH1E    1000.000 {sd=     0.001}      1.525

En ce qui concerne les angles de liaison, ils sont indispensables pour calculer la structure. Bien sur, on les prend constants et ils sont variables, mais 
on est oblige de faire avec... C'est pour cela que j'essaie de predire les variations par machine learning. Pour le recuit simule, les variations derangent moins bien sur. 

En ce qui concerne les distances et angles entre residus, dans topallhdg5.3.pro, il faut aller voir les patches presidue PEPT et presidue PEPP (pour les residus proline), 
aux lignes 1286 et 1308. Il y a des liaisons et des angles qui sont rajoutes. 

Pour les pseudo-atomes, ils sont toujours la pour les groupes methyl, en effet les trois hydrogenes des groupes methyles resonnent au meme deplacement chimique. 
Par ailleurs, ils peuvent etre la pour des groupes CH2 dont les deux H resonnent a la meme frequence. Comme les CH2 et les CH3 sont des objets rigides, on peut 
s'en tirer en remplacant les H par un pseudo-atome qui est le barycentre des 2 ou 3 H.
