doComputation = 1;
if doComputation
global quadrule_order quadrule_x quadrule_w;
quadrule_order = 10;
[quadrule_x,quadrule_w] = setupquad(quadrule_order);

loadconstants;
global E_omega E_sigma kappa;
E_inf   = 1.8;  % see note on AH thesis page 57 (65 of PDF)
E_omega = 2;
E_omega_highdiel = 15;
E_vac   = 1;
E_sigma = 80;
lambdavec = [2 11 15];  % 23 A is referenced by AH
for ki=1:length(lambdavec)
lambda=lambdavec(ki);
bigLambda = lambda * sqrt(E_inf/E_sigma);
kappa = 1/bigLambda;

prot_radius    = 13.4; % all following Gong, Hocky, and Freed (PNAS 2008)
ion_radius     = 1.4;
ion_separation = 2*ion_radius;  % simplest possible model
distance = 0:0.5:(prot_radius - ion_radius);

protFilename='../geometry/protsphere/protsphere_0.25.srf';
[meshBase, rootDir] = readsrf(protFilename);
meshData = readmesh(meshBase,1);
[centroids,normals,areas] = genmeshcolloc(meshData);
[A] = collocation_mesh(meshData,centroids,normals,areas);
[A_Y] = collocation_mesh_Yukawa(meshData,centroids,normals,areas,kappa);

% salt bridge perpendicular to radial direction
ionPairPQR = '../geometry/protsphere/ionPair.pqr';
ionPairSRF = '../geometry/protsphere/ionPair10.srf';
pqrData = readpqr(ionPairPQR);
E_nl_part1_ionpair(ki) = solvEnergyNonlocal(ionPairPQR,ionPairSRF,E_0,E_inf,E_omega,E_sigma,kappa);
E_l_part1_ionpair(ki) = solvEnergyLocal(ionPairPQR,ionPairSRF,E_0,E_inf,E_omega,E_sigma);
E_l_part1_highdiel_ionpair(ki) = solvEnergyLocal(ionPairPQR,ionPairSRF,E_0,E_inf,E_omega_highdiel,E_sigma);
for i=1:length(distance)
  pqrData.xyz(1,3) = distance(i);
  pqrData.xyz(2,3) = distance(i);

  [B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrData);


  [Anl,Bnl,Cnl]=generate_nonlocal_matrices(A,B,C,A_Y,areas,E_0, ...
														 E_inf,E_omega,E_sigma);
  xnl = gmres(Anl,Bnl*pqrData.q,[],1e-6,100);
  Enl_pair(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cnl * xnl);       % kJ/mol
  
  onecharge = [1; 0];
  xnlself = gmres(Anl,Bnl*onecharge,[],1e-6,100);
  Enl_self(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (onecharge' * Cnl * xnlself);       % kJ/mol

  [Al,Bl,Cl]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega,E_sigma);
  xl = gmres(Al,Bl*pqrData.q,[],1e-6,100);
  El_pair(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cl * xl);         % in kJ/mol

  xlself = gmres(Al,Bl*onecharge,[],1e-6,100);
  El_self(i,ki) = 0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (onecharge' * Cl * xlself);       % kJ/mol
  
  [Aref,Bref,Cref]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega,E_vac);
  xref = gmres(Aref,Bref*pqrData.q,[],1e-6,100);
  Eref_pair(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cref * xref);     % in kJ/mol

  xrefself = gmres(Aref,Bref*onecharge,[],1e-6,100);
  Eref_self(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (onecharge' * Cref * xrefself);     % in kJ/mol

  [Al,Bl,Cl]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega_highdiel,E_sigma);
  xl = gmres(Al,Bl*pqrData.q,[],1e-6,100);
  El_pair_highdiel(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cl * xl);         % in kJ/mol

  xlself = gmres(Al,Bl*onecharge,[],1e-6,100);
  El_self_highdiel(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (onecharge' * Cl * xlself);         % in kJ/mol

  [Aref,Bref,Cref]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega_highdiel,E_vac);
  xref = gmres(Aref,Bref*pqrData.q,[],1e-6,100);
  Eref_pair_highdiel(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cref * xref);     % in kJ/mol

  xrefself = gmres(Aref,Bref*onecharge,[],1e-6,100);
  Eref_self_highdiel(i,ki) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (onecharge' * Cref * xrefself);     % in kJ/mol
  
  fprintf('Done %i positions...\n',i);

end
end
save fig2

end
