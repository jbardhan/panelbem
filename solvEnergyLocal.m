function [E_transfer,A_l,B_l,C_l] = ...
    solvEnergyLocal(pqrData,meshData,epsIn,epsOut)
% E_transfer is the energy of transferring from a low-dielectric
% local solvent to the high-dielectric non-local solvent
global E_0

printOn = 1;

if printOn
  fprintf('Generating discretized Laplace operators...');
end

[centroids,normals,areas] = genmeshcolloc(meshData);

% create Laplace operators
[A] = collocation_mesh(meshData,centroids,normals,areas);
[B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrData);

dummy_E_inf = 1.0;
[A_l, B_l, C_l] = generate_local_matrices(A, B, C, areas, E_0, dummy_E_inf, ...
					  epsIn, epsOut);

E_transfer = 0.5 * (Na/1000) * (q^2/E_0) * 1e10 * (pqrData.q' * C_l ...
						  * (A_l \ (B_l * ...
						  pqrData.q)))/ ...
    joulesPerCalorie;  % comes out in kcal/mol
