srfFile = '/Users/jbardhan/repos/mala-sphere-bibee/meshes/sphere_0.5.srf';
[meshBase,rootDir] = readsrf(srfFile);
meshData = readmesh(meshBase,1);
[centroids,normals,areas] = genmeshcolloc(meshData);

R = 10;
min_distance = 2.5;
epsIn = 2;
epsOut = 80;

[q, xyz] = makeSphereChargeDistribution(R, min_distance);
pqrData = struct('q',q,'xyz',xyz);

[A] = collocation_mesh(meshData,centroids,normals,areas);
[B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrData);

[A_qual,B_qual,C_qual] = generate_ecfqual_matrices(A, B, C, areas, ...
						  E_0, 1.0, epsIn, epsOut);
for i=1:length(pqrData.q)
  q_temp = 0 * pqrData.q; q_temp(i) = 1;
  AinverseB(:,i) = gmres(A_qual, B_qual*q_temp, [], 1e-6, 100);
end

L_bem = C_qual * AinverseB;