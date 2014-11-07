function [results, panel_data] = computePanelBEM(srfFile, epsIn, ...
						 epsOut, pqrData, varargin)
global E_0 conv_factor
panel_data_defined = 0;
panel_data = 0;
P          = 'eye';

for i=1:2:length(varargin)
  if strcmp(varargin{i},'panel_data') == 1
    panel_data_defined = 1;
    panel_data = varargin{i+1};
  elseif strcmp(varargin{i},'diagP') == 1
    P = 'diag';
  elseif strcmp(varargin{i}, 'explicitP') == 1
    P = varargin{i+1};
  else 
    fprintf('Error: unknown optional argument %s\n', varargin{i});
    results = 0;
    panel_data = 0;
    return
  end
end

if panel_data_defined == 0
  [meshBase,rootDir] = readsrf(srfFile);
  meshData = readmesh(meshBase,1);
  [centroids,normals,areas] = genmeshcolloc(meshData);
  
  [A] = collocation_mesh(meshData,centroids,normals,areas);
  [B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrData);

  [A_qual,B_qual,C_qual] = generate_ecfqual_matrices(A, B, C, areas, ...
						  E_0, 1.0, epsIn, ...
						  epsOut);
  panel_data = struct('numElements',length(areas),'numCharges', ...
		      length(pqrData.q), 'meshData',meshData, ...
		      'centroids',centroids, 'normals',normals, ...
		      'areas',areas,'SurfaceToSurface',A, ...
		      'ChargeToSurface',B,'SurfaceToCharge',C, ...
		      'A',A_qual,'B',B_qual,'C',C_qual);
end

if strcmp(P,'eye') == 1
  P = eye(panel_data.numElements);
elseif strcmp(P,'diag') == 1
  P = diag(diag(panel_data.A));
end

for i=1:length(pqrData.q)
  q_temp = 0 * pqrData.q; q_temp(i) = 1;
  [AinverseB(:,i),junk1,junk2,junk3,resvec] = gmres(panel_data.A, ...
						    panel_data.B*q_temp, ...
						    [], 1e-6, 100, P);
  calculations(i) = struct('resvec',resvec,'solution',AinverseB(:,i));
end

L_bem = 4 * pi * conv_factor * panel_data.C * AinverseB;
[V, D] = eig(L_bem);

results = struct('calculation', calculations, 'L', L_bem, 'V', V, ...
		 'D', diag(D));