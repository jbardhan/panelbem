function [results, point_data] = computePointBEM(R, numSurfPoints, ...
						 epsIn, epsOut, ...
						 pqrData, varargin)
global conv_factor
P = 'eye';

origin = [0 0 0];
lapDiel = makeSphereOperators(origin, R, numSurfPoints);
[chargeToSurface,surfaceToCharge] = makeSphereChargeOperators(origin, ...
						  R, numSurfPoints, ...
						  lapDiel, pqrData, epsIn);

I = eye(numSurfPoints);
Z = zeros(numSurfPoints);

% A and B have been multiplied by epsIn
A_point = lapDiel.K' * diag(lapDiel.weights);
A_point = A_point + 0.5 * ((epsIn+epsOut)/(epsIn-epsOut))*diag(lapDiel.weights);
B_point = -surfaceToCharge.doubleLayerPotential';
C_point = surfaceToCharge.singleLayerPotential/epsIn;

point_data = struct('numElements', numSurfPoints,'numCharges', ...
		    length(pqrData.q),'meshData',0, ...
		    'centroids',lapDiel.points,'normals',lapDiel.normals, ...
		    'areas',lapDiel.weights,'A',A_point,'B', ...
		    B_point,'C',C_point);


if strcmp(P,'eye') == 1
  P = eye(point_data.numElements);
elseif strcmp(P,'diag') == 1
  P = diag(diag(point_data.A));
end

for i=1:length(pqrData.q)
  q_temp = 0 * pqrData.q; q_temp(i) = 1;
  [AinverseB(:,i),junk1,junk2,junk3,resvec] = gmres(point_data.A, ...
						    point_data.B*q_temp, ...
						    [], 1e-6, 100, P);
  calculations(i) = struct('resvec',resvec,'solution',AinverseB(:,i));
end

L_bem = 4*pi*conv_factor * point_data.C * AinverseB;
[V, D] = eig(L_bem);

results = struct('calculation', calculations, 'L', L_bem, 'V', V, ...
		 'D', diag(D));
