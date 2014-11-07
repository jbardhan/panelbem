addpath('bem');
addpath('fileIO');
initializeGlobals

pqrData = readPqrFromPdb('../meshes/tripeptide/tripeptide.pdb');
srfFile = '../meshes/tripeptide/tripeptide_1.srf';
epsIn  = 4;    % epsilons used for ?_res1.txt. do not change
epsOut = 80;

[results,matlab_bem] = computePanelBEM(srfFile, epsIn, ...
						  epsOut, pqrData);

loadTripeptideResults

normA = norm(fftsvd_bem.A);
normB = norm(fftsvd_bem.B);
normC = norm(fftsvd_bem.C);
A_error = norm(matlab_bem.A/epsIn - fftsvd_bem.A);
B_error = norm(matlab_bem.B/epsIn - fftsvd_bem.B);
C_error = norm(matlab_bem.C*conv_factor*4*pi - fftsvd_bem.C);

A_relerror = A_error/normA;
B_relerror = B_error/normB;
C_relerror = C_error/normC;

fprintf('A: %f abs., %f rel.\n', A_error, A_relerror);
fprintf('B: %f abs., %f rel.\n', B_error, B_relerror);
fprintf('C: %f abs., %f rel.\n', C_error, C_relerror);

% apparent large abs error for C results from the fact that fftsvd_bem
% has been multiplied by (conversion factor * 4 * pi = 4173.4) to put
% us in kcal/mol/e for potential.  note that rel error is still
% also about 1e-3 to 1e-4