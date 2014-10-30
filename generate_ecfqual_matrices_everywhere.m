function [A,B,C] = ...
    generate_ecfqual_matrices_everywhere(surfaceToSurfaceOperators, ...
					 chargeToSurfaceOperators, ...
					 surfaceToChargeOperators, ...
					 areas, E_0,E_inf,E_in, E_out,epsilons)
% note that the matrices A and B are missing a factor of (1/E_in)
% compared to Altman05 eqs 18-20... it will after all cancel out.

% new note: epsilons for source charges are in B now, but no 1/E_in
% in A or C!

A = surfaceToSurfaceOperators.K' * diag(areas);
A = A + 0.5 * ((E_in+E_out)/(E_in-E_out)) *diag(areas);
B = -surfaceToChargeOperators.dlpToCharges' * diag(1./epsilons); 
C = surfaceToChargeOperators.slpToCharges;
