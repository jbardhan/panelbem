% The parallelized Octave version of makePanelBemEcfQualMatrices
function bem = OctaveMakePanelBemEcfQualMatrices(surf, pqr, epsIn, epsOut)
disp("beginning of OctaveMakePanelBemEcfQualMatrices");
epsHat = (epsIn+epsOut)/(epsIn-epsOut);
surfsurfop = OctaveMakeSurfaceToSurfaceLaplacePanelOperators(surf);
chargesurfop = makeSurfaceToChargePanelOperators(surf, pqr);

% the two factors below 
B = -1*    chargesurfop.dlpToCharges' / epsIn;
C = 4*pi*  chargesurfop.slpToCharges;

A = surfsurfop.K'*diag(surf.areas) + diag(surf.areas)*epsHat/2;

bem = struct('B', B, 'A', A, 'C', C,'surfsurfop',surfsurfop,'chargesurfop',chargesurfop);
disp("end of OctaveMakePanelBemEcfQualMatrices");
