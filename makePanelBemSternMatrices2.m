function bem = makePanelBemSternMatrices2(surf1, surf2, pqr1, pqr2, epsIn1, epsIn2, epsOut, ...
					 kappaArg)
global kappa;
kappa = kappaArg;

npDiel1 = length(surf1.dielBndy(1).areas);
npStern1 = length(surf1.sternBndy(1).areas);
npDiel2 = length(surf2.dielBndy(1).areas);
npStern2 = length(surf2.sternBndy(1).areas);
%npTotal = npDiel1 + npStern1 + npDiel2 + npStern2;

nCharges1 = length(pqr1.q);
nCharges2 = length(pqr2.q);

diel1ChargeOp = makeSurfaceToChargePanelOperators(surf1.dielBndy(1), ...
						 pqr1);
diel2ChargeOp = makeSurfaceToChargePanelOperators(surf2.dielBndy(1), ...
						 pqr2);
Idiel1 = eye(npDiel1);
Istern1 = eye(npStern1);
Idiel2 = eye(npDiel2);
Istern2 = eye(npStern2);



% src = diel1, target = diel1 
diel1Diel1Op = makeSurfaceToSurfaceLaplacePanelOperators(surf1.dielBndy(1));

% src = diel1, target = stern1
diel1ToStern1Op = makeSurfaceToSurfaceLaplacePanelOperators(surf1.dielBndy(1), ...
						  surf1.sternBndy(1));

% src = stern1, target = diel1
stern1ToDiel1Op = ...
    makeSurfaceToSurfaceLaplacePanelOperators(surf1.sternBndy(1), ...
					      surf1.dielBndy(1));
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
% src = diel2, target = diel2 
diel2Diel2Op = makeSurfaceToSurfaceLaplacePanelOperators(surf2.dielBndy(1));

% src = diel2, target = stern2
diel2ToStern2Op = makeSurfaceToSurfaceLaplacePanelOperators(surf2.dielBndy(1), ...
						  surf2.sternBndy(1));

% src = stern2, target = diel2
stern2ToDiel2Op = ...
    makeSurfaceToSurfaceLaplacePanelOperators(surf2.sternBndy(1), ...
					      surf2.dielBndy(1));
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAPLACE kernel for src = stern1, target = stern1                     
% LAPLACE kernel for src = stern2, target = stern2

stern1Stern1Laplace = makeSurfaceToSurfaceLaplacePanelOperators(surf1.sternBndy(1));
stern2Stern2Laplace = makeSurfaceToSurfaceLaplacePanelOperators(surf2.sternBndy(1));

% YUKAWA kernel for src = stern1, target = stern1
% YUKAWA kernel for src = stern2, target = stern2

if kappa > 1e-6
  stern1Stern1Yukawa = ...
      makeSurfaceToSurfaceYukawaPanelOperators(surf1.sternBndy(1), kappa);
  stern2Stern2Yukawa = ...
      makeSurfaceToSurfaceYukawaPanelOperators(surf2.sternBndy(1), kappa);
  stern1Stern2Yukawa = ...
      makeSurfaceToSurfaceYukawaPanelOperators(surf1.sternBndy(1), kappa, ...
      surf2.sternBndy(1));
  stern2Stern1Yukawa = ...
      makeSurfaceToSurfaceYukawaPanelOperators(surf2.sternBndy(1), kappa, ...
      surf1.sternBndy(1));
else
  stern1Stern1Yukawa = ...
      makeSurfaceToSurfaceLaplacePanelOperators(surf1.sternBndy(1));
  stern2Stern2Yukawa = ...
      makeSurfaceToSurfaceLaplacePanelOperators(surf2.sternBndy(1));
  stern1Stern2Yukawa = ...
      makeSurfaceToSurfaceLaplacePanelOperators(surf1.sternBndy(1), surf2.sternBndy(1));
  stern2Stern1Yukawa = ...
      makeSurfaceToSurfaceLaplacePanelOperators(surf2.sternBndy(1), surf1.sternBndy(1));
end

ZD1D2 = zeros(npDiel1, npDiel2);
ZD1S1 = zeros(npDiel1, npStern1);
ZD1S2 = zeros(npDiel1, npStern2);
ZD2D1 = zeros(npDiel2, npDiel1);
ZD2S1 = zeros(npDiel2, npStern1);
ZD2S2 = zeros(npDiel2, npStern2);
ZS1D1 = zeros(npStern1, npDiel1);
ZS1D2 = zeros(npStern1, npDiel2);
ZS1S2 = zeros(npStern1, npStern2);
ZS2D1 = zeros(npStern2, npDiel1);
ZS2D2 = zeros(npStern2, npDiel2);
ZS2S1 = zeros(npStern2, npStern1);
% could just use Z' but this is clearer

%A
A11 = (Idiel1/2 + diel1Diel1Op.K);
A12 = -diel1Diel1Op.V;
A13 = ZD1D2;
A14 = ZD1D2;
A15 = ZD1S1;
A16 = ZD1S1;
A17 = ZD1S2;
A18 = ZD1S2;

A21 = (Idiel1/2 - diel1Diel1Op.K);
A22_base = diel1Diel1Op.V;
A23 = ZD1D2;
A24 = ZD1D2;
A25 = stern1ToDiel1Op.K;
A26 = -stern1ToDiel1Op.V;
A27 = ZD1S2;
A28 = ZD1S2;

A31 = ZD2D1;
A32 = ZD2D1;
A33 = (Idiel2/2+diel2Diel2Op.K);
A34 = -diel2Diel2Op.V;
A35 = ZD2S1;
A36 = ZD2S1;
A37 = ZD2S2;
A38 = ZD2S2;

A41 = ZD2D1;
A42 = ZD2D1; 
A43 = (Idiel2/2-diel2Diel2Op.K);
A44_base = diel2Diel2Op.V;
A45 = ZD2S1; 
A46 = ZD2S1; 
A47 = stern2ToDiel2Op.K; 
A48 = -stern2ToDiel2Op.V; 

A51 = -diel1ToStern1Op.K;
A52_base = diel1ToStern1Op.V;
A53 = ZS1D2;
A54 = ZS1D2;
A55 = (Istern1/2+stern1Stern1Laplace.K);
A56 = -stern1Stern1Laplace.V;
A57 = ZS1S2;
A58 = ZS1S2;


A61 = ZS1D1;
A62 = ZS1D1;
A63 = ZS1D2;
A64 = ZS1D2;
A65 = (Istern1/2-stern1Stern1Yukawa.K);
A66_base = stern1Stern1Yukawa.V; 
A67 = -stern2Stern1Yukawa.K;
A68_base = stern2Stern1Yukawa.V;


A71 = ZS2D1;
A72 = ZS2D1;
A73 = -diel2ToStern2Op.K;
A74_base = diel2ToStern2Op.V;
A75 = ZS2S1;
A76 = ZS2S1;
A77 = (Istern2/2+stern2Stern2Laplace.K);
A78 = -stern2Stern2Laplace.V;


A81 = ZS2D1;
A82 = ZS2D1;
A83 = ZS2D2;
A84 = ZS2D2;
A85 = -stern1Stern2Yukawa.K;
A86_base = stern1Stern2Yukawa.V;
A87 = (Istern2/2-stern2Stern2Yukawa.K);
A88_base = stern2Stern2Yukawa.V;


A = [A11 A12 A13 A14 A15 A16 A17 A18;
     A21 (epsIn1/epsOut)*A22_base A23 A24 A25 A26 A27 A28;
     A31 A32 A33 A34 A35 A36 A37 A38;
     A41 A42 A43 (epsIn2/epsOut)*A44_base A45 A46 A47 A48;
     A51 (epsIn1/epsOut)*A52_base A53 A54 A55 A56 A57 A58;
     A61 A62 A63 A64 A65 (epsOut/epsOut)*A66_base A67 (epsOut/epsOut)*A68_base;
     A71 A72 A73 (epsIn2/epsOut)*A74_base A75 A76 A77 A78;
     A81 A82 A83 A84 A85 (epsOut/epsOut)*A86_base A87 (epsOut/epsOut)*A88_base];
     
B = [diel1ChargeOp.phiCoul/epsIn1;
     zeros(npDiel1,nCharges1);
     diel2ChargeOp.phiCoul/epsIn2;
     zeros(npDiel2,nCharges2);
     zeros(npStern1,nCharges1);
     zeros(npStern1,nCharges1);
     zeros(npStern2,nCharges2);
     zeros(npStern2,nCharges2)];

C = 4*pi*[-diel1ChargeOp.dlpToCharges diel1ChargeOp.slpToCharges ...
    -diel2ChargeOp.dlpToCharges diel2ChargeOp.slpToCharges ...
     zeros(nCharges1,npStern1) zeros(nCharges1, npStern1) ...    
     zeros(nCharges2,npStern2) zeros(nCharges2, npStern2)];

bem = struct('B', B, 'A', A, 'C', C,...
	     'diel1ChargeOp',diel1ChargeOp,...
	     'diel2ChargeOp',diel2ChargeOp,...
	     'diel1Diel1Op',diel1Diel1Op,...
	     'diel2Diel2Op',diel2Diel2Op,...
	     'diel1ToStern1Op',diel1ToStern1Op,...
	     'diel2ToStern2Op',diel2ToStern2Op,...
	     'stern1ToDiel1Op',stern1ToDiel1Op,...
	     'stern2ToDiel2Op',stern2ToDiel2Op,...
	     'stern1Stern1Laplace',stern1Stern1Laplace,...
	     'stern2Stern2Laplace',stern2Stern2Laplace,...
	     'stern1Stern1Yukawa',stern1Stern1Yukawa,...
	     'stern1Stern2Yukawa',stern1Stern2Yukawa,...
	     'stern2Stern1Yukawa',stern2Stern1Yukawa,...
	     'stern2Stern2Yukawa',stern2Stern2Yukawa,...
	     'A11',A11,'A12',A12,'A13',A13,'A14',A14,...
	     'A15',A15,'A16',A16,'A17',A17,'A18',A18,...
	     'A21',A21,'A22_base',A22_base,'A23',A23,'A24',A24,...
	     'A25',A25,'A26',A26,'A27',A27,'A28',A28,...
         'A31',A31,'A32',A32,'A33',A33,'A34',A34,...
	     'A35',A35,'A36',A36,'A37',A37,'A38',A38,...
	     'A41',A41,'A42',A42,'A43',A43,'A44_base',A44_base,...
	     'A45',A45,'A46',A46,'A47',A47,'A48',A48,...
	     'A51',A51,'A52_base',A52_base,'A53',A53,'A54',A54,...
	     'A55',A55,'A56',A56,'A57',A57,'A58',A58,...
	     'A61',A61,'A62',A62,'A63',A63,'A64',A64,...
	     'A65',A65,'A66_base',A66_base,'A67',A67,'A68_base',A68_base,...
         'A71',A71,'A72',A72,'A73',A73,'A74_base',A74_base,...
	     'A75',A75,'A76',A76,'A77',A77,'A78',A78,...
	     'A81',A81,'A82',A82,'A83',A83,'A84',A84,...
	     'A85',A85,'A86_base',A86_base,'A87',A87,'A88_base',A88_base);
