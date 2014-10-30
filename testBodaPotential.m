addpath('kirkwood');
addpath('gb');
addpath('bem');
addpath('pointbem');
addpath('preconditioning');
addpath('simplefd');
initializeGlobals

% physical parameters
R = 5;  
pqrSphere = struct('q',0,'xyz',[0 0 0],'r',[R]);
gridLine = (-8.0:0.2:8.0)'; Z = 0*gridLine;
pqrGrid = struct('q',zeros(length(gridLine),1),'xyz', [Z Z gridLine]);
BodaIndex = 61; % corresponds to charge at 0, 0, 4

epsIn  = 80;
epsOut = 2;

epsilons = epsOut*ones(length(pqrGrid.q),1);
for i=1:length(pqrGrid.q)
  if pointIsInsideMolecule(pqrSphere, pqrGrid.xyz(i,:))
    epsilons(i) = epsIn;
  end
end


% numerical parameters
Nmax = 20; % order of expansion (0 = monopole, etc)
srfFile = '../meshes/sphereBoda_2.srf';

%%%%%%%%%%%%%%%%%%%% COMPUTE REACTION OPERATORS

% get L from analytical harmonics
[results_exact,G] = tryAnalytical(R, epsIn, epsOut, pqrGrid, Nmax);
phiAnal = results_exact.L(:,BodaIndex);

% get L_bem_panel via ecfqual: track statistics on gmres convergence
[results_bem_panel, bem_panel_data] = computePanelBEMEverywhere(srfFile, ...
						  epsIn, epsOut, ...
						  pqrGrid, ...
						  'charge_epsilons',epsilons);
phiBem = results_bem_panel.L(:,BodaIndex);