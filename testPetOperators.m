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

pqrRaw = load('../meshes/test_pet.pqr');
pqrGrid = struct('q', pqrRaw(:,4),'xyz',pqrRaw(:,1:3));

epsIn  = 80.0;
epsOut = 4.0;

epsilons = epsOut*ones(length(pqrGrid.q),1);
for i=1:length(pqrGrid.q)
  if pointIsInsideMolecule(pqrSphere, pqrGrid.xyz(i,:))
    epsilons(i) = epsIn;
  end
end

srfFile = '../meshes/sphereBoda_2.srf';

[results_bem_panel, bem_data] = computePanelBEMEverywhere(srfFile, ...
						  epsIn, epsOut, ...
						  pqrGrid, ...
						  'charge_epsilons',epsilons);

% exact normal E field at centroid; 1point quad for integral
petB = load('B_exact_field.txt'); 
% FD normal E field, 1 point quad for integral
petB_fd = load('B_FD_field.txt');

Amat

