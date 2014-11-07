addpath('kirkwood');
addpath('gb');
addpath('bem');
addpath('pointbem');
addpath('fileIO');
addpath('preconditioning');
initializeGlobals

bootstrap_source = 'gb'; % either gb or bem
numVectorsList = [1, 2, 3, 4, 5, 9, 16];
pqrData = readPqrFromPdb('../meshes/tripeptide/tripeptide.pdb');
srfFile = '../meshes/tripeptide/tripeptide_1.srf';
epsIn  = 4;    % epsilons used for ?_res1.txt. do not change
epsOut = 80;

[results_bem_panel, bem_panel_data] = computePanelBEM(srfFile, epsIn, ...
						  epsOut, pqrData);
results_gb = computeGB(epsIn, epsOut, pqrData, ...
		       diag(results_bem_panel.L));


%%%%%%%%%%%%%%%%%   GENERATE AND TEST THE RB PRECONDITIONERS

if strcmp(bootstrap_source,'gb') == 1
  V_ref = results_gb.V;
elseif strcmp(bootstrap_source,'bem') == 1
  V_ref = results_bem_panel.V;
else
  fprintf('bootstrap_source %s is not possible here!\n', ...
	  bootstrap_source);
end

for i=1:length(numVectorsList)
  P_rb = computeReducedBasisPreconditioner(bem_panel_data, ...
					      V_ref(:,1: numVectorsList(i)));
  results_rb_panel(i) = computePanelBEM(srfFile, epsIn, epsOut, ...
					pqrData, 'panel_data', bem_panel_data,...
					'explicitP',P_rb);
end

%%%%%%%%%%%%%%%%%   AS A REFERENCE, THE SIMPLE DIAGONAL PRECONDITIONER

P_diag = diag(diag(bem_panel_data.A));
results_diag_panel = computePanelBEM(srfFile, epsIn, epsOut, ...
					pqrData, 'panel_data', ...
					bem_panel_data, 'explicitP',P_diag);

%%%%%%%%%%%%%%%%%   PLOT CONVERGENCE FOR UNPRECOND, DIAG, AND RB

for i=1:length(numVectorsList)
  
  rb_panel = results_rb_panel(i);

  figure; set(gca,'fontsize',16);
  for j=1:length(pqrData.q)
    semilogy(results_bem_panel.calculation(j).resvec,'b-','linewidth',1.8); hold on;
    semilogy(results_diag_panel.calculation(j).resvec,'m-.','linewidth',1.8);
    semilogy(rb_panel.calculation(j).resvec,'k--','linewidth',1.8);
  end
  titleString = sprintf('Convergence using %d-rank preconditioner', ...
			numVectorsList(i));
  title(titleString)
  xlabel('Iteration');
  ylabel('Residual');

end
