initializeGlobals

% physical parameters
R = 10;
minDistance = 2.5;
pqrData = makeSphereChargeDistribution(R, minDistance,...
				       'Cartesian-grid');
epsIn  = 2;
epsOut = 80;

% numerical parameters
Nmax = 20; % order of expansion (0 = monopole, etc)
srfFile = '../meshes/sphere_0.5.srf';

%%%%%%%%%%%%%%%%%%%% COMPUTE REACTION OPERATORS

% get L, L_cfa from analytical harmonics
[results_exact, results_CFA] = computeKirkwoodAndBIBEE(R, epsIn, ...
						  epsOut, pqrData, Nmax)
% get L_gb_cfa from L exact
results_gb = computeGB(epsIn, epsOut, pqrData, diag(results_exact.L));

% get L_bem_panel via ecfqual: track statistics on gmres convergence
[results_bem_panel, bem_panel_data] = computePanelBEM(srfFile, epsIn, ...
						  epsOut, pqrData);

% get L_bem_point via ecfqual: track statistics on gmres convergence

%%%%%%%%%%%%%%%%%%%% PLOT COMPARISONS OF REACTION OPERATORS


%%%%%%%%%%%%%%%%%%%% TRY DIFFERENT-RANK GB BEM PRECONDITIONERS

V_reference = V_exact;  % choose any of the Vs from above
numVectorsList = [1, 4, 9, 16];
for i=1:length(numVectorsList)
end

%%%%%%%%%%%%%%%%%%%% PLOT SPECTRA OF UN/PRECONDITIONED SYSTEMS

for i=1:length(numVectorsList)
  rb_panel = results_panel(i);
  rb_point = results_point(i);

  figure; set(gca,'fontsize',16); hold on;
  plot(results_bem_panel.D, 'bo', 'linewidth', 2);
  plot(results_bem_point.D, 'bs', 'linewidth', 2);
  plot(rb_panel.D , 'rx', 'linewidth', 2);
  plot(rb_point.D , 'r+', 'linewidth', 2);
  legend('BEM panel','BEM point', 'BEM panel w/ preconditioner',...
	 'BEM point w/ preconditioner');
  titleString = sprintf('Spectra with %d-rank preconditioner', numVectorsList(i));
  title(titleString);
  xlabel('Real (\lambda)');
  ylabel('Imag (\lambda)');
end


%%%%%%%%%%%%%%%%%%%% PLOT RESULTS OF CALCULATIONS USING GB PRECOND

for i=1:length(numVectorsList)
  
  rb_panel = results_rb_panel(i);
%  rb_point = results_point(i);

  figure; set(gca,'fontsize',16);
  for j=1:length(q)
    semilogy(results_bem_panel.calculation(j).resvec,'b-','linewidth',1.8); hold on;
%    semilogy(results_bem_point.calculation(j).resvec,'r-.','linewidth',1.8);
    semilogy(rb_panel.calculation(j).resvec,'k--','linewidth',1.8);
%    semilogy(rb_point.calculation(j).resvec,'m:','linewidth',1.8);
  end
  legend('BEM panel','BEM point','BEM panel w/ preconditioner',...
	 'BEM point w/ preconditioner');
  titleString = sprintf('Convergence using %d-rank preconditioner', ...
			numVectorsList(i));
  title(titleString)
  xlabel('Iteration');
  ylabel('Residual');

end

%%%%%%%%%%%%%%%%%%%%

% down the road:
% get L_bem_* (both panel and point) from: {Green, Hildebrandt_local}



