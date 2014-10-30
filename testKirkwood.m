addpath('kirkwood');
initializeGlobals

R = 10;
pqrData = struct('q', [1 0.5]', 'xyz',[9 0 0; 0 5 0]);
epsIn = 4;
epsOut = 80;

Nmax = [16 20 30 40 50 60 70 80 90 100];
for i = 1:length(Nmax)
  [results_exact, results_CFA] = computeKirkwoodAndBIBEE(R, epsIn, ...
						  epsOut, pqrData, ...
						  Nmax(i));
  dG_exact(i) = .5 * pqrData.q' * results_exact.L * pqrData.q;
  dG_CFA(i) = .5 * pqrData.q' * results_CFA.L * pqrData.q;
  fprintf('Done with %d of %d Nmax settings.\n',i,length(Nmax));
end

ref_dG_exact = [ -24.380127549269751
 -24.767738099242614
 -25.024929360984533
 -25.056178926902231
 -25.059976910273999
 -25.060438562530848
 -25.060494680852425
 -25.060501502838466
 -25.060502303974264
 -25.060502368128805]';

ref_dG_CFA = [ -18.170996939549788
 -18.379701424676682
 -18.517369821656438
 -18.533999469457850
 -18.536014275935401
 -18.536258702079206
 -18.536288374830765
 -18.536291978480424
 -18.536292401365234
 -18.536292435210079]';

error_exact = norm(ref_dG_exact-dG_exact);
error_CFA   = norm(ref_dG_CFA-dG_CFA);

fprintf('Error in exact calc: %f kcal/mol\n',error_exact);
fprintf('Error in CFA calc: %f kcal/mol\n',error_CFA);
