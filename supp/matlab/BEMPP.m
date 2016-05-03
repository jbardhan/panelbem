loadConstants;

radius = 5.0; 

pqrData = struct('xyz',[[0 0 0];[0 4 0];[0 3.5 0];[0 3 0];[0 2.5 0];[0 2.0 0];[0 1.5 0];[0 4.0 0];[0 0 0]],'q',[0;0;0;0;0;0;0;1;0])

epsIn      = 1;
epsOut     = 80;
epsInfty   = 1.8;
lambdaVec  = 1
localLambdaVec = 1e-10;

Nmax = 30;
[gridPqrData,Ygrid,Zgrid]=  addPqrGrid(radius, struct('q',[],'xyz',[]), 0.2);

maxPot = -1e10;
minPot = 1e10;

[results_nl, minPot, maxPot] = solveEverywhereNonlocal(radius, ...
						  epsIn, epsOut, ...
						  epsInfty, lambdaVec,...
						  pqrData, ...
						  gridPqrData, Nmax, ...
						  minPot, maxPot, ...
						  Ygrid);
  

contourColors = minPot:(maxPot-minPot)/20:maxPot;

leftstart = 0.06;
topstart = 0.54;
rightstart = 0.55;
botstart = 0.09;
imwidth = 0.30;
imheight = 0.45;

figure;

subplot(1,1,1);
contourf(Ygrid,Zgrid,results_nl.gridPotentials,contourColors, ...
         'linestyle', 'none');
axis equal; caxis([minPot maxPot]); colorbar;
set(gca,'fontsize',16,'position',[leftstart botstart imwidth imheight]);
title('(a)');
