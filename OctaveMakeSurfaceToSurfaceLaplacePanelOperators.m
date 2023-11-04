% The parallelized Octave version of makeSurfaceToSurfaceLaplacePanelOperators
function operstruct = OctaveMakeSurfaceToSurfaceLaplacePanelOperators(surf, ...
                surf2)
disp("beginning of OctaveMakeSurfaceToSurfaceLaplacePanelOperators \n");
disp("start Octave_laplace_timer \n");
Octave_laplace_timer = tic();
np = size(surf.meshData.face,1);

if nargin < 2
  % parallelized section
  [V, K, Kp] = pararrayfun(nproc, @(i) Octave_colloc_Laplace(surf.meshData, ... 
            surf.centroids, surf.normals, surf.areas, i, np), 1:np, ...
            "UniformOutput", false);  
else
  % TODO: this part not parallelized since the test driver doesn't enter this 
  % part of the if statement, so not sure how to develop and test
  [V, K, Kp] = colloc_Laplace(surf.meshData, surf2.centroids, ...
			      surf2.normals, surf2.areas, i);
end

% for each variable V, K, and Kp, pararrayfun gives cell of doubles
% but the operstruct variable needs to be a struct of arrays, so
% turn the cells of doubles into arrays before combining into one struct
newV = cell2mat(V);
newK = cell2mat(K);
newKp = cell2mat(Kp);
operstruct = struct('V', newV, 'K', newK, 'Kp', newKp);

disp("end of OctaveMakeSurfaceToSurfaceLaplacePanelOperators \n");
Octave_laplace_timerVal = toc(Octave_laplace_timer);
disp("Octave_laplace_timerVal: "), disp(Octave_laplace_timerVal), disp("\n");