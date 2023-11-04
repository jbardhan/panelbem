function operstruct = makeSurfaceToSurfaceLaplacePanelOperators(surf, ...
                surf2)
disp("beginning of makeSurfaceToSurfaceLaplacePanelOperators \n");
disp("start laplace_timer \n");
laplace_timer = tic();

if nargin < 2
  %disp("from if nargin < 2");
  [V, K, Kp] = colloc_Laplace(surf.meshData, surf.centroids, surf.normals, ...
			      surf.areas);       
else
  %disp("from else");
  [V, K, Kp] = colloc_Laplace(surf.meshData, surf2.centroids, ...
			      surf2.normals, surf2.areas);
end

operstruct = struct('V',V, 'K', K, 'Kp', Kp);
disp("end of makeSurfaceToSurfaceLaplacePanelOperators \n");
laplace_timerVal = toc(laplace_timer);
disp("laplace_timerVal: "), disp(laplace_timerVal), disp("\n");