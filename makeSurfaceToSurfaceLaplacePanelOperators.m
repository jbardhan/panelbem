function operstruct = makeSurfaceToSurfaceLaplacePanelOperators(surf)

[V, K, Kp] = colloc_Laplace(surf.meshData, surf.centroids, surf.normals, ...
			surf.areas);

operstruct = struct('V',V, 'K', K, 'Kp', Kp);