function operstruct = makeSurfaceToSurfaceYukawaPanelOperators(surf,kappa)

[V, K] = colloc_Yukawa(surf.meshData, surf.centroids, surf.normals, ...
			surf.areas,kappa);

operstruct = struct('V',V, 'K', K);
