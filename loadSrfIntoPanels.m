function surfdata = loadSrfIntoPanels(srf)

[meshBase,rootDir] = readsrf(srf);
meshData = readmesh(meshBase,1,1);
[centroids, normals, areas] = genmeshcolloc(meshData);

surfdata = struct('meshData',meshData,'areas',areas, 'centroids',centroids, ...
		  'normals',normals);
