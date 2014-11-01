function surfdata = loadSrfIntoPanels(srf)

areaThreshold = 1e-6;

[meshBase,rootDir] = readsrf(srf);
meshData = readmesh(meshBase,1,1);
[centroids, normals, areas] = genmeshcolloc(meshData);

I = find(areas > areaThreshold);
areas = areas(I);
centroids = centroids(I,:);
normals = normals(I,:);
meshData.face=meshData.face(I,:);
meshData.X = meshData.X(:,I);
meshData.Y = meshData.Y(:,I);
meshData.Z = meshData.Z(:,I);

surfdata = struct('meshData',meshData,'areas',areas, 'centroids',centroids, ...
		  'normals',normals);
