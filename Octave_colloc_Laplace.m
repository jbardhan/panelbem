% The parallelized Octave version of colloc_Laplace
function [V,K,Kp] = Octave_colloc_Laplace(meshData,centroids,normals,areas, i, np)
panel = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
[area,collocpt,Z,fss,fds,fess] = calcp(panel,centroids,normals);

V = fss'/4/pi;
K = fds'/4/pi;
Kp = fess'/4/pi;

if mod(i,25)==0
  fprintf('Done with %d of %d columns.\n',i,np);
end

end