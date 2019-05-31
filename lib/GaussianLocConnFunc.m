function rsim=GaussianLocConnFunc(n,distfunc,sigmaG)
% Calculates local connectivity matrix with gaussian
% ARGS:
% n = length and width of the grid;
% distfunc = @distSheet (to calculate euclidian dist)
% sigmaG = standard deviation of the gaussian where
% G(x)=exp(-x.^2./(2*sigmaG^2)). The unit is in minicolumns, i.e. 1 means
% 50 micrometre, 2 means 100 micrometre, etc.
% RETURNS:
% rsim = sparse connectivity matrix

% width of the square
nsub=n^2;


%lay out coordinates
[coordx,coordy] = meshgrid(1:n,1:n);
coorx=reshape(coordx,nsub,1);
coory=reshape(coordy,nsub,1);

% Calculate the locations of the sparse elements
%for i=1:nsub
parfor i=1:nsub
    % call distfinc (i.e. sheet, torus etc
    distM=(distfunc([coorx(i) coory(i)],[coorx coory],n));                    %gets all the distances from all other points to current point
    p=exp(-distM.^2./(2*sigmaG^2));%generate gaussian
    p(distM==0)=0;%no self connection
    p=p/max(p);%renormalise
    pindex=find(p>0.022);%cut off unconnected ones
    pV=rand(size(pindex));
    pConn=p(pindex)>pV;
    
    
    
    %recover only nonzero entries for building sparse matrix later
    tstruct(i).locs=[pindex(pConn); i];
    tstruct(i).dists=[1./distM(pindex(pConn)); 0];
    
end


% make sparse matrix
rsim=sparse(nsub,nsub);

% put the sparse elements into the sparse matrix
for i=1:nsub
    rsim(tstruct(i).locs,i) = tstruct(i).dists;
end










end

