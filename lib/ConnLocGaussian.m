function rsim=ConnLocGaussian(n,distfunc,sigmaG)
% Calculates whole connectivity matrix
%this is a version that is slightly quicker than the  GaussianLocConnFunc
%at higher n>100...does the same though!
% ARGS:
% n = length and width of the grid;
% distfunc = @distSheet (to calculate euclidian dist)
% steepness = steepness of the function #conn/#all.possible.conn over
% distance
% RETURNS:
% rsim = sparse connectivity matrix

%check if matlabpool is open and if not then open it
if  matlabpool('size') == 0
    matlabpool
end

% width of the square
nsub=n^2;


%lay out coordinates
[coordx,coordy] = meshgrid(1:n,1:n);
coorx=reshape(coordx,nsub,1);
coory=reshape(coordy,nsub,1);

% Calculate the locations of the sparse elements
%for i=1:nsub
% make sparse matrix

rsim=sparse(nsub,nsub);
for i=1:nsub
    % call distfinc (i.e. sheet, torus etc
    distM=(distfunc([coorx(i) coory(i)],[coorx coory]));                    %gets all the distances from all other points to current point
    p=exp(-distM.^2./(2*sigmaG^2));%generate gaussian
    p(distM==0)=0;%no self connection
    p=p/max(p);%renormalise
    pindex=find(p>0.022);%cut off unconnected ones
    pV=rand(size(pindex));
    pConn=p(pindex)>pV;
    
    
    indTo=pindex(pConn);
    indTo=indTo(indTo~=i);
    
    rsim=rsim+sparse(indTo,i,1,nsub,nsub);

    
end

rsim(rsim>0)=1;

end

