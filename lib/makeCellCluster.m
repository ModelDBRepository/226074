function [CellLocs,CellLocsVector] = makeCellCluster(clusterCoeff, Percent, n, varargin)
 % Function for obtaining a list of clustered locations 
 % ARGS:
 % clusterCoeff = clustering coefficient for bad cells. Must be between 0 and 1 (inclusive).
 %                0 for totally diffuse, 1 for totally clustered
 % Percent  = percent of bad cells
 % n            = Size of the sheet
 % can provide x and y location of seed in Varargin
 % RETURNS:
 % CellLocs  = Locations of the bad cells in a nxn matrix
 % CellLocsVector = list of locations in a sub2ind style


%Check that our number of bad cells isn't too high
if Percent >=1
  Percent = 1;
end

numCells=ceil(Percent*n^2);

%Sigmoid transform of cluster coeff increases realism of behaviour
%clusterCoeff = 1/(1 + exp(-11*clusterCoeff+4)) + 0.001;

%Matrix of bad cell locations
CellLocs = zeros(n,n);

%Populate bad cell locations with numBadCells bad cells.
%Note: this loop gets exponentially more inefficient as numBadCells -> n^2
%fprintf('Generating %d cells in %d total cells.\n',numCells,n^2);


%calculate the kernel of the probability function
sigma=1/clusterCoeff;
[XI,YI]=meshgrid(1:4*n,1:4*n);
gCenter=2*n;
mu=[2*n 2*n];
GKernel=Gaussian(XI,YI,mu,sigma);
clear XI YI;
GKernelCenter=GKernel(:,gCenter)>0.1*max(GKernel(:));
GKernelMask=meshgrid(GKernelCenter,GKernelCenter').*meshgrid(GKernelCenter,GKernelCenter')';
GKernelMask=logical(GKernelMask);
GKernel=GKernel(GKernelMask);
GKernel=reshape(GKernel,nnz(GKernelCenter),nnz(GKernelCenter));
%imagesc(GKernel);



if nargin>3
    xloc=varargin{1};
    yloc=varargin{2};
else
    xloc = ceil(n*rand);
    yloc = ceil(n*rand);
end
CellLocs(xloc,yloc) = 1;




for k=1:numCells-1

   %generate probability dist
   pDistr=conv2(CellLocs,GKernel,'same');
   pDistr(CellLocs==1)=0;
   pDistr=pDistr/sum(pDistr(:));
   %figure(1);
   %imagesc(pDistr);
   
   
   
   %Set the cell
   pVector=cumsum(pDistr(:));
   indNewLoc=find(pVector>=rand,1,'first');
   [xloc,yloc] = ind2sub([n,n], indNewLoc);
   CellLocs(xloc,yloc) = 1;
   


end
locs=1:n^2;
CellLocsVector=locs(reshape(CellLocs>0,n^2,1));


    
end