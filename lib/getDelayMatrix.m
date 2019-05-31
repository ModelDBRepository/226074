function DD=getDelayMatrix(C,n,distfunc,hdtime)

%check if matlabpool is open and if not then open it
%if  matlabpool('size') == 0
%    matlabpool open
%end

% width of the square
nsub=n^2;


%lay out coordinates
[coordx,coordy] = meshgrid(1:n,1:n);
coorx=reshape(coordx,nsub,1);
coory=reshape(coordy,nsub,1);


%elements that need euclidian distance calculating:
[nEx,nEy]=find(C>0);
nel=length(nEx);

% Calculate the locations of the sparse elements

%D=zeros(nel,3);
parfor k=1:nel
    
    i=nEx(k);
    j=nEy(k);
    % call distfinc (i.e. sheet, torus etc
    distM=(distfunc([coorx(i) coory(i)],[coorx(j) coory(j)],n));
    D(k,:)=[i,j,distM*50];%gives the distance in micrometers
    
end



nTS=ceil(max(D(:,3))/hdtime);
DD={};
for k=1:nTS
    toconvert=D(k-1<=D(:,3)/hdtime & D(:,3)/hdtime<k,:);
    length(toconvert);
    DD{k}=sparse(toconvert(:,1),toconvert(:,2),1,nsub,nsub);
    
end




