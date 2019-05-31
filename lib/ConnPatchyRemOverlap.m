function rRemote=ConnPatchyRemOverlap(n,nM,patchSize,numPatches,remoteRad,nOut,nOverlap,distfunc,CellLocFunc)
 % Function for obtaining patchy remote connectivity matrix
 % ARGS:
 % n^2 is the number of units in the sheet
 % nM^2 is the number of units in a macrocolumn
 % patchSize is the number of units in a remote patch
 % numPatches is the number of remote patches a local macrocolumn connects to
 % remoteRad is the max distance (in number of units) of the location of the patches
 % nOut is the number of outgoing connections to remote patches per unit
 % nOverlap is the number of remote patches a local macrocolumn shares with its neighbour
 % distfunc is the function used to calculate the distances (can be @distTorus or @distSheet for toroidal or zero-flux boundaries)
 % CellLocFunc is the function that generates patches of a certain size, @makeCellCluster is used for clusters that do not wrap around the boundary, @makeCellClusterToroidal is used for clusters that wrap around the boundary
 % RETURNS:
 % rRemote = is a sparse n^2 x n^2 connectivity matrix of the remote patchy connections



numMacroCol=n/nM;%number of macrocolumns
patchPercent=patchSize/n^2;

MCmap=zeros(numMacroCol^2,patchSize,numPatches);


[xMacro,yMacro]=meshgrid(1:numMacroCol,1:numMacroCol);
xMacro=reshape(xMacro,numMacroCol^2,1);
yMacro=reshape(yMacro,numMacroCol^2,1);


[xF,yF]=meshgrid(1:n,1:n);
xF=reshape(xF,n^2,1);
yF=reshape(yF,n^2,1);

expectedCl=round(pi*(remoteRad^2 - nM^2));%expected number of neighbours

for k=randperm(numMacroCol^2)
    
    distM=(distfunc([xMacro(k) yMacro(k)],[xMacro yMacro]));
    neighbourMacroIndex=find(distM<=1&distM>0);
    
    %get the middle of the macrocolumn in terms of the full sheet posistion
    xiF=(ceil(k/numMacroCol)-1)*nM+nM/2;
    remTemp=(rem(k,numMacroCol)-1);
    if remTemp==-1
        remTemp=numMacroCol-1;
    end
    yiF=remTemp*nM+nM/2;
    distMFull=(distfunc([xiF yiF],[xF yF]));

    %get the cells that are ok as patchposition starters
    neighbourFullIndex=find(distMFull>nM & distMFull<=remoteRad);
    cl=length(neighbourFullIndex);
    realNumPatch=round(cl/expectedCl*6);
    
    
    %check if any of the index points are filled in MCmap
    NeighbourfilledIndex=find(MCmap(neighbourMacroIndex,1,1)>0);
    PatchCoordinates=zeros(patchSize,numPatches);
    if ~isempty(NeighbourfilledIndex)
        
        chosenNeighbourIndex=NeighbourfilledIndex(randi(length(NeighbourfilledIndex),1));
        choosePatch=randperm(numPatches);
        choosePatch=choosePatch(1:nOverlap);
        PatchCoordinates(:,1:nOverlap)=MCmap(neighbourMacroIndex(chosenNeighbourIndex),:,choosePatch);
        
        PatchCoordinates = generatePatchesOverlap(numPatches,realNumPatch,cl,neighbourFullIndex,xF,yF,n,patchPercent,nOverlap,PatchCoordinates,CellLocFunc);
        
        MCmap(k,:,:)=PatchCoordinates;
    else

        %generate the patches
        PatchCoordinates = generatePatchesOverlap(numPatches,realNumPatch,cl,neighbourFullIndex,xF,yF,n,patchPercent,0,PatchCoordinates,CellLocFunc);
        MCmap(k,:,:)=PatchCoordinates;
    end
    
    
    
    
    

    
end




%write results to sparse matrix
rRemote=sparse(n^2,n^2);
indFull=reshape(1:n^2,n,n);
for xcood=1:n
    for ycood=1:n
        %find matching macro column
        MacroCNumber=(ceil(xcood/nM)-1)*numMacroCol+ceil(ycood/nM);
        PossiblePos=MCmap(MacroCNumber,:,:);
        PossiblePos=PossiblePos(:);
        randPos=randperm(patchSize*numPatches);
        WireTarget=PossiblePos(randPos(1:nOut));
        WireTarget=WireTarget(WireTarget>0);%get rid of -1s
        %wire current cell to a patch.
        cind=indFull(ycood,xcood);
        rRemote(WireTarget,cind)=1;
    end
end


%% 
%imagesc(reshape(rRemote(:,149*150+80),n,n))

end

