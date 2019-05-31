function PatchCoordinates = generatePatchesOverlap(numPatches,realNumPatch,cl,neighbourFullIndex,xF,yF,n,patchPercent,nOverlap,PatchCoordinates,CellLocFunc)
%subfunction used by ConnPatchyRemOverlap

        for l=nOverlap+1:numPatches
            xpatch=xF(neighbourFullIndex(randi(cl,1)));
            ypatch=yF(neighbourFullIndex(randi(cl,1)));
            [~,CellLocsVector] = CellLocFunc(1.5, patchPercent, n, ypatch, xpatch);%in this script the coordinate system is where x and y are actually in correct (x=horizontal, y = vertical) directions
            PatchCoordinates(:,l)=CellLocsVector;
        end
        
        
                
        alreadydeleted=find(PatchCoordinates(1,:)==-1);

        if numPatches>realNumPatch+length(alreadydeleted)
            numtodelete=numPatches-(realNumPatch+length(alreadydeleted));
            todelete=randperm(numPatches);
            if isempty(alreadydeleted)
                alreadydeleted=0;
            end
            todelete=todelete(~ismember(todelete,alreadydeleted));
            todelete=todelete(1:numtodelete);
            PatchCoordinates(:,todelete)=-1;
        end
        
end