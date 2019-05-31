function [mMacroCol,mMacroColLFP]=meanMacroCol(n,ns,Py,LFP)
% get the mean activity and LFP of a macrocolumn

fieldNumber=1:n^2;
fieldNumber=reshape(fieldNumber,n,n)';
mMacroCol=zeros(size(Py,1),(n/ns)^2,'single');
mMacroColLFP=zeros(size(Py,1),(n/ns)^2,'single');
counter=1;
for kk=1:n/ns
    for lk=1:n/ns
        accessM=fieldNumber((kk-1)*ns+1:kk*ns,(lk-1)*ns+1:lk*ns);
        accessV=accessM(:);
        mMacroCol(:,counter)=mean(Py(:,accessV),2);
        mMacroColLFP(:,counter)=mean(LFP(:,accessV),2);
        counter=counter+1;
    end
    
end

end