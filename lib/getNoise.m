function Pv=getNoise(nIt,n,varargin)

if nargin==3
    rng(varargin{1});
end

pulselocs=reshape([1:n^2],n,n)';
ns=10;%size of subsquares
inputV=randn((n/ns)^2,nIt);
Pv=zeros(n^2,nIt);
m=1;
for k=1:n/ns
    for l=1:n/ns
        klocs=pulselocs((k-1)*ns+1:k*ns,(l-1)*ns+1:l*ns);
        Pv(klocs(:),:)=repmat(inputV(m,:),length(klocs(:)),1);
        m=m+1;
    end
end

end