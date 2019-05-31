function distT=distTorus(a,b,varargin)
if nargin>2
    n=varargin{1};
    m=n;
else
m=sqrt(length(b));
n=m;
end

distx=min(abs(a(1)-b(:,1)),n-abs(a(1)-b(:,1)));
disty=min(abs(a(2)-b(:,2)),m-abs(a(2)-b(:,2)));
distT=sqrt(distx.^2+disty.^2);
end