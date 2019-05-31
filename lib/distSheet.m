function distT=distSheet(a,b,varargin)
%calculates teh euclidian distance from the b vector to the a point
distT=sqrt((abs(b(:,1)-a(1)).^2)+(abs(b(:,2)-a(2)).^2));
end


