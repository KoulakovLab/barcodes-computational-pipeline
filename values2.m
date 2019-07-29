%   
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%
function [val, abu] = values2(y)

    
    yy=sort(y(:));
    dy=diff(double(yy));
    ind=find([1; dy; 1]);
    val=yy(ind(1:(length(ind)-1)));
    abu=diff(ind);
    
return
