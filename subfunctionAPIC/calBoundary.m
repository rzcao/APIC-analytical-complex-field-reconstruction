function [ bound ] = calBoundary(center,mylength)
%calBoundary calculate the boundary for a box (side length L) that centers
%at point O. The first input should be the center and the second should be
%the side length of the box. When the size length is a scalar, it generate
%a square.


%   Detailed explanation goes here
% [xsize,ysize] = size(length);
% if xsize < ysize
%     length = length.';
% end
% 
% [xsize,ysize] = size(center);
% if xsize < ysize
%     center = center.';
% end

if max(size(mylength)) == 1
    mylength = ones(size(center)).*mylength;
end

bound = zeros(2,max(size(center)));
for ii = 1:max(size(center))
    len_c = floor(mylength(ii)/2+1);
    bound(1,ii) = center(ii)-len_c+1;
    bound(2,ii) = bound(1,ii)+mylength(ii)-1;
end

end

