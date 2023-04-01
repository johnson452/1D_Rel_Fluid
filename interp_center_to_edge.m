function [x_interp] = interp_center_to_edge(x)
% interp_center_to_edge() Procedure
%[  A   B   C   D  ] - Centered Values
%    \ / \ / \ /
%[P   X   Y   Z   P] - Edge Values
% P = 0  - Needs to be updated elsewhere 
%   = (OLD) (A + D)/2 

%interpolate the edge cells to centered
Nx = max(size(x));
if Nx ~= 1
x_interp = zeros(1,Nx+1);
for i = 1:Nx-1
    x_interp(i+1) = (x(i) + x(i+1))/2;
end

%Handle the edge case: Periodic, but 
% is overwritten by BC is non-periodic
x_interp(1) = 0; %(x(1)+x(Nx))/2;
x_interp(Nx+1) = 0; %(x(1)+x(Nx))/2;
else
    x_interp = x;
end


end