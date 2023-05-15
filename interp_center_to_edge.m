function [y_interp] = interp_center_to_edge(y)
% interp_center_to_edge() Procedure
%[  A   B   C   D  ] - Centered Values
%    \ / \ / \ /
%[P   X   Y   Z   P] - Edge Values
% P = 0  - Needs to be updated elsewhere 
%   = (OLD) (A + D)/2 


% Pad y with perioid conditions
NX = max(size(y));
y = [y(NX),y,y(1)];
Nx = max(size(y));
x = linspace(0,1,Nx);
dx = x(2)-x(1);
x2 = linspace(0+dx/2,1-dx/2,Nx-1);
y_interp = interp1(x,y,x2,'spline');


% %interpolate the edge cells to centered
% Nx = max(size(x));
% if Nx ~= 1
% x_interp = zeros(1,Nx+1);
% for i = 1:Nx-1
%     x_interp(i+1) = (x(i) + x(i+1))/2;
% end
% 
% %Handle the edge case: Periodic, but 
% % is overwritten by BC if non-periodic
% x_interp(1) = 0; %(x(1)+x(Nx))/2;
% x_interp(Nx+1) = 0; %(x(1)+x(Nx))/2;
% else
%     x_interp = x;
% end


end