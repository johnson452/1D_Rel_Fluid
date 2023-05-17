function [y_interp] = interp_edge_to_center(y)
% interp_edge_to_center() Procedure
%[A   B   C   D   E] - Edge Values (x)
%  \ / \ / \ / \ /
%[  X   Y   Z   S  ] - Centered Values (x2)


%interpolate the edge cells to centered
NX = max(size(y));
y = [y(NX),y,y(1)];
Nx = max(size(y));
x = linspace(0,1,Nx);
dx = x(2)-x(1);
x2 = linspace(0+3*dx/2,1-3*dx/2,Nx-3);
y_interp = interp1(x,y,x2,'spline');

% Manual interp (linear) just turrible
% Nx = max(size(x));
% x_interp = zeros(1,Nx-1);
% for i = 1:Nx-1
%     x_interp(i) = (x(i) + x(i+1))/2;
% end

end