function [y_interp] = interp_edge_to_center(y,grid)
% interp_edge_to_center() Procedure
%[A   B   C   D   E] - Edge Values (x)
%  \ / \ / \ / \ /
%[  X   Y   Z   S  ] - Centered Values (x2)


%interpolate the edge cells to centered
if grid.BC_type == "Periodic"
    NX = max(size(y));
    y = [y(NX),y,y(1)];
    Nx = max(size(y));
    x = linspace(0,1,Nx);
    dx = x(2)-x(1);
    x2 = linspace(0+3*dx/2,1-3*dx/2,Nx-3);
    y_interp = interp1(x,y,x2,'spline');
elseif grid.BC_type == "Propagation into a plasma wave beach"
    Nx = max(size(y));
    x = linspace(0,1,Nx);
    dx = x(2)-x(1);
    x2 = linspace(0+dx/2,1-dx/2,Nx-1);
    y_interp = interp1(x,y,x2,'linear');
else
    Nx = max(size(y));
    x = linspace(0,1,Nx);
    dx = x(2)-x(1);
    x2 = linspace(0+dx/2,1-dx/2,Nx-1);
    y_interp = interp1(x,y,x2,'spline'); %spline
end
end