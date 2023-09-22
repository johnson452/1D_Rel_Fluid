function [Ex,Ey,Ez,Bx,By,Bz] = external_field_fun(Ex,Ey,Ez,Bx,By,Bz,grid)

%Specify external fields on FDTD:
if grid.laser_envelope_model == "true" && grid.BC_type == "WFA"
    t = grid.time;
    m = grid.m0;
    q = grid.e0;
    c = grid.c;
    z0 = grid.laser1.position;               % This point is on the laser plane
    tau = grid.laser1.profile_duration;        % The duration of the laser (in s)
    a0 = 2.66;
    if grid.solve_type_field == "FDTD"
        x = grid.x2;
    else
        x = grid.x1;
    end
    a_env = a0*exp( -((x - c*t -z0).^2)/(c*c*tau*tau) );
    rhs = (m.*(x - c.*t -z0).*a_env.^2)./(q.*tau.*tau.*sqrt(1+0.5*a_env.^2));
    Ex = Ex + rhs;
end

end