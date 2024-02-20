function [Ex,Ey,Ez] = muscl_field_source(Ex,Ey,Ez,Jx,Jy,Jz,grid)

% MUSCL field source:
option = "CD";

if option == "CD"
    %Retrieve Constants:
    dt = grid.dt;
    c_sq_dt = dt*grid.c*grid.c;
    mu_0_c_sq_dt = c_sq_dt*grid.mu_0;

    % Save the inital state
    Ex0 = Ex;
    Ey0 = Ey;
    Ez0 = Ez;

    %Compute the applied J
    t0 = grid.time;
    t1 = grid.time + grid.dt;
    t2 = grid.time + grid.dt/2;
    Jy_ext_0 = Jy*0.0;
    Jy_ext_1 = Jy*0.0;
    Jy_ext_2 = Jy*0.0;
    Jy_ext_0 = Jy_laser(Jy_ext_0,t0,grid);
    Jy_ext_1 = Jy_laser(Jy_ext_1,t1,grid);
    Jy_ext_2 = Jy_laser(Jy_ext_2,t2,grid);


    %Update Ex, Ey, Ez
    % SSP - Stage 1:
    Ex1 = Euler_step(Ex0,mu_0_c_sq_dt,Jx,0);
    Ey1 = Euler_step(Ey0,mu_0_c_sq_dt,Jy,Jy_ext_0);
    Ez1 = Euler_step(Ez0,mu_0_c_sq_dt,Jz,0);
    
    % SSP - Stage 2:
    Ex2 = (3/4)*Ex0 + 0.25*Euler_step(Ex1,mu_0_c_sq_dt,Jx,0);
    Ey2 = (3/4)*Ey0 + 0.25*Euler_step(Ey1,mu_0_c_sq_dt,Jy,Jy_ext_1);
    Ez2 = (3/4)*Ez0 + 0.25*Euler_step(Ez1,mu_0_c_sq_dt,Jz,0);

    % SSP - Stage 3:
    Ex = (1/3)*Ex0 + (2/3)*Euler_step(Ex2,mu_0_c_sq_dt,Jx,0);
    Ey = (1/3)*Ey0 + (2/3)*Euler_step(Ey2,mu_0_c_sq_dt,Jy,Jy_ext_2);
    Ez = (1/3)*Ez0 + (2/3)*Euler_step(Ez2,mu_0_c_sq_dt,Jz,0);


end

end


function [E_new] = Euler_step(E,mu_0_c_sq_dt,J,Jy_ext)
    E_new = E - mu_0_c_sq_dt*J - mu_0_c_sq_dt*Jy_ext;
end



function Jy = Jy_laser(Jy,t,grid)

    %Laser amplitude (t + dt/2)
    
    %ad_hoc_factor = 16*10*2.4e5;
    Jy_laser = (1/grid.dx)*sin(2*pi*grid.c*t/grid.laser1.wavelength)*(2.0/(grid.mu_0*grid.c))*grid.laser1.E_max * exp(- ((t - grid.laser1.profile_t_peak)^2) / (grid.laser1.profile_duration^2));

    %1th order interpolation to Jy
    reim = mod(grid.Nx*(grid.laser1.position-grid.xmin)/(grid.xmax - grid.xmin),1);
    nearest_i = floor( grid.Nx*(grid.laser1.position-grid.xmin)/(grid.xmax - grid.xmin) );
    if nearest_i >= 1
        Jy(nearest_i) = (1-reim)*Jy_laser + Jy(nearest_i);
        Jy(nearest_i+1) = (reim)*Jy_laser + Jy(nearest_i+1);
    end

end