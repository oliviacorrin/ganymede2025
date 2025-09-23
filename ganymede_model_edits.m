%% Constants
pmass = 1.67*10^-27; % (kg) proton mass
mu_0 = 4*pi*10^-7; % (H/m) vacuum permeability
emass = 9.109 * 10^-31; % (kg) electron mass
q = -1.602*10^-19; % (C) electron charge
flowspeed = 139*10^3; % (m/s) speed of Jovian plasma 
kb = 1.38*10^-23; % (JK-1) Boltzmann's constant
gamma = 5/3; % (dimensionless) used in sonic mach number calculation

%% Ganymede Parameters
R = 2.634*10^6;
b_g = 719*10^-9; % (T) Ganymede equatorial B field strength 
MP = 2.2; % magnetopause standoff distance in planetary radii
limitplotting = 8; % Random guess for where we want the paraboloid to end
surface_den = 0.002 * 10^6 * pmass; % need ganymede value 
density_jov = 0.03 * 10^6 * pmass; % 
alfven_num = 0.62; % Alfven Mach number

upstream_Bmag = 1e-6;
upstream_Bx = 10e-200;
upstream_By = 10e-200;
upstream_Bz = upstream_Bmag;

%%  Set up volume
% all grid spaces are in units of one planetary radii
n = 101; % resolution
xi = linspace(-MP-1,limitplotting,n); 
yi = linspace(-8,8,n);
zi = linspace(-8,8,n);
[xg,yg,zg] = meshgrid(xi,yi,zi);

%% Create paraboloid as magnetopause surface
[yii,zii] = meshgrid(yi,zi);
paraboloid = 1/(2*MP)*(yii.^2 + zii.^2) - MP;

% Create paraboloids on either side of "magnetopause" to map the 3D grid
% values for solar wind and magnetosphere to. This approach is necessary so
% the curved grid spaces line up with each other in each paraboloid. 
% Angela made this but i guess it works to test the values either side of the
% boundary - since magnetopause is smaller at ganymede maybe we should
% reduce 0.3 value?
mp1 = MP - 0.3;
mp2 = MP + 0.3;
p1 = 1/(2*mp1)*(yii.^2 + zii.^2) - mp1;
p2 = 1/(2*mp2)*(yii.^2 + zii.^2) - mp2;

yp = yii;
zp = zii;

%% Define grid spaces where solar wind and magnetosphere are located.
grid = zeros(n,n,n);
for i = 1:n
    for j = 1:n
        for k = 1:n
            if xg(1,j,1) < paraboloid(i,k)
                grid(i,j,k) = 1; % solar wind space
            else
                grid(i,j,k) = 0; % magnetosphere space
            end
        end
    end
end


%% Designate density and magnetic field for magnetosphere and solar wind
% Initialize vectors
jovian_bx = zeros(n,n,n);
jovian_by = zeros(n,n,n);
jovian_bz = zeros(n,n,n);
jovian_btot = zeros(n,n,n);
jovian_den = zeros(n,n,n);

l = zeros(n,n,n);
drape_bx = zeros(n,n,n);
drape_by = zeros(n,n,n);
drape_bz = zeros(n,n,n);
drape_tot = zeros(n,n,n);

ganymede_bx = zeros(n,n,n);
ganymede_by = zeros(n,n,n);
ganymede_bz = zeros(n,n,n);
ganymede_den = zeros(n,n,n);

% Set initial upstream Jovian values in 3D space
for i = 1:1:n
    for j = 1:1:n
        for k = 1:1:n
            if grid(i,j,k) == 0 % if magnetosphere, set Jovian values to zero at these spots in grid
                jovian_bz(i,j,k) = 0;
                jovian_bx(i,j,k) = 0;
                jovian_by(i,j,k) = 0;
                jovian_den(i,j,k) = 0;
            else % if Jovian field, set Jovian values in grid
                jovian_bz(i,j,k) = upstream_Bz;
                jovian_by(i,j,k) = upstream_By;
                jovian_bx(i,j,k) = upstream_Bx;
                jovian_btot(i,j,k) = sqrt(jovian_bz(i,j,k)^2 + jovian_by(i,j,k)^2 + jovian_bx(i,j,k)^2);
                jovian_den(i,j,k) = density_jov;
            end
        end
    end
end


% Employ approach from Cooling et al (2001) to drape the IMF field lines
A = 2; % I'm not sure if this works for a subsonic flow, but the paper says A is usually 2
for i = 1:1:n
    for j = 1:1:n
        for k = 1:1:n
            l(i,j,k) = xg(i,j,k);
            drape_bx(i,j,k) = -A*(-jovian_bx(i,j,k)*(1-MP/(2*l(i,j,k))) + jovian_by(i,j,k)*-1*yg(i,j,k)/l(i,j,k) + jovian_bz(i,j,k)*zg(i,j,k)/l(i,j,k) );
            drape_by(i,j,k) =  A*(-jovian_bx(i,j,k)*-1*yg(i,j,k)/(2*l(i,j,k)) + jovian_by(i,j,k)*(2-(-1*yg(i,j,k))^2/(l(i,j,k)*MP) ) - jovian_bz(i,j,k)*-1*yg(i,j,k)*zg(i,j,k)/(l(i,j,k)*MP));
            drape_bz(i,j,k) = A*(-jovian_bx(i,j,k)*zg(i,j,k)/(2*l(i,j,k)) + jovian_by(i,j,k)*-1*yg(i,j,k)*zg(i,j,k)/(l(i,j,k)*MP) + jovian_bz(i,j,k)*(2 - (zg(i,j,k))^2/(l(i,j,k)*MP)));

            % Modify for space where x > 0
            if xg(i,j,k) > 0
                drape_bz(i,j,k) = -1*drape_bz(i,j,k);
                drape_bx(i,j,k) = -1*drape_bx(i,j,k);
                drape_by(i,j,k) = -1*drape_by(i,j,k);
            end

            drape_tot(i,j,k) = sqrt(drape_bx(i,j,k)^2 + drape_by(i,j,k)^2 + drape_bz(i,j,k)^2);

            % catch for if the draped field total is somehow larger than
            % the original field
            if drape_tot(i,j,k) > jovian_btot(i,j,k)
                drape_bx(i,j,k) = (drape_bx(i,j,k)/drape_tot(i,j,k)) * jovian_btot(i,j,k);
                drape_by(i,j,k) = (drape_by(i,j,k)/drape_tot(i,j,k)) * jovian_btot(i,j,k);
                drape_bz(i,j,k) = (drape_bz(i,j,k)/drape_tot(i,j,k)) * jovian_btot(i,j,k);
            end
        end
    end
end

% Set initial Jovian values in 3D space
% This is probably not right because there's an induced dipole but just as
% a placeholder for now
ganymede_bx = ((-3*b_g*xg.*zg))./((sqrt(xg.^2+yg.^2+zg.^2)).^5);
ganymede_by = ((-3*b_g*yg.*zg))./((sqrt(xg.^2+yg.^2+zg.^2)).^5);
ganymede_bz = ((b_g*(-2*zg.^2+xg.^2+yg.^2)))./((sqrt(xg.^2+yg.^2+zg.^2)).^5);
% this is a mask so the values are zero if not inside magnetopause (had to
% ask chat gpt how to mask in matlab :( )
ganymede_bx(grid~=0) = 0;
ganymede_by(grid~=0) = 0;
ganymede_bz(grid~=0) = 0;

%% Assign 3D density and magnetic field values to the 2D paraboloid surface
% initialize matrices only using values that lie right along the magnetopause border
bjovian_den = zeros(n,n);
bjovian_bz = zeros(n,n);
bjovian_by = zeros(n,n);
bjovian_bx = zeros(n,n);
bjovian_btot = zeros(n,n);
bganymede_den = zeros(n,n);
bganymede_bx = zeros(n,n);
bganymede_by = zeros(n,n);
bganymede_bz = zeros(n,n);
test_p1 = zeros(n,n,n);
test_p2 = zeros(n,n,n);
bjovian_bzd = zeros(n,n);
bjovian_byd = zeros(n,n);
bjovian_bxd = zeros(n,n);
bjovian_btotd = zeros(n,n);

for i = 1:n
    for j = 1:n
        for k = 1:n

            % difference between position in x direction and paraboloids
            % bordering magnetopause
            test_p1(i,j,k) = abs(xg(1,j,1) - p1(i,k));
            test_p2(i,j,k) = abs(xg(1,j,1) - p2(i,k));

            difference = 1.50;

            if grid(i,j,k) == 0 % magnetosphere space
                bganymede_den(i,k) = surface_den; % constant in space inside the magnetosphere
                if test_p1(i,j,k) < difference
                    bganymede_bx(i,k) = ganymede_bx(i,j,k);
                    bganymede_by(i,k) = ganymede_by(i,j,k);
                    bganymede_bz(i,k) = ganymede_bz(i,j,k);
                end

            elseif grid(i,j,k) == 1 % solar wind space
                bjovian_den(i,k) = density_jov;
                if test_p2(i,j,k) < difference
                    bjovian_bzd(i,k) = drape_bz(i,j,k);
                    bjovian_byd(i,k) = drape_by(i,j,k);
                    bjovian_bxd(i,k) = drape_bx(i,j,k);
                    bjovian_btotd(i,k) = sqrt(bjovian_bxd(i,k)^2 + bjovian_byd(i,k)^2 + bjovian_bzd(i,k)^2);
                end

            end

        end

    end
end

%% Calculate the left hand side of diamagnetic drift condition (theta_b)
% initialize
b1dotb2 = zeros(n,n);
bmag_ganymede = zeros(n,n);
bmag_upstream = zeros(n,n);
leftside = zeros(n,n);

for j = 1:n
    for k = 1:n
        
        b1dotb2(j,k) = bganymede_bx(j,k)*bjovian_bx(j,k) + bganymede_by(j,k)*bjovian_by(j,k) + bganymede_bz(j,k)*bjovian_bz(j,k);

        bmag_ganymede(j,k) = sqrt(bganymede_bx(j,k)^2 + bganymede_by(j,k)^2 + bganymede_bz(j,k)^2);
        bmag_upstream(j,k) = sqrt(bjovian_bx(j,k)^2 + bjovian_by(j,k)^2 + bjovian_bz(j,k)^2);


        leftside(j,k) = b1dotb2(j,k)/(bmag_ganymede(j,k) * bmag_upstream(j,k));
    end
end

%% Calculate the right hand side of diamagnetic drift condition

% same thing but for plasma beta on each side 

