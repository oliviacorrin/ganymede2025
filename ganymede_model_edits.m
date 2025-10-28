% Analytical model to evaluate reconnection onset conditions at the
% magnetopause of Ganymede based on those outlined in Masters, 2014 
% ("Magnetic reconnection at Uranus’ magnetopause")
% Code developed by Olivia Corrin and Katelin Donaldson at the University of Oregon, 
% with help from prior work from Angela (last name?) 

% Current todos: 
% more accurate plasma density values (in Jbook?)
% more accurate guess for paraboloid end value
% add flaring angle to paraboloid?
% look into whether or not we can assume A = 2 (from Cooling paper)
% add variable Jovian field & density!

%% Constants
pmass = 1.67*10^-27; % (kg) proton mass
mu_0 = 4*pi*10^-7; % (H/m) vacuum permeability
emass = 9.109 * 10^-31; % (kg) electron mass
q = -1.602*10^-19; % (C) electron charge
flowspeed = 139*10^3; % (m/s) speed of Jovian plasma 
kb = 1.38*10^-23; % (JK-1) Boltzmann's constant
gamma = 5/3; % (dimensionless) used in sonic mach number calculation

%% Ganymede Parameters
R = 2.634*10^6; % (m) Ganymede radius 
b_ganymede_eq = 719*10^-9; % (T) Ganymede equatorial dipole field strength 
MP = 2.2; % magnetopause standoff distance in planetary radii
limitplotting = 8; % Random guess for where we want the paraboloid to end
ganymede_surface_den = 0.002 * 10^6 * pmass; % (kg/m^3) need Ganymede value. plasma mass density inside Ganymede's magnetosphere
alfven_num = 0.62; % Alfven Mach number
ganymede_ionlength = sqrt(pmass/(ganymede_surface_den*mu_0*q^2));
ganymede_beta = 1; % placeholder
 
%% Upstream Parameters
% Assume B_z is only nonzero component of upstream field (in plasma sheet)
upstream_Bmag = 1e-6; % (T) Jovian external field magnitude in plasma sheet
upstream_Bz = upstream_Bmag; 
% Set other 2 components essentially to 0
upstream_Bx = 10e-200; % 0
upstream_By = 10e-200; % 0
jupiter_den = 0.03 * 10^6 * pmass; % (kg/m^3) come back to this. plasma mass density on the Jovian side of the magnetopause
upstream_ionlength = sqrt(pmass/(jupiter_den*mu_0*q^2));
upstream_beta = 1; % placeholder 

%%  Set up volume using meshgrid
% meshgrid takes 2 input vectors & generates 2 2D matrices representing 
% coordinates of a grid w/ each element corresponding to a specific (x,y) 
% all grid spaces are in units of one planetary radii
n = 101; % resolution (number of points along each axis)

% Create 1D vectors of x,y,z coordinates to define spatial extent of the
% volume grid 
xi = linspace(-MP-1,limitplotting,n); 
yi = linspace(-8,8,n);
zi = linspace(-8,8,n);
% Create the corresponding arrays  
[xg,yg,zg] = meshgrid(xi,yi,zi); % creates three 3D matrices xg,yg,zg containing x,y,z coordinates for every pt. in the volume

%% Create paraboloid as magnetopause surface using meshgrid
[yii,zii] = meshgrid(yi,zi); % yii,zii are 2D matrices representing Y,Z coordinates of paraboloid surface 
paraboloid = 1/(2*MP)*(yii.^2 + zii.^2) - MP; % represent magnetopause surface as 2D surface in y-z plane with standoff distance MP calculated via pressure balance 

% Create paraboloids on either side of "magnetopause" to map the 3D grid
% values for Jovian wind and magnetosphere to. This approach is necessary so
% the curved grid spaces line up with each other in each paraboloid. 
% Initial code from Angela (past grad student) for use at Uranus: works to test values on either
% side of the boundary.
% Q: since magnetopause is smaller at Ganymede maybe we should reduce 0.3 value?
mp1 = MP - 0.3; % bounding paraboloid standoff distance (0.3 Ganymede radii in one direction from MP)
mp2 = MP + 0.3; % bounding paraboloid standoff distance (0.3 Ganymede radii in other direction from MP)
p1 = 1/(2*mp1)*(yii.^2 + zii.^2) - mp1; % first bounding paraboloid
p2 = 1/(2*mp2)*(yii.^2 + zii.^2) - mp2; % second bounding paraboloid

%% Define grid spaces where Jovian field &  magnetosphere are located.
% All spaces where grid = 0 are inside the magnetosphere 

grid = zeros(n,n,n); % create 3D matrix of zeros for each point inside the grid

% Set all grid points located outside the magnetosphere in the upstream
% Jovian field to 1
% iterate over every point (i,j,k) inside the volume 
for i = 1:n
    for j = 1:n
        for k = 1:n
            % if the x-coordinate of the grid < x-coordinate of the
            % paraboloid then you are in Jovian regime
            if xg(1,j,1) < paraboloid(i,k)
                grid(i,j,k) = 1; % Jupiter space

            % otherwise the x-coordinate of the grid => x-coordinate of the
            % paraboloid and you are in Ganymede regime 
            else
                grid(i,j,k) = 0; % Ganymede space
            end
        end
    end
end


%% Designate density and magnetic field for magnetosphere and Jovian field.
% Initialize vectors to store Jovian B-field components and density
jovian_bx = zeros(n,n,n); 
jovian_by = zeros(n,n,n);
jovian_bz = zeros(n,n,n);
jovian_btot_2d = zeros(n,n,n);
jovian_den_2d = zeros(n,n,n);

% Initialize vectors to store Ganymede B-field components and density
ganymede_bx = zeros(n,n,n);
ganymede_by = zeros(n,n,n);
ganymede_bz = zeros(n,n,n);
ganymede_den = zeros(n,n,n);

% Set initial upstream Jovian field values in 3D space
for i = 1:1:n
    for j = 1:1:n
        for k = 1:1:n
            if grid(i,j,k) == 0 % if magnetosphere, set Jovian values to zero at these spots in grid
                jovian_bz(i,j,k) = 0;
                jovian_bx(i,j,k) = 0;
                jovian_by(i,j,k) = 0;
                jovian_den_2d(i,j,k) = 0;
            else % if Jovian field, set Jovian values in grid
                jovian_bz(i,j,k) = upstream_Bz;
                jovian_by(i,j,k) = upstream_By;
                jovian_bx(i,j,k) = upstream_Bx;
                jovian_btot_2d(i,j,k) = sqrt(jovian_bz(i,j,k)^2 + jovian_by(i,j,k)^2 + jovian_bx(i,j,k)^2);
                jovian_den_2d(i,j,k) = jupiter_den;
            end
        end
    end
end

% Set initial Ganymede field values in 3D space
% This is probably not right because there's an induced dipole but just as
% a placeholder for now 
% 9/23: induced dipole is pretty small (pretty much just changing y component) 
% compared to overall dipole & not in z - start without and go from there
% (from Carol)
ganymede_bx = ((-3*b_ganymede_eq*xg.*zg))./((sqrt(xg.^2+yg.^2+zg.^2)).^5);
ganymede_by = ((-3*b_ganymede_eq*yg.*zg))./((sqrt(xg.^2+yg.^2+zg.^2)).^5);
ganymede_bz = ((b_ganymede_eq*(-2*zg.^2+xg.^2+yg.^2)))./((sqrt(xg.^2+yg.^2+zg.^2)).^5);
% this is a mask so the values are zero if not inside magnetopause 
ganymede_bx(grid~=0) = 0;
ganymede_by(grid~=0) = 0;
ganymede_bz(grid~=0) = 0;

%% Employ approach from Cooling et al. (2001) to drape the Jovian field lines
% over the magnetopause nose
l = zeros(n,n,n); % will temporarily store x-coordinate of every point in the grid
% Initialize B-field draping vectors to represent the upstream field values
% directly outside the magnetopause being bent around the magnetopause
% boundary
drape_bx = zeros(n,n,n);
drape_by = zeros(n,n,n);
drape_bz = zeros(n,n,n);
drape_tot = zeros(n,n,n);
A = 2; % not sure if this works for a subsonic flow, but the paper says A is usually 2

% Iterate over the grid to drape the Jovian field over the magnetopause
% boundary as in Cooling et al.
for i = 1:1:n
    for j = 1:1:n
        for k = 1:1:n
            l(i,j,k) = xg(i,j,k);
            drape_bx(i,j,k) = -A*(-jovian_bx(i,j,k)*(1-MP/(2*l(i,j,k))) + jovian_by(i,j,k)*-1*yg(i,j,k)/l(i,j,k) + jovian_bz(i,j,k)*zg(i,j,k)/l(i,j,k) );
            drape_by(i,j,k) =  A*(-jovian_bx(i,j,k)*-1*yg(i,j,k)/(2*l(i,j,k)) + jovian_by(i,j,k)*(2-(-1*yg(i,j,k))^2/(l(i,j,k)*MP) ) - jovian_bz(i,j,k)*-1*yg(i,j,k)*zg(i,j,k)/(l(i,j,k)*MP));
            drape_bz(i,j,k) = A*(-jovian_bx(i,j,k)*zg(i,j,k)/(2*l(i,j,k)) + jovian_by(i,j,k)*-1*yg(i,j,k)*zg(i,j,k)/(l(i,j,k)*MP) + jovian_bz(i,j,k)*(2 - (zg(i,j,k))^2/(l(i,j,k)*MP)));

            % Modify magnitude for space where x > 0 
            if xg(i,j,k) > 0
                drape_bz(i,j,k) = -1*drape_bz(i,j,k);
                drape_bx(i,j,k) = -1*drape_bx(i,j,k);
                drape_by(i,j,k) = -1*drape_by(i,j,k);
            end
            
            % Calculate norm of the draped field vector
            drape_tot(i,j,k) = sqrt(drape_bx(i,j,k)^2 + drape_by(i,j,k)^2 + drape_bz(i,j,k)^2);

            % Catch for if the draped field total is somehow larger than
            % the original upstream field
            if drape_tot(i,j,k) > jovian_btot_2d(i,j,k)
                drape_bx(i,j,k) = (drape_bx(i,j,k)/drape_tot(i,j,k)) * jovian_btot_2d(i,j,k);
                drape_by(i,j,k) = (drape_by(i,j,k)/drape_tot(i,j,k)) * jovian_btot_2d(i,j,k);
                drape_bz(i,j,k) = (drape_bz(i,j,k)/drape_tot(i,j,k)) * jovian_btot_2d(i,j,k);
            end
        end
    end
end


%% Assign 3D density and magnetic field values to the 2D paraboloid surface
% Because the reconnection criteria happen only across the current sheet,
% we need to isolate those values that lie only on either side of this
% sheet.
% Initialize matrices only using boundary values (values that lie right
% along the magnetopause border) 

% 2D boundary arrays mapping the field strength at MP surface only.
jovian_den_2d = zeros(n,n); % 2D Jovian plasma density array
jovian_bz_2d = zeros(n,n); % 2D Jovian B field z component array
jovian_by_2d = zeros(n,n); % 2D Jovian B field y component array
jovian_bx_2d = zeros(n,n); % 2D Jovian B field x component array
jovian_btot_2d = zeros(n,n); % 2D Jovian B field magnitude array

ganymede_den_2d = zeros(n,n); % 2D Ganymede plasma density array
ganymede_bx_2d = zeros(n,n); % 2D Ganymede B field x component array
ganymede_by_2d = zeros(n,n); % 2D Ganymede B field y component array
ganymede_bz_2d = zeros(n,n); % 2D Ganymede B field z component array

% Test variables to store absolute difference from each grid point to p1 and p2
% This helps to choose which 3D grid points actually lie along the 2D paraboloid 
% so we can define from the thin sheet across which reconnection happens
test_p1 = zeros(n,n,n);
test_p2 = zeros(n,n,n);

% Iterate over each 3D grid point
for i = 1:n
    for j = 1:n
        for k = 1:n

            % Difference between grid x-coord. and the offset paraboloids
            % p1,p2
            test_p1(i,j,k) = abs(xg(1,j,1) - p1(i,k));
            test_p2(i,j,k) = abs(xg(1,j,1) - p2(i,k));

            % Set difference value so we are only taking B-field values
            % from a thin layer close to the magnetopause 
            % NEED to change this - 1.50 is way too big of a difference to
            % see reconnection effects
            difference = 1.50

            if grid(i,j,k) == 0 % magnetosphere space
                ganymede_den_2d(i,k) = ganymede_surface_den; % constant density in space inside the magnetosphere
                % If we're close enough to the p1 then set the
                % boundary values to Ganymede B-field
                if test_p1(i,j,k) < difference
                    ganymede_bx_2d(i,k) = ganymede_bx(i,j,k);
                    ganymede_by_2d(i,k) = ganymede_by(i,j,k);
                    ganymede_bz_2d(i,k) = ganymede_bz(i,j,k);
                end
            
            % Note: the goal is to eventually test multiple Jovian density & field
            % values depending on location w.r.t. the plasma sheet. This is
            % an initial implementation that assumes constant density in
            % space.
            elseif grid(i,j,k) == 1 % Jupiter space
                jovian_den_2d(i,k) = jupiter_den; % Jovian boundary density can be set equal to Jovian density
                % If we're close enough to p2 then set the
                % boundary values to Jovian B-field
                if test_p2(i,j,k) < difference
                    jovian_bz_2d(i,k) = drape_bz(i,j,k);
                    jovian_by_2d(i,k) = drape_by(i,j,k);
                    jovian_bx_2d(i,k) = drape_bx(i,j,k);
                    jovian_btot_2d(i,k) = sqrt(jovian_bx_2d(i,k)^2 + jovian_by_2d(i,k)^2 + jovian_bz_2d(i,k)^2);
                end

            end

        end

    end
end

%% Calculate the left hand side of diamagnetic drift condition (theta_b)
% Reconnection onset requires a "sub-Alfvénic relative diamagnetic drift
% between ions & electrons within the current sheet in the diretion of
% reconnection" 
% initialize
b1dotb2 = zeros(n,n); % dot product of 2 boundary fields
bmag_ganymede = zeros(n,n); % Ganymede B-field mag.
bmag_upstream = zeros(n,n); % Jovian upstream B-field mag.
leftside_1 = zeros(n,n); % vector to store magnetic shear values (theta_b)

for j = 1:n
    for k = 1:n
        
        b1dotb2(j,k) = ganymede_bx_2d(j,k)*jovian_bx_2d(j,k) + ganymede_by_2d(j,k)*jovian_by_2d(j,k) + ganymede_bz_2d(j,k)*jovian_bz_2d(j,k);
        leftside_1(j,k) = acos(b1dotb2(j,k));
        
        % Old way of calculating theta_b: revised to above on 9/30  
        % bmag_ganymede(j,k) = sqrt(bganymede_bx(j,k)^2 + bganymede_by(j,k)^2 + bganymede_bz(j,k)^2);
        % bmag_upstream(j,k) = sqrt(jovian_bx_2d(j,k)^2 + jovian_by_2d(j,k)^2 + jovian_bz_2d(j,k)^2);
        % leftside_1(j,k) = b1dotb2(j,k)/(bmag_ganymede(j,k) * bmag_upstream(j,k));
    end
end

%% Calculate the right hand side of diamagnetic drift condition
% The drift condition requires that the magnetic shear is greater than the following:
% diamagnetic_drift = 2*atan(ion_inertial_length*delta_beta)/(2*ion_inertial_length))

delta_beta = abs(upstream_beta - ganymede_beta); % change in plasma beta values between Jovian & Ganymede magnetic environments
ion_inertial_length = zeros(n,n); % ion inertial length 

for j = 1:n
    for k = 1:n 
    end
end



%% Calculate the LHS of plasma beta drift condition
% Reconnection onset requires "a sub-Alfvénic flow shear across the current
% sheet in the direction of reconnection outflow"

% Calculate flow shear in direction of reconnection outflow
flowspeed_ganymede = zeros(n,n); % plasma velocity inside the magnetosphere
flowspeed_upstream = zeros(n,n); % upstream Jovian plasma velocity 
leftside_2 = zeros(n,n); % vector to store flow shear values 





%% Calculate RHS of plasma beta drift condition

