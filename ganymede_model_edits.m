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
jovian_beta = 1; % placeholder 

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
% We can use these paraboloids to represent either side of the
% "infinitesimally" small current sheet over which reconnection takes
% place.
% Q: since magnetopause is smaller at Ganymede maybe we should reduce 0.3 value?
mp1 = MP - 0.3; % bounding paraboloid standoff distance (0.3 Ganymede radii in one direction from MP)
mp2 = MP + 0.3; % bounding paraboloid standoff distance (0.3 Ganymede radii in other direction from MP)
p1 = 1/(2*mp1)*(yii.^2 + zii.^2) - mp1; % first bounding paraboloid located slightly inside MP
p2 = 1/(2*mp2)*(yii.^2 + zii.^2) - mp2; % second bounding paraboloid loacted slightly outside MP

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
jovian_bmag_2d = zeros(n,n,n);
jovian_den_2d = zeros(n,n,n);

% Initialize vectors to store Ganymede B-field components and density
% Check to see where/if these are used or duplicated 
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
                jovian_bmag_2d(i,j,k) = sqrt(jovian_bz(i,j,k)^2 + jovian_by(i,j,k)^2 + jovian_bx(i,j,k)^2);
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
            if drape_tot(i,j,k) > jovian_bmag_2d(i,j,k)
                drape_bx(i,j,k) = (drape_bx(i,j,k)/drape_tot(i,j,k)) * jovian_bmag_2d(i,j,k);
                drape_by(i,j,k) = (drape_by(i,j,k)/drape_tot(i,j,k)) * jovian_bmag_2d(i,j,k);
                drape_bz(i,j,k) = (drape_bz(i,j,k)/drape_tot(i,j,k)) * jovian_bmag_2d(i,j,k);
            end
        end
    end
end


%% Assign 3D density and magnetic field values to the 2D paraboloid surface
% Because the reconnection criteria happen only across the current sheet,
% we need to isolate those values that lie only on either side of this
% sheet.

% Initialize 2D boundary arrays only using boundary values (values that lie right
% along the magnetopause border) 
jovian_den_2d = zeros(n,n); % 2D Jovian plasma density array
jovian_bz_2d = zeros(n,n); % 2D Jovian B field z component array
jovian_by_2d = zeros(n,n); % 2D Jovian B field y component array
jovian_bx_2d = zeros(n,n); % 2D Jovian B field x component array
jovian_bmag_2d = zeros(n,n); % 2D Jovian B field magnitude array

ganymede_den_2d = zeros(n,n); % 2D Ganymede plasma density array
ganymede_bx_2d = zeros(n,n); % 2D Ganymede B field x component array
ganymede_by_2d = zeros(n,n); % 2D Ganymede B field y component array
ganymede_bz_2d = zeros(n,n); % 2D Ganymede B field z component array
ganymede_bmag_2d = zeros(n,n); % 2D Ganymede B field magnitude array

% We need to find the closest spots on p1 and p2 to bound the current sheet analog we are
% defining.
% Calculate distance from each x coordinate to target and find the minimum
for i = 1:n
        for k = 1:n
            [~, p1_x] = min(abs(xi - p1(i,k))); % store the index of the minimum in p1_x_min - this is one x bound of the current sheet
            [~, p2_x] = min(abs(xi - p2(i,k))); % store the index of the minimum in p1_x_min - this is one x bound of the current sheet

            % Now fill in the 2D boundary arrays with the 3D field values
            % at the boundary points (each point at index i, p1_x, k)
            % Ganymede arrays
            ganymede_den_2d(i,k) = ganymede_surface_den; % constant plasma density assumed inside the magnetosphere
            ganymede_bx_2d(i,k) = ganymede_bx(i, p1_x, k);
            ganymede_by_2d(i,k) = ganymede_by(i, p1_x, k);
            ganymede_bz_2d(i,k) = ganymede_bz(i, p1_x, k);
            % Calculate Ganymede B field magnitude
            ganymede_bmag_2d(i,k) = sqrt(ganymede_bx_2d(i,k)^2 + ganymede_by_2d(i,k)^2 + ganymede_bz_2d(i,k)^2);
       
            % Jovian arrays
            jovian_den_2d(i,k) = jupiter_den; % will vary this later
            jovian_bx_2d(i,k) = drape_bx(i, p2_x, k);
            jovian_by_2d(i,k) = drape_by(i, p2_x, k);
            jovian_bz_2d(i,k) = drape_bz(i, p2_x, k);
            % Calculate Jovian B field magnitude
            jovian_bmag_2d(i,k) = sqrt(jovian_bx_2d(i,k)^2 + jovian_by_2d(i,k)^2 + jovian_bz_2d(i,k)^2);
        end
end


%% Calculate the left hand side of diamagnetic drift condition (theta_b)
% Reconnection onset requires a "sub-Alfvénic relative diamagnetic drift
% between ions & electrons within the current sheet in the diretion of
% reconnection" 

% Take the dot product of the two B field vectors and store in b1dotb2_2d
b1dotb2_2d = ganymede_bx_2d .* jovian_bx_2d + ganymede_by_2d .* jovian_by_2d + ganymede_bz_2d .* jovian_bz_2d;
b1dotb2_normalize = ganymede_bmag_2d .* jovian_bmag_2d; % need the magnitude to normalize the resulting vector

cos_theta_b = b1dotb2_2d ./ b1dotb2_normalize;
% If the product of the 2 fields is 0, make the magnitude vector
% antiparallel
cos_theta_b(b1dotb2_normalize == 0) = -1;
% Take arccos of the resulting angle 
diamagnetic_leftside = acos(cos_theta_b);
        
% Old way of calculating theta_b: revised to above on 9/30  
% bmag_ganymede(j,k) = sqrt(bganymede_bx(j,k)^2 + bganymede_by(j,k)^2 + bganymede_bz(j,k)^2);
% bmag_upstream(j,k) = sqrt(jovian_bx_2d(j,k)^2 + jovian_by_2d(j,k)^2 + jovian_bz_2d(j,k)^2);
% leftside_1(j,k) = b1dotb2(j,k)/(bmag_ganymede(j,k) * bmag_upstream(j,k));

%% Calculate the right hand side of diamagnetic drift condition
% The drift condition requires that the magnetic shear is greater than the following:
% diamagnetic_drift = 2*atan(ion_inertial_length*delta_beta)/(2*ion_inertial_length))

delta_beta = abs(jovian_beta - ganymede_beta); % change in plasma beta values between Jovian & Ganymede magnetic environments
ion_inertial_length = zeros(n,n); % ion inertial length 
diamagnetic_rightside = zeros(n,n); % vector to store diamagnetic drift values 

% For loop to iterate over 2d space and caculate ion inertial length
for j = 1:n
    for k = 1:n 
        % if we are in Ganymede space, calculate ion inertial length using
        % Ganymede plasma density
        ion_inertial_length(j,k) = sqrt((pmass)/(ganymede_den_2d(j,k)*mu_0*q^2));
        % otherwise, we are in upstream Jupiter field and need to calculate
        % ion inertial length using Jupiter plasma density 
        ion_inertial_length(j,k) = sqrt((pmass)/(jovian_den_2d(j,k)*mu_0*q^2));
    end
end

% For loop to iterate over 2D space and calculate the diamagnetic drift
for j = 1:n
    for k = 1:n
        diamagnetic_rightside(j,k) = 2*atan(ion_inertial_length(j,k)*delta_beta)/(2*ion_inertial_length(j,k));
    end
end

%% Calculate the LHS of plasma beta drift condition
% Reconnection onset requires "a sub-Alfvénic flow shear across the current
% sheet in the direction of reconnection outflow"

% Calculate flow shear in direction of reconnection outflow
flowspeed_ganymede = zeros(n,n); % plasma velocity inside the magnetosphere
flowspeed_upstream = zeros(n,n); % upstream Jovian plasma velocity 
flowshear_leftside = zeros(n,n); % vector to store flow shear values 



%% Calculate RHS of plasma beta drift condition
flowshear_rightside = zeros(n,n);

% The following calculation should use variables for Ganymede_bmag and Jovian_bmag that should be only the reconnecting components of
% the magnetic fields. Until I talk to Carol about which components these should be, using the total mag. as a placeholder.
% Start loop to calculate components of the RHS
for i = 1:n
    for k = 1:n
        flowshear_rightside(i,k) = sqrt((jovian_bmag_2d(i,k)^2 + jovian_bmag_2d(i,k)*ganymede_bmag_2d(i,k))/jovian_den_2d*mu_0)
    end
end

%% Plot regions on projection of Ganymede's surface where reconnection is  either allowed or disallowed based on the 2 conditions


