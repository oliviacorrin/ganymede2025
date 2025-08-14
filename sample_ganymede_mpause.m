% Description: An analytical model to evaluate the possibility of
% reconnection onset at Ganymede's magnetopause with variation in upstream
% Jovian plasma sheet conditions (plasma velocity, plasma density, magnetic 
% field strength)

% Initial Goal: Use a constant set of upstream parameters
% Eventual Goal: Vary upstream parameters (at least 5 iterations: above,
% intermediate above, inside, intermediate below, and below Jovian plasma
% sheet)

clc
close all 
clear all

% define constants 
pmass = 1.67*10^-27; % (kg) proton mass
emass = 9.109 * 10^-31; % (kg) electron mass
q = -1.602*10^-19; % (C) electron charge
mu_0 = 4*pi*10^-7; % (H/m) vacuum permeability
k_b = 1.38*10^-23; % (JK-1) Boltzmann's constant
gamma = 5/3; % (dimensionless) used in sonic mach number calculation

% define upstream parameters
v_j = 139*10^3; % (m/s) Jovian plasma velocity
% v_dir_j = 3 component Jovian plasma velocity vector
alfven_num = 0.62 % (dimless) Alfven mach number
pden_j = 56 % (amu/cm^3) Jovian plasma density - number density * pmass 
% b_j = [] % (T) Jovian magnetic field strength
% b_dir_j = [] % 3 component vector of Jovian magnetic field direction (in GPhiO?)
% plasma_beta_j % [dimless] plasma beta value at Jupiter

% define moon parameters
% right now these are constant; later some will be varied to match up with
% different orientations of the moon in the PS.
r_g = 2.634*10^6 % (m) Ganymede radius
b_g = 719*10^-9 % (T) Ganymede equatorial B field strength 
% b_dir_g = % 3 component Ganymede B field direction vector
% rho_g = % (amu/cm^3) Ganymede plasma density 
% v_g = Ganymede plasma velocity
% v_dir_g = 3 dimensional Ganymede plasma velocity vector 
% plasma_beta_g = % [dimless] plasma beta value at Ganymede

% This will eventually be calculated in the code for varying upstream
% parameters(?), but for now, using a constant
standoff_distance = 2.2 % (in r_g) magnetopause standoff distance IN plasma sheet using jovianden = 56 amu/cm^3

%% Volume 
% Currently using a plasma density value of 56 amu/cm^3 corresponding to a
% position inside the Jovian plasma sheet. Eventually plan to add other
% values corresponding to different locations wrt the PS.

% get (x, y, z) GPhiO coordinates and use to create paraboloid meshgrid
% goal: to extract these from data (using rough estimates for now)
n = 101 % spacing of y,z vectors - increase for a finer mesh
x_init = linspace(-standoff_distance-1,10,n)
y_init = linspace(-10,10,n)
z_init = linspace(-10,10,n)
% set up volume
[x_vol,y_vol,z_vol] = meshgrid(x_init,y_init,z_init)

% paraboloid magnetopause surface - from Donaldson et al. 2024, derived
% from Cooling et al. 2001:
% Y^2 + Z^2 = 2mp(X+mp)
% rearrange to find X:
% (Y^2 + Z^2 - mp)/2mp = X
[y_parab,z_parab] = meshgrid(y_init,z_init)
paraboloid = (y_parab.^2 + z_parab.^2 - standoff_distance)/2*standoff_distance

% paraboloids on either side of mpause 
% mp1 = mp - 0.3
% mp2 = mp + 0.3
% p1 = (Y2.^2 + Z2.^2 - mp1)/2*mp1
% p2 = (Y2.^2 + Z2.^2 - mp2)/2*mp2

% plot using surf(x,y,z) function, which represents each value in matrix Z as the height of a surface above x-y plane.
figure
magnetopause = surf(paraboloid,y_parab,z_parab)
title('Paraboloid model of Ganymede magnetopause with standoff distance 2.2 R_G')
xlabel('X (R_G)')
ylabel('Y (R_G)')
zlabel('Z (R_G)')
grid on

% find normal vector to the paraboloid
% [nx,ny,nz] = surfnorm(paraboloid,Y,Z)

% Now check each conditional to determine if reconnection is possible or
% prohibited.
% Onset condition 1: diamagnetic drift condition. 

% Before calculations, double-check values so we don't divide by zero
% if v_g == 0 || v_j == 0 || b_g == 0 || b_j == 0 
%   disp('Cannot divide by zero')
%   return
% end

% Find the ion inertial length to use in calculations.
% First convert plasma density (amu/cm^-3) to protons number density (m^-3)
% assuming H+ (can change later)
hplus_num_density = rho_g*1.66*10^-21/1.00727647 % [m^-3]
% oplus_num_density = [] Goal: add in O+?
ion_inertial_length = sqrt(pmass/(hplus_num_density*mu_0*q^2)) 

% Find the change in plasma beta to use in calculations.
% delta_beta = abs(plasma_beta_j - plasma_beta_g)

% Calculate magnetic shear angle between Jovian & Ganymede fields using dot
% product
% cos_theta_b = dot(b_dir_j/norm(b_dir_j),b_dir_g/norm(b_dir_j)) % returns
% cosine of the shear angle in radians
% magnetic_shear = acos(cos_theta_b)

%% Onset condition 1: Diamagnetic drift condition from Masters 2014:
% drift_condition = 2*atan(ion_inertial_length*delta_beta)/(2*ion_inertial_length))
% diamagnetic_drift = magnetic_shear > drift_condition     
% if diamagnetic_drift 
%     disp('Condition 1 met')
% else
%     disp('Condition 1 not met')
% end

%% Onset condition 2: flow shear onset condition from Masters 2014: