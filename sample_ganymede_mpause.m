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
v_j = 400 *10^3; % (m/s) Jovian plasma velocity
rho_j = [] % (kg/m^3) Jovian plasma density - number density * pmass 
b_j = [] % (T) Jovian magnetic field strength
b_dir_j = [] % Jovian magnetic field direction (in GPhiO?)

% define moon parameters
% right now these are constant; later some will be varied to match up with
% different orientations of the moon in the PS.
r_g = 2.634*10^6 % (m) Ganymede radius
b_eq_g = 719*10^-9 % (T) Ganymede equatorial B field strength 
b_dir_g = [1 1 1] % Ganymede B field direction vector
rho_g = % (kg/m^3) Ganymede plasma density 
plasma_beta_g = % plasma beta value at Ganymede

% This will eventually be calculated in the code for varying upstream
% parameters(?), but for now, using a constant
mp = 2.2 % (in r_g) magnetopause standoff distance IN plasma sheet using jovianden = 56 amu/cm^3

%% Volume 
% Currently using a plasma density value of 56 amu/cm^3 corresponding to a
% position inside the Jovian plasma sheet. Eventually plan to add other
% values corresponding to different locations wrt the PS.

% get 2D (y,z) GPhiO coordinates and use to create paraboloid meshgrid
% goal: to extract these from data (using rough estimates for now)
n = 101
x = linspace(-mp-1,10,n)
y = linspace(-10,10,n)
z = linspace(-10,10,n)
[X1,Y1,Z1] = meshgrid(x,y,z)

% paraboloid magnetopause surface - from Donaldson et al. 2024, derived
% from Cooling et al. 2001:
% Y^2 + Z^2 = 2mp(X+mp)
% rearrange to find X:
% (Y^2 + Z^2 - mp)/2mp = X
[Y2,Z2] = meshgrid(x,y,z)
X2 = (Y2.^2 + Z2.^2 - mp)/2*mp

% paraboloids on either side of mpause 
mp1 = mp - 0.3
mp2 = mp + 0.3
p1 = (Y2.^2 + Z2.^2 - mp1)/2*mp1
p2 = (Y2.^2 + Z2.^2 - mp2)/2*mp2

% plot using surf(x,y,z) function, which represents each value in matrix Z as the height of a surface above x-y plane.
figure
paraboloid = surf(X,Y,Z)
title('Paraboloid model of Ganymede magnetopause with standoff distance 2.2 R_G')
xlabel('X (R_G')
ylabel('Y (R_G)')
zlabel('Z (R_G)')
grid on

% find normal vector to the paraboloid
[nx,ny,nz] = surfnorm(paraboloid,Y,Z)

% Now check each conditional to determine if reconnection is possible or
% prohibited.
% Onset condition 1: diamagnetic drift condition. 

% Onset condition 2: flow shear onset condition. 