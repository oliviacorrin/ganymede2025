%% Graph spacecraft trajectory of G8 
clc
close all
clear all

% File name G8 with delimeterIn = , 
% headerlinesIn = 1 (start at 2nd column so we don't import the dates)
G8 = importdata('Galileo_O8.csv',',',2)

% create trajectory array with each column representing an x,y,z component of
% the spacecraft trajectory
trajectory = zeros(length(G8.data),3)
for i = 1:3
    trajectory(:,i) = G8.data(:,i+4)
end

% spacecraft-measured B field values
B_field = zeros(length(G8.data),3)
for i = 1:3
    B_field(:,i) = G8.data(:,i)
end

%% create 3d plot of trajectory 
plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3),'-b')
xlabel("x (in R_g)")
ylabel("y (in R_g)")
zlabel("z (in R_g)")
title("G8 trajectory practice graph")
hold on
% add Ganymede at origin 
r_factor = 0.038 % multiply by R_g/R_j
[X,Y,Z] = sphere
X2 = X*r_factor
Y2 = Y*r_factor
Z2 = Z*r_factor 
% assuming R_j 
% G1_plot = surf(X2,Y2,Z2)

% assuming R_g
G8_plot = surf(X,Y,Z)
G8_plot.FaceColor = [0.25 0.8 0.7]
legend("trajectory", "Ganymede")
axis equal

%% add quiver plot to show magnitude/direction of B field 
quiver3(trajectory(:,1),trajectory(:,2),trajectory(:,3),B_field(:,1),B_field(:,2),B_field(:,3))
hold off

%% plot B field time series
% plot the components (in GPhiO)
plot(B_field(:,1))
hold on
plot(B_field(:,2))
plot(B_field(:,3))
% add the magnitude 
plot(G8.data(:,4))
xlabel("time elapsed (# of data points out of total taken)")
ylabel("B field strength (nT)")
title("G8 flyby magnetometer measurements")

%% add apparent boundaries of magnetopause corresponding to sharp changes in
%% B field (all eyeballed by yours truly)
% first apparent boundary
xline(2284,'-k','LineWidth',2)
xline(4586,'-k','LineWidth',2)
% possible bounds to the entrance crossing 
xline(2132,'--','Color',[0 .8 .8],'LineWidth',2)
xline(2680,'--','Color',[0 .8 .8],'LineWidth',2)
% possible bounds to the exit crossing 
xline(4728,'--','Color',[0 .8 .8],'LineWidth',2)
xline(4159,'--','Color',[0 .8 .8],'LineWidth',2)
legend("B_x","B_y","B_z","B (mag)",'mpause entrance crossing','mpause exit crossing','possible min/max crossing bounds')



