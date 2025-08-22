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
plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3),'-k','LineWidth',2)
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
G8_plot.EdgeColor = 'none'
axis equal

%% add quiver plot to show magnitude/direction of B field 
quiver3(trajectory(:,1),trajectory(:,2),trajectory(:,3),B_field(:,1),B_field(:,2),B_field(:,3),'Color',[0.1 0.9 0])

%% add the magnetopause boundaries to more clearly visualize 
% define z values
z = -2:0.1:2

% create constant vectors for both crossing x and y values with the same number of elements as z
x_cross1 = -1.48 *ones(size(z))
y_cross1 = -1.22 * ones(size(z))

x_cross2 = -1.33*ones(size(z))
y_cross2 = 1.29*ones(size(z))

% plot lines
plot3(x_cross1, y_cross1, z, '--', 'Color',[.8 0 .2])
plot3(x_cross2, y_cross2, z, '--', 'Color',[.8 0 .2])

legend("trajectory", "Ganymede", "B field quiver plot", "magnetopause boundary")
grid on

hold off