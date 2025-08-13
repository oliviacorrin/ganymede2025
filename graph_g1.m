% Graph spacecraft trajectory of G1 

% File name G1 with delimeterIn = , 
% headerlinesIn = 1 (start at 2nd column so we don't import the dates)
G1 = importdata('Galileo_O1.csv',',',1)

% create trajectory array with each column representing an x,y,z component of
% the spacecraft trajectory
trajectory = zeros(length(G1.data),3)
for i = 1:3
    trajectory(:,i) = G1.data(:,i+4)
end

% create 3d plot of trajectory 
plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3),'-b')
xlabel("x (in R_j)")
ylabel("y (in R_j)")
zlabel("z (in R_j)")
title("G1 trajectory practice graph")
hold on
% add Ganymede at origin 
r_factor = 0.038 % multiply by R_g/R_j
[X,Y,Z] = sphere
X2 = X*r_factor
Y2 = Y*r_factor
Z2 = Z*r_factor 
% assuming R_j 
G1_plot = surf(X2,Y2,Z2)

% assuming R_g
% G1_plot = surf(X,Y,Z)
G1_plot.FaceColor = [0.25 0.8 0.7]
legend("trajectory", "Ganymede")
axis equal