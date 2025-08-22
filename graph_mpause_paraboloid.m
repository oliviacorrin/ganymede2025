%% Create magnetopause paraboloid graph.
% This is a very basic version of the paraboloid model, found using 3
% points: 
% 1: standoff distance calculated via pressure balance (-2.22, 0) 
% 2 & 3: found using G8 flyby data - matched B field fluctuations to
% spatial data and found crossings estimated at (-1.48, -1.22) and (-1.33,
% 1.29)

% create y & x vectors
y = linspace(-5,5,100)
x = 0.504*y.^2 + 0.0245.*y - 2.22

% create plot 
plot(x,y)
grid on
xlim([-3 1.5])