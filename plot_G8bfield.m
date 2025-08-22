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