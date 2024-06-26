load('resultHierarchical.mat');
load('resultSingleLayer.mat');
totalHours =  24;

figure(1);
subplot(2,4,1);
plot(0:1/12:totalHours-1/12, resH.x(1,:), 'LineWidth', 1.2); hold on;
plot(0:1/4:totalHours-1/4, resS.x(1,:), 'LineWidth', 1.2); title('Weight in kg m^{-2}');
subplot(2,4,2);
plot(0:1/12:totalHours-1/12, resH.x(2,:), 'LineWidth', 1.2); hold on;
plot(0:1/4:totalHours-1/4, resS.x(2,:), 'LineWidth', 1.2); title('Indoor CO2 in kg m^{-3}');%hold on;
% plot(  solH.value(xsphat(2,1))*ones(Nsim+1,1),'--', 'LineWidth', 1.2); title('Indoor CO2 in kg m^{-3}');
subplot(2,4,3);
plot(0:1/12:totalHours-1/12, resH.x(3,:), 'LineWidth', 1.2); hold on;
plot(0:1/4:totalHours-1/4, resS.x(3,:), 'LineWidth', 1.2);  title('Air Temp in  ^oC'); %hold on;
% plot( solH.value(xsphat(1,1))*ones(Nsim+1,1), '--', 'LineWidth', 1.2); title('Air Temp in  ^oC');
subplot(2,4,4);
plot(0:1/12:totalHours-1/12, resH.x(4,:), 'LineWidth', 1.2); hold on;
plot(0:1/4:totalHours-1/4, resS.x(4,:), 'LineWidth', 1.2);  title('RH C_{H2O} in %');
subplot(2,4,5);
plot(0:1/12:totalHours-1/12, resH.u(1,:), 'LineWidth', 1.2); hold on;
plot(0:1/4:totalHours-1/4, resS.u(1,:), 'LineWidth', 1.2); title('CO2 in mg*m^{-2}*s^{-1}');
subplot(2,4,6);
plot(0:1/12:totalHours-1/12, resH.u(2,:), 'LineWidth', 1.2); hold on;
plot(0:1/4:totalHours-1/4, resS.u(2,:), 'LineWidth', 1.2); title('vent. rate in mm s^{-1}');
subplot(2,4,7);
plot(0:1/12:totalHours-1/12, resH.u(3,:), 'LineWidth', 1.2); hold on;
plot(0:1/4:totalHours-1/4, resS.u(3,:), 'LineWidth', 1.2);  title(' heat supply in W*m^{-2}');
subplot(2,4,8);
plot(0:1/12:totalHours-1/12, resH.d(1,:), 'LineWidth', 1.2); hold on;
plot(0:1/4:totalHours-1/4, resS.d(1,:), 'LineWidth', 1.2);  title(' Incoming radiation in W m^{-2}');
legend('Hierarhical','SingleLayer');