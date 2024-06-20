figure(1);
subplot(2,4,1);hold on; %, 'Color', [0.4660 0.6740 0.1880]
plot( x_cl_buf(1,:), 'LineWidth', 1.2); title('Weight in kg m^{-2}');
subplot(2,4,2);hold on;
plot( x_cl_buf(2,:), 'LineWidth', 1.2); title('Indoor CO2 in kg m^{-3}');%hold on;
% plot(  solH.value(xsphat(2,1))*ones(Nsim+1,1),'--', 'LineWidth', 1.2); title('Indoor CO2 in kg m^{-3}');
subplot(2,4,3);hold on;
plot( x_cl_buf(3,:), 'LineWidth', 1.2);  title('Air Temp in  ^oC'); %hold on;
% plot( solH.value(xsphat(1,1))*ones(Nsim+1,1), '--', 'LineWidth', 1.2); title('Air Temp in  ^oC');
subplot(2,4,4);hold on;
plot( x_cl_buf(4,:), 'LineWidth', 1.2); title('RH C_{H2O} in %');
subplot(2,4,5);
plot( u_cl_buf(1,:), 'LineWidth', 1.2); title('CO2 in mg*m^{-2}*s^{-1}');
subplot(2,4,6);
plot( u_cl_buf(2,:), 'LineWidth', 1.2); title('vent. rate in mm s^{-1}');
subplot(2,4,7);
plot( u_cl_buf(3,:), 'LineWidth', 1.2); title(' heat supply in W*m^{-2}');
subplot(2,4,8);
plot( d_buf(1,:), 'LineWidth', 1.2); title(' Incoming radiation in W m^{-2}');