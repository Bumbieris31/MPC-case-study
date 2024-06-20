clear, clc
import casadi.*
% load parameters (p)
load('vanHentenParams.mat');
% load weather (weather)
load('weatherConditions.mat');

% x1: dry weight kg m^{-2}
% x2: indoor CO2 in kg m^{-3}
% x3: air temperature ^oC
% x4: humidity in kg m^{-3}

% y1: Weight in g m^{-2}
% y2: Indoor CO2 in ppm 10^{3}
% y3: Air Temp in  ^oC
% y4: RH C_{H2O} in %

% u1: Supply rate of CO2 in mg*m^{-2}*s^{-1}
% u2: Ventilation rate through the vents in mm s^{-1}
% u3: Energy supply by the heating system in W*m^{-2}

% d1: Incoming radiation in W m^{-2}
% d2: Outside CO2 in kg m^{-3}
% d3: Outdoor tempetarure in ^oC
% d4: Outdoor humidity content C_{H2O} in kg m^{-3}

% time step
Delta = 3*300; %  15 minutes
N = 20;

% initial conditions
t0 = 0;
x0 = [0.0035; 0.001; 15; 0.008];

%% LOW LEVEL EMPC  
% Downsample weather from 5 mins to 15 mins
d_cl = [weather.iGlob, weather.co2, weather.tOut, weather.hum];
d_cl = interp1(d_cl,1:(Delta/300):length(d_cl))';

% Define model from rk4
f = @(x,u,d) f_rk4(x,u,d,Delta,p);

% Define cost parameters
c_co2 = 0.42*Delta; % per kg/s of CO2
c_q = 6.35E-9*Delta; % per W of heat
c_dw = -16; % price per kg of dry weight

% Define optimization problem
opti = casadi.Opti();

% State and input trajectory
xhat = opti.variable(4,N+1);
uhat = opti.variable(3,N);
% Initial state
x0_p = opti.parameter(4);
% Weather in trajectory
d_p = opti.parameter(4,N);

cost = 0;
opti.subject_to(xhat(:,1)==x0_p);
for k = 1:N
    % Enforce model
    opti.subject_to(xhat(:,k+1)==f(xhat(:,k),uhat(:,k),d_p(:,k)))
    
    % Input constraints
    opti.subject_to(0 <= uhat(1,k) <= 1.2);
    opti.subject_to(0 <= uhat(2,k) <= 7.5);
    opti.subject_to(0 <= uhat(3,k) <= 150);
    
    % State/output constraints
    yhat = h_meas(xhat(:,k),p);
    opti.subject_to(0 <= yhat(1));
    opti.subject_to(0 <= yhat(2) <= 1.6);
    opti.subject_to(10 <= yhat(3) <= 25);
    opti.subject_to(0 <= yhat(4) <= 80);
    
    % Add stage cost
    cost = cost + c_dw*(xhat(1,k+1)-xhat(1,k)) ...
        + c_co2*uhat(1,k)*10^(-6) ...
        + c_q*uhat(3,k);
end
% Define objective function
opti.minimize(cost)

% Solver options
solver_options = struct();
% Suppress output
solver_options.ipopt.print_level = 0;
solver_options.print_time = 0;
solver_options.ipopt.sb = 'yes';
% Select solver
opti.solver('ipopt',solver_options);

% Closed-loop simulation
Nsim = 4*24; % Number of simulation time steps

% Initialize closed-loop state and input variables
x_cl = zeros(4,Nsim+1);
y_cl = zeros(4,Nsim+1);
u_cl = zeros(3,Nsim);

% Set initial condition
x_cl(:,1) = x0;
y_cl(:,1) = h_meas(x0,p);
disp("Simulation step:")
total_cost = 0;

u_cl_buf = [];
x_cl_buf = [];
y_cl_buf = [];
d_buf    = [];

for k=1:Nsim
    disp(k)
    % Set current state at time k as initial state in MPC optimization
    opti.set_value(x0_p,x_cl(:,k));
    opti.set_value(d_p,d_cl(:,t0+k:(t0+k+N-1)));
    % Solve optimization problem
    sol = opti.solve();
    % Save warm-start
    opti.set_initial(sol.value_variables());
    % Get optimal input trajectory
    u_opt = sol.value(uhat);
    % Take only the *first* input
    u_cl(:,k) = u_opt(:,1);
    % Simulate system
    x_cl(:,k+1) = f(x_cl(:,k),u_cl(:,k),d_cl(:,t0+k));
    y_cl(:,k+1) = h_meas(x_cl(:,k+1),p);
    
    total_cost = total_cost + c_dw*(x_cl(1,k+1)-x_cl(1,k)) ...
        + c_co2*u_cl(1,k)*10^(-6) ...
        + c_q*u_cl(3,k);

 % Store 
    u_cl_buf = [u_cl_buf u_cl(:,k)];
    x_cl_buf = [x_cl_buf x_cl(:,k+1)];
    y_cl_buf = [y_cl_buf y_cl(:,k+1)];
    d_buf    = [d_buf d_cl(:,t0+k)];
end

%%
Plotting

resS.x = x_cl_buf;
resS.y = y_cl_buf;
resS.u = u_cl_buf;
resS.d = d_buf;
save('resultSingleLayer.mat','resS');
% figure(1);
% subplot(1,4,1);hold on; %, 'Color', [0.4660 0.6740 0.1880]
% plot(0:Nsim, y_cl(1,:), 'LineWidth', 1.2); title('Weight in g m^{-2}');
% subplot(1,4,2);hold on;
% plot(0:Nsim, y_cl(2,:), 'LineWidth', 1.2); title('Indoor CO2 in ppm 10^{3}');
% subplot(1,4,3);hold on;
% plot(0:Nsim, y_cl(3,:), 'LineWidth', 1.2); title('Air Temp in  ^oC');
% subplot(1,4,4);hold on;
% plot(0:Nsim, y_cl(4,:), 'LineWidth', 1.2); title('RH C_{H2O} in %');
%%
% Function to perform explicit rk4
function xplus = f_rk4(x,u,d,Delta,p)
n_rk4 = 1;
delta_rk4 = Delta/n_rk4;
for i=1:n_rk4
    k_1 = vanHentenODEs(x, u, d, p);
    k_2 = vanHentenODEs(x+0.5*delta_rk4*k_1, u, d, p);
    k_3 = vanHentenODEs(x+0.5*delta_rk4*k_2, u, d, p);
    k_4 = vanHentenODEs(x+delta_rk4*k_3, u, d, p);
    x = x + (1/6)*(k_1+2*k_2+2*k_3+k_4)*delta_rk4;
end
xplus = x;
end

% Function to convert states to outputs
function y = h_meas(x,p)
y_1 = 10^(3)*x(1);                                                    % Weight in g m^{-2}
y_2 = p(12)*(x(3)+p(13))/(p(14)*p(15)).*x(2)*10^(3);                 % Indoor CO2 in ppm 10^{3}
y_3 = x(3);                                                               % Air Temp in  ^oC
y_4 = p(12)*(x(3) + p(13))./(11*exp(p(27)*x(3)./(x(3)+p(28)))).*x(4)*10^(2); % RH C_{H2O} in %
y = [y_1, y_2, y_3, y_4];
end
