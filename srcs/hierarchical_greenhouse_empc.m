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

    % initial conditions
t0 = 0;
x0 = [0.0035; 0.001; 15; 0.008];

u_cl_buf = [];
x_cl_buf = [];
y_cl_buf = [];
d_buf = [];

for i=1:30*24
    %% HIGH LEVEL EMPC
    % time step
    Delta = 12*300; %  15 minutes
    N = 24*7;
    % Downsample weather from 5 mins to 1hour
    d_cl = [weather.iGlob, weather.co2, weather.tOut, weather.hum];
    d_cl = interp1(d_cl,1:(Delta/300):length(d_cl))';
    
    % Define model from rk4
    f = @(x,xsp,d) cfruit_rk4(x,xsp,d,Delta,p);
    
    % Define cost parameters
    c_co2 = 0.42*Delta; % per kg/s of CO2
    c_q = 6.35E-9*Delta; % per W of heat
    c_dw = -16; % price per kg of dry weight
    c_tair = 1; % penalty for tAir tracking
    c_tco2 = 5e3; % penalty for co2 tracking
    
    % Define optimization problem
    opti = casadi.Opti();
    
    % State and input trajectory
    xhat = opti.variable(1,N+1);
    xsphat = opti.variable(2,N);
    % Initial state
    x0_p = opti.parameter(1);
    % Weather in trajectory
    d_p = opti.parameter(4,N);
    
    cost = 0;
    opti.subject_to(xhat(:,1)==x0_p);
    for k = 1:N
        % Enforce model
        opti.subject_to(xhat(:,k+1)==f(xhat(:,k),xsphat(:,k),d_p(:,k)))
        
        % Input constraints
        iglob = d_cl(1,i:(i+N-1));
        if iglob(k) > 10
            d = 1;
        else 
            d = 0;
        end
        opti.subject_to(10 +5*d <= xsphat(1,k) <= 30);
        opti.subject_to(0 <= xsphat(2,k));% <= 7.5);
        
        % State/output constraints
        yhat = 10^(3)*xhat(:,k);
        opti.subject_to(0 <= yhat(1));
        % opti.subject_to(0 <= yhat(2) <= 1.6);
        % opti.subject_to(10 <= yhat(3) <= 25);
        % opti.subject_to(0 <= yhat(4) <= 80);
        
        % Add stage cost
        cost = cost + c_dw*(xhat(1,k+1)-xhat(1,k)) ... %[0, 48E-5]
             + xsphat(2,k)*10^(-3);% ... % [0, 2.4E-5]
             %+ 10^(-6)*xsphat(1,k);
    
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
    
    % Set current state at time k as initial state in MPC optimization
    opti.set_value(x0_p,x0(1,1));
    opti.set_value(d_p,d_cl(:,i:(i+N-1)));
    % Solve optimization problem
    solH = opti.solve();   
    
    
    %% LOW LEVEL EMPC  
    % time step
    Delta = 300; %  5 minutes
    N = 12;
    solH.value(xsphat(1,1))
    xsp = [zeros(N,1), solH.value(xsphat(2,1))*ones(N,1), solH.value(xsphat(1,1))*ones(N,1), zeros(N,1)]';
    % Downsample weather from 5 mins to 15 mins
    d_cl = [weather.iGlob, weather.co2, weather.tOut, weather.hum];
    d_cl = interp1(d_cl,1:(Delta/300):length(d_cl))';
    
    % Define model from rk4
    f = @(x,u,d) f_rk4(x,u,d,Delta,p);
    
    % Define cost parameters
    c_co2 = 0.42*Delta; % per kg/s of CO2
    c_q = 6.35E-9*Delta; % per W of heat
    c_dw = -16; % price per kg of dry weight
    c_tair = 1; % penalty for tAir tracking
    c_tco2 = 5e3; % penalty for co2 tracking
    
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
        % cost = cost + c_dw*(xhat(1,k+1)-xhat(1,k)) ...
        %     + c_co2*uhat(1,k)*10^(-6) ...
        %     + c_q*uhat(3,k);
        cost = cost + c_tco2*(xsp(2,k) - xhat(2,k))^2 + ...
            c_tair*(xsp(3,k) - xhat(3,k))^2 + ...
            c_co2*uhat(1,k)*10^(-6) + ...
            uhat(2,k) + ...
            c_q*uhat(3,k);
    
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
    Nsim = 12;%4*24; % Number of simulation time steps
    
    % Initialize closed-loop state and input variables
    x_cl = zeros(4,Nsim+1);
    y_cl = zeros(4,Nsim+1);
    u_cl = zeros(3,Nsim);
    
    % Set initial condition
    x_cl(:,1) = x0;
    y_cl(:,1) = h_meas(x0,p);
    disp("Simulation step:")
    total_cost = 0;
    for k=1:Nsim
        % disp(k)
        % Set current state at time k as initial state in MPC optimization
        opti.set_value(x0_p,x_cl(:,k));
        opti.set_value(d_p,d_cl(:,t0+k:(t0+k+N-1)));
        % Solve optimization problem
        sol = opti.solve(); %sol.value(xsphat);
        % Save warm-start
        opti.set_initial(sol.value_variables());
        % Get optimal input trajectory
        u_opt = sol.value(uhat);
        % Take only the *first* input
        u_cl(:,k) = u_opt(:,1);
        % Simulate system
        x_cl(:,k+1) = f(x_cl(:,k),u_cl(:,k),d_cl(:,t0+k));
        y_cl(:,k+1) = h_meas(x_cl(:,k+1),p);

        % Store 
        u_cl_buf = [u_cl_buf u_cl(:,k)];
        x_cl_buf = [x_cl_buf x_cl(:,k+1)];
        y_cl_buf = [y_cl_buf y_cl(:,k+1)];
        d_buf    = [d_buf d_cl(:,t0+k)];
        
        total_cost = total_cost + c_dw*(x_cl(1,k+1)-x_cl(1,k)) ...
            + c_co2*u_cl(1,k)*10^(-6) ...
            + c_q*u_cl(3,k);
    end
    t0 = i*12;
    x0 = x_cl(:,k+1); %[0.0035; 0.001; 15; 0.008];
end
Plotting

resH.x = x_cl_buf;
resH.y = y_cl_buf;
resH.u = u_cl_buf;
resH.d = d_buf;
save('resultHierarchical.mat','resH');
%%

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

function xplus = cfruit_rk4(x,xsp,d,Delta,p)
n_rk4 = 1;
delta_rk4 = Delta/n_rk4;
for i=1:n_rk4
    k_1 = cfruitODE(x, xsp, d, p);
    k_2 = cfruitODE(x+0.5*delta_rk4*k_1, xsp, d, p);
    k_3 = cfruitODE(x+0.5*delta_rk4*k_2, xsp, d, p);
    k_4 = cfruitODE(x+delta_rk4*k_3, xsp, d, p);
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

function dx1dt = cfruitODE(x, xsp, d, p)
p1 = p(1); p2 = p(2); p3 = p(3); p4 = p(4); p5 = p(5); p6 = p(6);
p7 = p(7); p8 = p(8); p9 = p(9); p10 = p(10); p11 = p(11); p12 = p(12);
p13 = p(13); p14 = p(14); p15 = p(15); p16 = p(16); p17 = p(17);
p18 = p(18); p19 = p(19); p20 = p(20); p21 = p(21); p22 = p(22);
p23 = p(23); p24 = p(24); p25 = p(25); p26 = p(26); p27 = p(27);
p28 = p(28);

w1 = 0; w2 = 0; w3 = 0; w4 = 0; % No disturbances

x1 = x(1); %x2 = x(2); x3 = x(3); x4 = x(4);

xsp1        = xsp(1); % Tair setpoint {^oC}
xsp2        = xsp(2); % co2 setpoint in kg m^{-3}

d1        = d(1); % Incoming radiation in W m^{-2}

phi       = p4*d1 + (-p5*xsp1.^2 + p6*xsp1 - p7)*(xsp2 - p8);
PhiPhot_c = (1-exp(-p3*x1))*(p4*d1*(-p5*xsp1.^2 + p6*xsp1 - p7)*(xsp2-p8))/phi;         % gross canopy phootsynthesis rate

dx1dt = (p1*PhiPhot_c - p2*x1*2^(xsp1/10 - 5/2))*(1+w1);
end
