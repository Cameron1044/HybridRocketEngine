clc; clear; close all;

%% Inputs
% Conversion Factors
psi2Pa = 6894.757293;       % [Pa/psi]
in2m = 0.0254;              % [m/in]
lbm2kg = 0.45359237;        % [kg/lbm]
lbf2N = 4.448221615;        % [N/lbf]

% Constants
Ru = 8.3145;                % [J/(mol*K)]
P_amb = 80.8E3;             % [Pa]

% Oxidizer Tank
V_t = 3E-3;                 % [m^3] Tank Volume
P_o_i = 3000 * psi2Pa;      % [Pa] Initial Tank Pressure
rho_o = 1226;               % [kg/m^3] Liquid Nitrous Density @ 273 K

% Injector
% n_i = 4;                    % Number of Injector Holes
% d_i = 1/8 * in2m;           % [m] Diameter of Injector Holes
% A_i = n_i * pi * (d_i/2)^2; % [m^2] Area of Injector
A_i = 1.02E-05;
C_d = 0.4;                  % Coefficient of Discharge Injector

% Fuel Grain
L = 12 * in2m;              % [m]
fuel_OD = 1.688 * in2m;     % [m]
fuel_ID = 1.35 * in2m;      % [m]
rho_f = 2767.99;            % [kg/m^3]
    % Combustion Chamber Gas Chemical Constants
T_c = 5583.19;              % [K]
M_c = 39.992 / 1000;        % [g/mol --> kg/mol] Effective Molecular Weight
gam = 1.1587;               % [1]
    % St. Roberts Law
n = 0.8;                    % [1]   Fuel Regression Constant 1
a_m = 0.0007;               %  Fuel Regression Constant 2
convert_a = (0.0254^(1+2*(n)))*(0.453592^(-n)); % [m^(1+2n)kg^(-n)s^(n-1)]
a = convert_a * a_m;

P_c_i = P_amb;              % [Pa] Initial Chamber Pressure

% Nozzle
d_t = 0.5 * in2m;           % [m] Diameter of Nozzle throat
A_t = pi*(d_t/2)^2;         % [m^2] Area of Nozzle throat

%% Initializing Calculations
% Calculating Oxidizer Tank Ullage/Fuel Ratios
u_frac = 0.2;               % [1] Assuming Ullage is 20% of tank
V_u_i = V_t * u_frac;       % [m^3] Ullage Volume
V_o_i = V_t * (1 - u_frac); % [m^3] Oxidizer Volume
m_o_i = V_o_i * rho_o;      % [kg] Mass of Liquid Oxidizer

% Calculating Fuel Grain Volume
V_c_tot = pi*(fuel_OD/2)^2*L;   % [m^3] Outer Diameter
V_c_i = pi*(fuel_ID/2)^2*L;     % [m^3] Inner Diameter
V_f_i = V_c_tot - V_c_i;        % [m^3] Outer - Inner Final Volume
m_f_i = rho_f*pi*(fuel_OD^2-fuel_ID^2)/4; % [kg] Fuel Grain Mass
port_i = fuel_ID/2;             % [m] Initial Port radius
R = Ru/M_c;                     % Combustion Gas constant [m^2 / (K * s^2)]

tspan = [0 8];
var_i = [P_c_i; P_o_i; m_f_i; m_o_i; port_i; 0];

%% Run ode45

opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

[t, var] = ode45(@(t, var) NakkaHybrid(t, var, Ru, gam, T_c, M_c, L, rho_f, rho_o, a, n, A_i, C_d, A_t, P_o_i, V_u_i, V_t, P_amb, psi2Pa, fuel_OD), tspan, var_i,opts);

%% Recover Data
P_c = var(:,1);
P_o = var(:,2);
m_f = var(:,3);
m_o = var(:,4);
port = var(:,5);
I = var(:,6);

der = zeros(length(t),6);

for i = 1:length(t)
    der(i,:) = NakkaHybrid(t(i), var(i,:), Ru, gam, T_c, M_c, L, rho_f, rho_o, a_m, n, A_i, C_d, A_t, P_o_i, V_u_i, V_t, P_amb, psi2Pa, fuel_OD);
end

dPc_dt = der(:,1);
dPo_dt = der(:,2);
mdot_f = der(:,3);
mdot_o = der(:,4);
F = der(:,6);

%% Plot Data
figure()
grid minor
hold on
plot(t, P_c/psi2Pa, 'LineWidth', 2)
plot(t, P_o/psi2Pa, 'LineWidth', 2)
title("Rocket Pressures Over Time")
legend("Chamber Pressure", "Oxidizer Tank Pressure")
legend("Location", "best")
xlabel("Time (s)")
ylabel("Pressure (psi)")

figure()
grid minor
hold on
% Changed for intuitively positive mass flow rate
plot(t, mdot_f/lbm2kg * -1, 'LineWidth', 2)
plot(t, mdot_o/lbm2kg * -1, 'LineWidth', 2)
title("Rocket Mass Flow Rate over Time")
legend("Fuel Mass Flow", "Oxidizer Mass Flow")
legend("Location", "best")
xlabel("Time (s)")
ylabel("Mass Flow Rate (lbm/s)")

figure()
grid minor
hold on
plot(t, port/in2m, 'LineWidth', 2)
yline(fuel_OD/2/in2m, '--r', 'LineWidth', 2)
legend("Fuel Grain Port Radius", "Fuel Grain Outer Diameter")
legend("Location", "best")
title("Port Radius vs Time")
xlabel("Time (s)")
ylabel("Port Radius (in)")

figure()
subplot(1,2,1)
grid minor
hold on
plot(t, F/lbf2N, 'LineWidth', 2)
title("Thrust Curve")
xlabel("Time (s)")
ylabel("Thrust (lbf)")
hold off
subplot(1,2,2)
grid minor
hold on
plot(t, I/lbf2N, 'LineWidth', 2)
title("Impulse Curve")
xlabel("Time (s)")
ylabel("Total Impulse Produced (lbf-s)")
hold off

figure()
grid on
grid minor
hold on
plot(t, mdot_o./mdot_f)
title("Mixture Ratio vs Time")
xlabel("Time (s)")
ylabel("O/F Ratio")

I_tot = I(end)/lbf2N % [lbf-s]
I_sp = I_tot/((m_f_i+m_o_i)/lbm2kg)

%% Function
    % ODE45 function used to talk in design inputs and output system
    % variables and especially thrust at any time step.
function var_dot = NakkaHybrid(t, var, Ru, gam, T_c, M_c, L, rho_f, rho_o, a, n, A_i, C_d, A_t, P_o_i, V_u_i, V_t, P_amb, psi2Pa, fuel_OD)
    P_c = var(1);       % Chamber Pressure
    P_o = var(2);       % Oxidizer Tank Pressure
    m_f = var(3);       % Fuel grain mass
    m_o = var(4);       % Liquid Oxidizer mass
    port = var(5);      % Port radius
    
    R = Ru/M_c; 
    A_b = 2*pi*port*L;      % [m^2] Burning Surface Area 
    V_c = pi*port.^2*L;     % [m^3] Chamber Volume
    rho_c = P_c/(R*T_c);    % [kg/m^3] Chamber Density
    P_e = P_amb;% I question this but if you're going for a totally optimal design then sure
        % Perfect expansion will address later

    % Events at Each Time step
    if port >= fuel_OD/2
        % If the port radius burns more than the Outer Fuel Grain Radius
        A_b = 0;
        a = 0;
    end

    if m_o <= 0
        % When Liquid Oxidizer mass runs out
        A_i = 0;
    end

    if P_c <= P_amb
        % When the Chamber pressure is less than ambient pressure
        A_t = 0;
        P_c = P_amb;
    end
    
    % Burn Rate 
    r = a*(P_c/psi2Pa)^n;

    % Fuel mass flow rate
    mdot_f = -A_b*rho_f*r; 
    
    % Oxidizer mass flow rate
    mdot_o = -A_i*C_d*sqrt(2*rho_o*(P_o-P_c));

    % Tank pressure derivative
    dPo_dt = -P_o_i*V_u_i/(V_t - m_o/rho_o)^2 * A_i*C_d*sqrt(2*(P_o-P_c)/rho_o); 

    % Chamber pressure derivative
    dPc_dt = R*T_c./V_c .* (A_b*r*(rho_f-rho_c) - mdot_o - P_c*A_t*sqrt(gam/(R*T_c))*(2/(gam+1))^((gam+1)/(2*(gam-1))));   

    % Thrust Calculation
    F = A_t*P_c*sqrt(2*gam^2/(gam-1) * (2/(gam+1))^((gam+1)/(gam-1)) * (1-(P_e/P_c)^((gam-1)/gam))); % + (P_e-P_amb)*A_e;
    
    %t
    var_dot = [dPc_dt; dPo_dt; mdot_f; mdot_o; r; F];
end
