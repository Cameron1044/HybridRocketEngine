clc; clear; close all;

m_loaded = 19.32933; % N2O mass initially loaded into tank [kg]: Test 1
Ti = 286.5; % initial temperature [K]: Test 1
V = 0.0354 ; % total tank volume [m^3]
R = 8314.3 ; % universal gas constant [J/(kmol*K)]
MW2 = 44.013; % molecular weight of N2O
G1 = 96.512 ; % vapor pressure of N2O [Pa] coefficients
G2 = -4045 ; % valid for Temp range [182.3 K - 309.57 K]
G3 = -12.277 ;
G4 = 2.886e-5 ;
G5 = 2 ; 
Q1 = 2.781; % molar specific volume of liquid N2O [m^3/kmol] coefficients
Q2 = 0.27244;
Q3 = 309.57;
Q4 = 0.2882;

n_to = m_loaded/MW2; % initial total N2O in tank [kmol]
Vhat_li = Q2^(1+(1-Ti/Q3)^Q4)/Q1; % molar volume of liquid N2O [m^3/kmol]
To = Ti ; % initial temperature [K]
P_sato = exp(G1 + G2/To + G3*log(To) + G4*To^G5); % initial vapor pressure of N20 [Pa]
n_go = P_sato*(V - Vhat_li*n_to) / (-P_sato*Vhat_li + R*To); % initial N2O gas [kmol]
n_lo = (n_to*R*To - P_sato*V) / (-P_sato*Vhat_li + R*To); % initial N2O liquid [kmol]

% Initial conditions
y0 = [Ti; n_go; n_lo];

% Time span
tspan = [0 10];
options = odeset('Events', @eventFunction);
% Solve the ODE system using ode45
[t, y] = ode45(@(t,y) odefun(t,y), tspan, y0, options);

% Extract the results
T = y(:, 1);
n_g = y(:, 2);
n_l = y(:, 3);

Pe = zeros(length(t), 1);
P = zeros(length(t), 1);
for i = 1:length(t)
    [~, ~, ~, Pe(i), P(i)] = tank_system(t(i), T(i), n_g(i), n_l(i));
end

% Plotting code
figure(1); grid on;
plot(t,T,'r','LineWidth',2)
title('Temperature vs. Time')
xlabel('Time [s]')
ylabel('Temperature [K]')

figure(2); grid on; hold on;
plot(t,n_g.*MW2,'b','LineWidth',2)
plot(t,n_l.*MW2,'g','LineWidth',2)
title('Mass of N20 vs. Time')
xlabel('Time [s]')
ylabel('Mass of N2O [kg]')
legend('Mass of N2O gas','Mass of N2O liquid')

figure(3); grid on; hold on;
plot(t,P,'m','LineWidth',2)
plot(t,Pe,'c','LineWidth',2)
title('Pressure vs. Time')
xlabel('Time [s]')
ylabel('Pressure [Pa]')
legend('tank pressure','chamber pressure')

function [value, isterminal, direction] = eventFunction(t, y)
    To = y(1);
    n_go = y(2);
    n_lo = y(3);
    [~, ~, ~, Pe, P] = tank_system(t, To, n_go, n_lo);

    % Initialize default values
    value = [Pe - P; n_lo];
    isterminal = [0; 0];
    direction = [0; 0]; % Detect zero crossing in any direction

    if Pe >= P
        isterminal(1) = 1; % Stop the integration if Pe >= P
    end

    if n_lo <= 0
        isterminal(2) = 1; % Stop the integration if n_lo <= 0
    end
end

function dy = odefun(t,y)
    To = y(1);
    n_go = y(2);
    n_lo = y(3);
    [dT, dn_g, dn_l, ~, ~] = tank_system(t, To, n_go, n_lo);
    dy = [dT; dn_g; dn_l];
end

function [dT, dn_g, dn_l, Pe, P] = tank_system(t, To, n_go, n_lo)
    Pe = -2924.42*t^6 + 46778.07*t^5 - 285170.63*t^4 + 813545.02*t^3 - 1050701.53*t^2 + 400465.85*t + 1175466.2; % Test 1

    % Given constants
    n_He = 0; % helium gas [kmol]
    Cd = 0.425; % discharge coefficient: Test 1
    Ainj = 0.0001219352; % injector area [m^2]
    V = 0.0354 ; % total tank volume [m^3]
    m_T = 6.4882; % tank mass [kg]
    R = 8314.3 ; % universal gas constant [J/(kmol*K)]
    MW2 = 44.013; % molecular weight of N2O
    G1 = 96.512 ; % vapor pressure of N2O [Pa] coefficients
    G2 = -4045 ; % valid for Temp range [182.3 K - 309.57 K]
    G3 = -12.277 ;
    G4 = 2.886e-5 ;
    G5 = 2 ; 
    Tc = 309.57 ; % critical temperature of N2O [K]
    J1 = 2.3215e7 ; % heat of vaporization of N2O [J/kmol] coefficients
    J2 = 0.384 ; % valid for Temp range [182.3 K - 309.57 K]
    J3 = 0 ;
    J4 = 0 ;
    C1 = 0.2079e5 ; % heat capacity of He at constant pressure [J/(kmol*K)] coefficients
    C2 = 0 ; % valid for Temp range [100 K - 1500 K]
    C3 = 0 ;
    C4 = 0 ;
    C5 = 0 ;
    D1 = 0.2934e5 ; % heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients
    D2 = 0.3236e5 ; % valid for Temp range [100 K - 1500 K]
    D3 = 1.1238e3 ;
    D4 = 0.2177e5 ;
    D5 = 479.4 ;
    E1 = 6.7556e4 ; % heat capacity of N2O liquid at constant pressure [J/(kmol*K)] coefficients
    E2 = 5.4373e1 ; % valid for Temp range [182.3 K - 200 K]
    E3 = 0 ;
    E4 = 0 ;
    E5 = 0 ;
    Q1 = 2.781; % molar specific volume of liquid N2O [m^3/kmol] coefficients
    Q2 = 0.27244;
    Q3 = 309.57;
    Q4 = 0.2882;

    % 
    Vhat_l = Q2^(1+(1-To/Q3)^Q4)/Q1;
    CVhat_He = C1 + C2*To + C3*To^2 + C4*To^3 + C5*To^4 - R;
    CVhat_g = D1 + D2*((D3/To)/sinh(D3/To))^2 + D4*((D5/To)/cosh(D5/To))^2 - R;
    CVhat_l = E1 + E2*To + E3*To^2 + E4*To^3 + E5*To^4;
    Tr = To/Tc;
    delta_Hv = J1*(1 - Tr) ^ (J2 + J3*Tr + J4*Tr^2);
    P_sat = exp(G1 + G2/To + G3*log(To) + G4*To^G5);
    dP_sat = (-G2/(To^2) + G3/To + G4*G5*To^(G5-1)) * exp(G1 + G2/To + G3*log(To) + G4*To^G5);
    Cp_T = (4.8 + 0.00322*To)*155.239;
    P = (n_He + n_go)*R*To / (V - n_lo*Vhat_l);

    %
    a = m_T*Cp_T + n_He*CVhat_He + n_go*CVhat_g + n_lo*CVhat_l; 
    b = P*Vhat_l;
    e = -delta_Hv + R*To; 
    f = -Cd*Ainj*sqrt(2/MW2)*sqrt((P-Pe)/Vhat_l);
    j = -Vhat_l*P_sat;
    k = (V - n_lo*Vhat_l)*dP_sat;
    m = R*To;
    q = R*n_go;
    
    %
    Z=(-f*(-j*a + (q-k)*b)) / (a*(m+j) + (q-k)*(e-b));
    W=(-Z*(m*a + (q-k)*e)) / (-j*a + (q-k)*b);
    
    dT = (b*W+e*Z)/a ;
    dn_g = Z;
    dn_l = W;
end




