% Masters Thesis *PAGE 127*
% Ideal Model
% Explicit Method - Forward Difference
clc; clear; close all;

% Given constants
n_He = 0; % helium gas [kmol]
m_loaded = 19.32933; % N2O mass initially loaded into tank [kg]: Test 1
% m_loaded = 16.23298; % Test 2
% m_loaded = 14.10076; % Test 3
% m_loaded = 23.62427; % Test 4
Ti = 286.5; % initial temperature [K]: Test 1
% Ti = 278.5; % Test 2
% Ti = 271.5; % Test 3
% Ti = 291.3; % Test 4
R = 8314.3 ; % universal gas constant [J/(kmol*K)]
Cd = 0.425; % discharge coefficient: Test 1
% Cd = 0.365; % Test 2 and 3
% Cd = 0.09; % Test 4
Ainj = 0.0001219352; % injector area [m^2]
MW2 = 44.013; % molecular weight of N2O
V = 0.0354 ; % total tank volume [m^3]
m_T = 6.4882; % tank mass [kg]

% Perry's Chemical Engineers' Handbook Property Equations
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

% Initial Conditions
n_to = m_loaded/MW2; % initial total N2O in tank [kmol]
Vhat_li = Q2^(1+(1-Ti/Q3)^Q4)/Q1; % molar volume of liquid N2O [m^3/kmol]
To = Ti ; % initial temperature [K]
P_sato = exp(G1 + G2/To + G3*log(To) + G4*To^G5); % initial vapor pressure of N20 [Pa]
n_go = P_sato*(V - Vhat_li*n_to) / (-P_sato*Vhat_li + R*To); % initial N2O gas [kmol]
n_lo = (n_to*R*To - P_sato*V) / (-P_sato*Vhat_li + R*To); % initial N2O liquid [kmol]

% Forward Difference Time Loop
tf=100; % final time [s]
tstep=0.0005; % time step [s]
i_i=0;
i_f=tf/tstep;

for i=i_i:i_f
    t = i*tstep;
    % Curve fitted combustion chamber pressure [Pa]:
    Pe = -2924.42*t^6 + 46778.07*t^5 - 285170.63*t^4 + 813545.02*t^3 - 1050701.53*t^2 + 400465.85*t + 1175466.2; % Test 1
    % Pe = 95.92*t^6 - 2346.64*t^5 + 21128.78*t^4 - 87282.73*t^3 + ...
    % 186675.17*t^2 - 335818.91*t + 3029190.03; % Test 2
    % Pe = 58.06*t^6 - 1201.90*t^5 + 8432.11*t^4 - 22175.67*t^3 + ...
    % 21774.66*t^2 - 99922.82*t + 2491369.68; % Test 3
    % Pe = -4963.73*t + 910676.22; % Test 4

    % Given functions of temperature:
    Vhat_l = Q2^(1+(1-To/Q3)^Q4)/Q1;
    %molar specific volume of liquid N2O [m^3/kmol]
    CVhat_He = C1 + C2*To + C3*To^2 + C4*To^3 + C5*To^4 - R;
    %specific heat of He at constant volume [J/(kmol*K)]
    CVhat_g = D1 + D2*((D3/To)/sinh(D3/To))^2 + D4*((D5/To)/cosh(D5/To))^2 - R;
    %specific heat of N2O gas at constant volume [J/(kmol*K)]
    CVhat_l = E1 + E2*To + E3*To^2 + E4*To^3 + E5*To^4;
    %specific heat of N2O liquid at constant volume, approx. same as at constant pressure [J/(kmol*K)]
    Tr = To/Tc; % reduced temperature
    delta_Hv = J1*(1 - Tr) ^ (J2 + J3*Tr + J4*Tr^2); % heat of vaporization of N2O [J/kmol]
    P_sat = exp(G1 + G2/To + G3*log(To) + G4*To^G5); % vapor pressure of N20 [Pa]
    dP_sat = (-G2/(To^2) + G3/To + G4*G5*To^(G5-1)) * exp(G1 + G2/To + G3*log(To) + G4*To^G5);
    %derivative of vapor pressure with respect to temperature
    Cp_T = (4.8 + 0.00322*To)*155.239; % specific heat of tank, Aluminum [J/(kg*K)]

    % Simplified expression definitions for solution
    P = (n_He + n_go)*R*To / (V - n_lo*Vhat_l) ;
    a = m_T*Cp_T + n_He*CVhat_He + n_go*CVhat_g + n_lo*CVhat_l;
    b = P*Vhat_l;
    e = -delta_Hv + R*To;
    f = -Cd*Ainj*sqrt(2/MW2)*sqrt((P-Pe)/Vhat_l);
    j = -Vhat_l*P_sat;
    k = (V - n_lo*Vhat_l)*dP_sat;
    m = R*To;
    q = R*n_go;

    Z=(-f*(-j*a + (q-k)*b)) / (a*(m+j) + (q-k)*(e-b));
    W=(-Z*(m*a + (q-k)*e)) / (-j*a + (q-k)*b);

    % Derivative Functions
    dT = (b*W+e*Z)/a ;
    dn_g = Z;
    dn_l = W;

    % Record variables for each time step in an array
    T(i+1,1) = t;
    T(i+1,2) = To;
    n_g(i+1,1) = t;
    n_g(i+1,2) = n_go;
    n_l(i+1,1) = t;
    n_l(i+1,2) = n_lo;
    Pres(i+1,1) = t;
    Pres(i+1,2) = P;
    PE(i+1,1) = t;
    PE(i+1,2) = Pe;

    % Forward Difference Method
    To = To + dT*tstep;
    n_go = n_go + dn_g*tstep;
    n_lo = n_lo + dn_l*tstep;

    % Physical stops to kick out of loop
    if Pe>=P
        break
    end
    if n_lo<=0
        break
    end
end

% Plot results
figure(1); grid on;
plot(T(:,1),T(:,2),'r','LineWidth',2)
title('Temperature vs. Time')
xlabel('Time [s]')
ylabel('Temperature [K]')

figure(2); grid on; 
plot(n_g(:,1),n_g(:,2),'b',n_l(:,1),n_l(:,2),'g','LineWidth',2)
title('kmol of N20 vs. Time')
xlabel('Time [s]')
ylabel('kmol of N2O [kmol]')
legend('kmol of N2O gas','kmol of N2O liquid')

figure(3); grid on;
plot(Pres(:,1),Pres(:,2),'m',PE(:,1),PE(:,2),'c','LineWidth',2)
title('Pressure vs. Time')
xlabel('Time [s]')
ylabel('Pressure [Pa]')
legend('tank pressure','chamber pressure')