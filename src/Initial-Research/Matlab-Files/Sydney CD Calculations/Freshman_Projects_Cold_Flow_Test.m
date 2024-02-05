clear; clc; close all

data = importdata("Freshmen_cold flow data.txt");

time_bad = data((9910:10546),1)/1000; % sec
weight_bad = data((9910:10546),2) / 2.205; % kg
pressure_bad = data((9910:10546),3) * 6895; % pa

g = 9.81;

N = 130;

time = time_bad(1:N,1);
mass = weight_bad(1:N,1)/g;
pressure = pressure_bad(1:N,1);

%% Getting initial and Final Mass
mass_init = abs(mass(1,1));
mass_end = abs(mass(end,1));

%% Calculate Mass Flow Rate
delta_time = diff(time);
delta_mass = diff(mass);

mass_flow_rate = (mass_end - mass_init) ./ delta_time;

p = polyfit(time(1:end-1,1),mass(1:end-1,1),2);
val = polyval(p,time(1:end-1,1));

dxval = p(1,1)*2.*time(1:end-1,1) + p(1,2);

%% Theoretical Mass Flow Rate
num_holes = 5;
diameter = 0.03 / 39.37; % Diameter of one hole in m
radius = diameter / 2; % Radius of one hole in m
area_one_hole = pi * (radius)^2; % Area of one hole in square m
orifice_area = num_holes * area_one_hole; % Total orifice area in square m

rho = 0.04429207176 * 27680; % kg/m3
atmospheric_pressure = 14.93 * 6894.76; % pa

theoretical_mass_flow_rate = sqrt(2 * rho .* (pressure(1:end-1) - atmospheric_pressure)) * orifice_area;

%% Estimate of Coefficient of Discharge
coefficient_of_discharge = dxval(:,1) ./ theoretical_mass_flow_rate(:,1);

cd_est = (-dxval(end,1)+dxval(1,1)) ./ (theoretical_mass_flow_rate(end,1)-theoretical_mass_flow_rate(1,1));
% ^ This one gives 3x the cd_estimate form just using the mean of the
% equation so idk if that is better or worse for this context.

cd_estimate = mean((coefficient_of_discharge));

test_cd = mass_flow_rate ./ theoretical_mass_flow_rate;
test_mean = mean(test_cd);

% g = polyfit(time(1:end-1,1),test_cd(:,1),1);
% valg = polyval(g,time(1:end-1,1));

% gval = g(1,1)*time(1:end-1,1) + g(1,2);

% coeff_of_discharge_test = gval(:,1) ./ theoretical_mass_flow_rate(:,1);
% cd_estimate_test = mean((coeff_of_discharge_test));

disp('Cd estimate:');
disp(cd_estimate);

% disp('Cd estimate:');
% disp(cd_estimate_test);

%% Analyze and Validate Data
figure(); hold on
grid minor
plot(time(1:end-1,1), coefficient_of_discharge(:,1), "Color","m");
%plot(time(1:end-1,1), coeff_of_discharge_test(:,1), "Color","g");
%plot(time(1:end-1,1), test_cd(:,1), "Color","r")
plot(time(1:end-1,1), mass_flow_rate(:,1), "Color","b")
plot(time(1:end-1,1), dxval(:,1), "Color","g")
plot(time(1:end-1,1),theoretical_mass_flow_rate(:,1), "Color","r")
legend('Coefficient of Discharge','Mass flow rate','derivative of mass polyfit','Theoretical mass flow rate');
xlabel('Time (s)');
ylabel('Coefficient of Discharge');
title('Coefficient of Discharge over Time');
hold off

% figure(2); hold on
% grid minor
% plot(time(1:end-1,1),-delta_mass(:,1),'r');
% plot(time(1:end-1,1),val(:,1),'b')
% xlabel('Time (s)');
% ylabel('Mass (kg)');
% title('Mass over Time');
% hold off
% 
% figure(3); hold on
% grid minor
% plot(time,pressure);
% xlabel('Time (s)');
% ylabel('Pressure (Pa)');
% title('Pressure over Time');
% hold off




