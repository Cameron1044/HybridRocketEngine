clear; clc; close all

data = importdata("Freshmen_cold flow data.txt");

time_bad = data((9910:10546),1)/1000; % sec
weight_bad = data((9910:10546),2); % lbs
pressure_bad = data((9910:10546),3); % psi
assumed_density = 0.04429207176;

time = time_bad; %(1:305,1);
weight = weight_bad; %(1:305,1);
pressure = pressure_bad; %(1:305,1);

% Injector properties (already given)
num_holes = 5;
diameter_inch = 0.03; % Diameter of one hole in inches
radius_inch = diameter_inch / 2; % Radius of one hole in inches
area_one_hole_in2 = pi * (radius_inch)^2; % Area of one hole in square inches
total_orifice_area_in2 = num_holes * area_one_hole_in2; % Total orifice area in square inches

burn_time = time / 1000; % Convert from ms to s
burn_weight = weight;
burn_pressure = pressure;

mass_flow_rate = diff(burn_weight) ./ diff(burn_time);

fluid_density = assumed_density; % You need to provide this value
volumetric_flow_rate = mass_flow_rate / fluid_density; 

% Calculate theoretical flow rate using Bernoulli's equation
theoretical_flow_rate = total_orifice_area_in2 .* sqrt(2 .* burn_pressure .* fluid_density * 32.2 * 12);

test = mass_flow_rate ./ theoretical_flow_rate(1:end-1);
averageTest = mean(isfinite(test));

calc_flow = averageTest * theoretical_flow_rate;

p = polyfit(theoretical_flow_rate((1:end-1),1), volumetric_flow_rate, 1);
cd_estimated = p(1);

% Display the estimated coefficient of discharge
% disp(['Estimated Coefficient of Discharge (Cd): ', num2str(cd_estimated)]);

figure(1); hold on
plot(burn_time,-burn_weight);
% xlim(220,236);
% ylim(2412,30);
xlabel("time");
ylabel("weight");
grid on

figure(2); hold on
plot(burn_time,burn_pressure);
xlabel("time");
ylabel("pressure"); 

figure(3); hold on
plot(burn_time, theoretical_flow_rate)
plot(burn_time, calc_flow);
plot(burn_time(1:end-1), mass_flow_rate)


