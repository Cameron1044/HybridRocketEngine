clear; clc; close all;

% Load data
data = importdata("Characteristic Fire 1_30_6_7/Freshmen_cold flow data.txt");

time = data(:,1) / 1000; % Convert ms to seconds
weight = data(:,2); % Assuming weight is in pounds (lb)
pressure = data(:,3); % Assuming pressure is in psi

% Constants for unit conversion
g = 32.2; % Gravity in ft/s^2
lb_to_kg = 0.45359237; % Pounds to kilograms conversion
psi_to_pa = 6894.757; % PSI to Pascal conversion
in2_to_ft2 = 1 / 144; % Square inches to square feet conversion
in2_to_m2 = 0.00064516; % Square inches to square meters conversion

% Injector properties (given)
num_holes = 5;
diameter_inch = 0.03; % Diameter of one hole in inches
radius_inch = diameter_inch / 2; % Radius of one hole in inches
area_one_hole_in2 = pi * (radius_inch)^2; % Area of one hole in square inches
total_orifice_area_in2 = num_holes * area_one_hole_in2; % Total orifice area in square inches
total_orifice_area_m2 = total_orifice_area_in2 * in2_to_m2; % Convert to square meters

% Filter data for burn time
burn_time = time((9835:10546)); % Already in seconds
burn_weight = weight((9835:10546));
burn_pressure = pressure((9835:10546));

% Remove any negative pressures
burn_pressure(burn_pressure < 0) = 0;

% Convert weight loss to mass flow rate (kg/s)
mass_flow_rate = diff(burn_weight) * lb_to_kg ./ diff(burn_time);

% Convert pressure to Pascals
burn_pressure_pa = burn_pressure * psi_to_pa;

% Calculate volumetric flow rate (m^3/s)
% Assuming fluid density is in slugs/ft^3 and converting to kg/m^3
assumed_density_lb_s3 = 0.04429207176; % Density in slugs per cubic foot
assumed_density_kg_m3 = assumed_density_lb_s3 * 27680; % Convert slugs/ft^3 to kg/m^3
volumetric_flow_rate = mass_flow_rate / assumed_density_kg_m3;

fluid_density = assumed_density_kg_m3; % You need to provide this value
theoretical_flow_rate = total_orifice_area_in2 .* sqrt(2 .* burn_pressure .* fluid_density);

test = mass_flow_rate ./ theoretical_flow_rate(1:711);
averageTest = mean(isfinite(test));
disp('Average Test for Cd');
disp(averageTest);

calc_flow = averageTest * theoretical_flow_rate;

p = polyfit(theoretical_flow_rate((1:711),1), volumetric_flow_rate, 1);
cd_estimated = p(1);
disp('Average Test for Cd');
disp(cd_estimated);

figure(1); hold on
plot(burn_time,-burn_weight);
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
plot(burn_time(1:711), mass_flow_rate)

