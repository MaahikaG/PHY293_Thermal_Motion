%Part 1
%Constant Values
conversion_factor = 0.1155;
T = 296.5;
T_unc = 0.5;
viscosity = 1.00*10^-3;
viscosity_unc = 0.05*10^-3;
radius = (1.9*10^-6) / 2; %in meters
radius_unc = 0.1*10^-16;
stokes_drag = 6 * pi * viscosity * radius;
k_actual = 1.38*10^-23;


%loading file and accessing x_vals and y_vals
filenames = {'p.1.3.txt', 'p1.2.txt', 'p2.1.txt', 'p2.2.txt', 'p2.3.txt', ...
             'p3.1.txt', 'p3.3.txt', 'p4.1.txt', 'p4.2.txt', ...
             'p4.3.txt', 'p4.4.txt'};
k_calc_vals = zeros(length(filenames), 1);
percent_error_vals = zeros(length(filenames), 1);
k_calc_uncertainties = zeros(length(filenames), 1);
r_vals = [];

for j = 1:length(filenames)
    filename = filenames{j};
    data = readmatrix(filename);
    x_vals = data(:,1);
    y_vals = data(:,2);
    
    %calculating squared displacement
    squared_displacement = zeros(length(x_vals), 1);
    squared_displacement(1) = 0;
    for i = 1:(length(x_vals)-1)
        x_displacement = (x_vals(i+1) - x_vals(i)) * conversion_factor;
        y_displacement = (y_vals(i+1) - y_vals(i)) * conversion_factor;
        step_size = sqrt (x_displacement^2 + y_displacement^2);
        squared_displacement(i+1) = squared_displacement(i) + x_displacement^2 + y_displacement^2;
        r_vals = [r_vals; step_size];
    end
    
    time_intervals = (1:length(squared_displacement)) / 2; %2 images per second
    
    %plotting squared displacement vs time
    % figure;
    % plot (time_intervals, squared_displacement);
    
    %Using a linear fit model
    mdl = fitlm(time_intervals', squared_displacement);
    
    % Access coefficients and R-squared value
    coefficients = mdl.Coefficients;
    r2 = mdl.Rsquared.Adjusted;
    D = coefficients.Estimate(2);
    D_SI_units = D * 10^-12;
    k_calc_vals(j) = (D_SI_units*stokes_drag)/T;
    percent_error_vals(j) = (abs(k_calc_vals(j)-k_actual)/k_actual) * 100;
    r_squared = mdl.Rsquared.Ordinary;  % Regular R-squared
    r_squared_adjusted = mdl.Rsquared.Adjusted;  % Adjusted R-squared
    D_standard_error = mdl.Coefficients.SE(2) * 10^-12; 
    k_calc_uncertainties(j) = sqrt ((((6*pi*radius*viscosity)/T)*D_standard_error)^2 + ...
                                (((D_SI_units*6*pi*viscosity)/T)*radius_unc)^2 + ...
                                (((D_SI_units*6*pi*radius)/T)*viscosity_unc)^2 + ...
                                (((-D_SI_units*6*pi*radius*viscosity)/T^2)*T_unc)^2);
    
    figure;
    plot(mdl);
    title(['Filename: ', filenames{j}, ' Percent Error: ', num2str(percent_error_vals(j)), ' K = ', num2str(k_calc_vals(j)), ' +- ', num2str(k_calc_uncertainties(j))]);
    xlabel ('Time');
    ylabel ('Distance Travelled');
    % disp(mdl);
    disp (['The k value from ', filenames{j}, ' is ', num2str(k_calc_vals(j))]);
    disp (['The percent uncertainty from ', filenames{j}, ' is ', num2str(percent_error_vals(j))]);

end
% disp (k_calc_vals);
% disp (percent_error_vals);
disp (['The average percent error is for the linear fit model prediction is ', num2str(mean(percent_error_vals))]);


%Part 2
%Distance travelled
N = length(r_vals);
bins = (round (sqrt (N)));

figure;
histfit(r_vals, bins, 'rayleigh'); 
title('Histogram with Rayleigh Fit');
xlabel('Step Size');
ylabel('Probability Density');

pd = fitdist(r_vals, 'Rayleigh');  % 'pd' is a probability distribution object
sigma_estimated = pd.B;  % 'B' is the sigma (scale) parameter of the Rayleigh distribution

% Calculate the diffusion coefficient D
D_estimated = sigma_estimated^2 / (2 * 0.5);
D_estimated_SI_units = D_estimated * 10^-12;

% Calculate the Boltzmann constant k using the Stokes-Einstein equation
k_B_estimated = (D_estimated_SI_units * stokes_drag) / T;
percent_error_part2 = (abs(k_B_estimated-k_actual)/k_actual) * 100;
disp (['The estimated k value from the Rayleigh distribution is ', num2str(k_B_estimated)]);
disp (['The percent error is ', num2str(percent_error_part2)]);















