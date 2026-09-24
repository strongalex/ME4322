clc
clear
close all
tstart = tic;
A = [0,0,0];
B0 = [-.002,3.81,0];
C0 = [-14.758,0,0];
L_AB = norm(B0-A);
L_BC = norm(C0-B0);

% Run conditions
rpm    = 13;   % crank speed, tachometer / Hall timing
e      = 0;      % slider offset from crank centerline
branch = -1;     % -1 puts the slider on the -x side

omega2 = rpm*2*pi/60;
theta  = (0:359)';
th     = deg2rad(theta);
t      = theta/(360*rpm/60);   % uniform timestep over one revolution

% Crank pin
Bx = L_AB*cos(th);
By = L_AB*sin(th);

% Circle-line intersection: circle of radius L_BC about B, line y = e
disc = L_BC^2 - (e - By).^2;
if any(disc < 0)
    error('No intersection. Check L_AB, L_BC, or e.')
end
Cx = Bx + branch*sqrt(disc);
Cy = e*ones(size(Cx));

% Connecting rod angle, measured B to C
th3 = atan2(Cy - By, Cx - Bx);

% Velocity
omega3 = -L_AB*cos(th)*omega2 ./ (L_BC*cos(th3));
vC     = -L_AB*sin(th)*omega2 - L_BC*sin(th3).*omega3;

% Acceleration, constant crank speed
alpha3 = (L_AB*sin(th)*omega2^2 + L_BC*sin(th3).*omega3.^2) ./ (L_BC*cos(th3));
aC     = -L_AB*cos(th)*omega2^2 - L_BC*(cos(th3).*omega3.^2 + sin(th3).*alpha3);

% Checks
fprintf('L_AB = %.4f, L_BC = %.4f\n', L_AB, L_BC);
fprintf('Stroke = %.4f (2*L_AB = %.4f)\n', max(Cx)-min(Cx), 2*L_AB);
fprintf('Cycle time = %.4f s, peak |omega3| = %.4f rad/s\n', t(end)+t(2), max(abs(omega3)));

% Export for the RMSE comparison
results = table(t, theta, Cx*10, rad2deg(th3), omega3, vC*10, alpha3, aC*10, ...
    'VariableNames', {'Time_s','Theta2_deg','Cx_mm','Theta3_deg', ...
                      'Omega3_rads','vC_mmps','Alpha3_radss','aC_mmpss'});
writetable(results,'CrankSliderTheoretical.csv');

figure
subplot(2,1,1)
plot(t, omega3,'LineWidth',1.2)
grid on
xlabel('Time (s)'), ylabel('\omega_3 (rad/s)')
title('Connecting rod angular velocity')

subplot(2,1,2)
plot(t, Cx*10,'LineWidth',1.2)
grid on
xlabel('Time (s)'), ylabel('C_x (mm)')
title('Slider position')

toc(tstart)