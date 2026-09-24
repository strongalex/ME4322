%% ME/RBE 4322 Lab 3 - Bathroom Scale Lumped Parameter Model
% Alex Fenwick
% Josh Buttrick
% Torsional model referenced to the dial (pinion) angle theta.

clear; clc; close all;

%% Load and dial calibration
lb2N  = 4.44822;
W_lb  = 130;            % baseline test load
W     = W_lb*lb2N;      % N
m_w   = W/9.81;         % kg
deg_per_lb = 360/300;   % 300 lb per dial revolution

%% Geometry (m)
L1 = 0.2575;    % long lever, A to B (Fig. 26)
L2 = 0.1528;    % short lever, D to C (Fig. 26)
L3 = 0.008;     % bell crank arm to B (Fig. 20)
L4 = 0.018;     % bell crank arm to rack joint S (Fig. 21)
L5 = 0.02;      % A to platform load point (est.)
L7 = 0.128;     % A to C
L8 = L5*L2/L7;  % D to platform load point, level platform
Rp = 0.00175;   % pinion pitch radius, 10 teeth

%% Springs  k = G d^4 / (8 D^3 Na)
G = 79.3e9;     % spring steel, Pa
springk = @(d,OD,Na) G*d^4/(8*(OD-d)^3*Na);

K  = springk(0.9e-3, 8e-3, 8);        % one plate spring (two total)
K3 = springk(2.3e-3, 18.5e-3, 5-2);   % main spring at B, closed ends
K4 = springk(0.469e-3, 6e-3, 40);     % rack return spring

%% Masses (kg) and inertias (est.)
m_L = 0.080;    % one long lever
m_S = 0.050;    % one short lever
m_B = 0.030;    % lever platform B
m_c = 0.005;    % bell crank
m_r = 0.005;    % rack
m_p = 1.000;    % top plate
m_d = 0.020;    % dial disk + pinion
r_d = 0.070;    % dial radius

J_L = m_L*L1^2/3;            % slender rod about A
J_S = m_S*L2^2/3;            % slender rod about D
J_c = m_c*(L3^2+L4^2)/3;     % two short arms about P
J_d = 0.5*m_d*r_d^2;         % disk about its axle

%% Damping, viscous at each joint (est.)
D_A  = 1e-3;    % long lever pivot, N*m*s/rad
D_D  = 1e-3;    % short lever pivot, N*m*s/rad
D_P  = 1e-4;    % bell crank pivot, N*m*s/rad
D_r  = 1e-2;    % rack guide, N*s/m
D_th = 3e-4;    % dial axle, N*m*s/rad

%% Kinematic ratios (d/dtheta of each coordinate)
r_S   = Rp;             % rack translation
r_psi = r_S/L4;         % bell crank rotation
r_B   = L3*r_psi;       % point B translation
r_phA = r_B/L1;         % long lever rotation
r_C   = L7*r_phA;       % point C translation
r_phD = r_C/L2;         % short lever rotation
r_p   = L5*r_phA;       % platform translation

fprintf('Platform consistency L5*phiA / (L8*phiD) = %.3f\n', r_p/(L8*r_phD));

%% Lumped parameters at the dial
Keq = 2*K*r_p^2 + K3*r_B^2 + K4*r_S^2;

Jeq = J_d + J_c*r_psi^2 + m_r*r_S^2 + m_B*r_B^2 ...
    + 2*J_L*r_phA^2 + 2*J_S*r_phD^2 + (m_p + m_w)*r_p^2;

Deq = D_th + D_P*r_psi^2 + D_r*r_S^2 + 2*D_A*r_phA^2 + 2*D_D*r_phD^2;

Teq = W*r_p;

wn   = sqrt(Keq/Jeq);
zeta = Deq/(2*sqrt(Keq*Jeq));
ts   = 4/(zeta*wn);

fprintf('K = %.1f  K3 = %.0f  K4 = %.2f  N/m\n', K, K3, K4);
fprintf('Keq  = %.4e N*m/rad\n', Keq);
fprintf('Jeq  = %.4e kg*m^2  (dial share %.1f%%)\n', Jeq, 100*J_d/Jeq);
fprintf('Deq  = %.4e N*m*s/rad\n', Deq);
fprintf('Teq  = %.4e N*m\n', Teq);
fprintf('wn   = %.2f rad/s   zeta = %.3f   ts = %.2f s\n', wn, zeta, ts);

%% Step response  Jeq*th'' + Deq*th' + Keq*th = Teq
tEnd  = 3;
tspan = linspace(0, tEnd, 3000);
f = @(t,x) [x(2); (Teq - Deq*x(2) - Keq*x(1))/Jeq];
[t, x] = ode45(f, tspan, [0; 0]);
reading = rad2deg(x(:,1))/deg_per_lb;

th_ss = Teq/Keq;
fprintf('Steady state  %.2f deg  ->  %.2f lb (applied %.1f lb)\n', ...
    rad2deg(th_ss), rad2deg(th_ss)/deg_per_lb, W_lb);

figure('Color','w','Position',[100 100 800 450]); hold on
plot(t, reading, 'LineWidth', 1.6, 'DisplayName', 'Model reading');
yline(W_lb, ':', 'LineWidth', 1.4, 'DisplayName', 'Applied load');
xlabel('Time (s)');
ylabel('Dial reading (lbf)');
title(sprintf('Scale dial response to a %g lbf step load', W_lb));
legend('Location', 'southeast');
grid on

%% Static response, 0 to 300 lbf
W_in  = 0:10:300;                        % applied load, lbf
th_in = W_in*lb2N*r_p/Keq;               % steady-state dial angle, rad
W_out = rad2deg(th_in)/deg_per_lb;       % indicated load, lbf

fprintf('Static gain  %.3f lbf indicated per lbf applied\n', W_out(end)/W_in(end));

figure('Color','w','Position',[150 150 600 500]); hold on
plot(W_in, W_out, 'o-', 'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', 'Model');
plot(W_in, W_in, '--', 'LineWidth', 1.2, 'DisplayName', 'Ideal scale');
xlabel('Applied load (lbf)');
ylabel('Indicated load (lbf)');
title('Static response of the scale model');
legend('Location', 'northwest');
xlim([0 300]); ylim([0 300]);
grid on