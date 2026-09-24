%% ME/RBE 4322 Lab 3 - Bathroom Scale Lumped Parameter Model
% Alex Fenwick
% Torsional model referenced to the dial (pinion) angle theta.
% Values marked TODO are estimates. Replace with measurements.

clear; clc; close all;

%% Load
lb2N  = 4.44822;
W_lb  = 5;              % laptop test load
W     = W_lb*lb2N;      % N
m_w   = W/9.81;         % kg

%% Geometry (m)
L1 = 0.2575;    % long lever length, A to B (Fig. 26)
L2 = 0.1528;    % short lever length, D to C (Fig. 26)
L3 = 0.008;     % pin to pivot, bell crank arm to B (Fig. 20)
L4 = 0.018;     % pivot to rack joint S (Fig. 21)
L5 = 0.02;     % A to platform load point, long lever   
L7 = 0.128;     % A to C, short lever tip on long lever   
L8 = 0.02;     % D to platform load point, short lever  
Rp = 0.00175;   % pinion pitch radius, 10 teeth  

%% Springs  k = G d^4 / (8 D^3 N)
G = 79.3e9;     % spring steel, Pa
springk = @(d,OD,N) G*d^4/(8*(OD-d)^3*N);

K  = springk(0.9e-3, 8e-3, 8);    % one plate spring (two total)  
K3 = springk(2.0e-3, 18.5e-3, 5);   % main spring at B             
K4 = springk(0.469e-3, 6e-3, 40);   % rack return spring     

%% Masses (kg) and inertias
m_L = 0.080;    % one long lever        TODO
m_S = 0.050;    % one short lever       TODO
m_B = 0.030;    % lever platform B      TODO
m_c = 0.005;    % bell crank            TODO
m_r = 0.005;    % rack                  TODO
m_p = 1.000;    % top plate             TODO
m_d = 0.020;    % dial disk + pinion    TODO
r_d = 0.070;    % dial radius           TODO

J_L = m_L*L1^2/3;            % slender rod about A
J_S = m_S*L2^2/3;            % slender rod about D
J_c = m_c*(L3^2+L4^2)/3;     % two short arms about P
J_d = 0.5*m_d*r_d^2;         % disk about its axle

%% Damping (assumed viscous at each joint)
D_A  = 1e-3;    % long lever pivot, N*m*s/rad     TODO
D_D  = 1e-3;    % short lever pivot, N*m*s/rad    TODO
D_P  = 1e-4;    % bell crank pivot, N*m*s/rad     TODO
D_r  = 1e-2;    % rack guide, N*s/m               TODO
D_th = 3e-4;    % dial axle, N*m*s/rad            TODO

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
Keq = 2*K*r_p^2 + K3*r_B^2 + K4*r_S^2;   % plate springs act on top plate

Jeq = J_d + J_c*r_psi^2 + m_r*r_S^2 + m_B*r_B^2 ...
    + 2*J_L*r_phA^2 + 2*J_S*r_phD^2 + (m_p + m_w)*r_p^2;

Deq = D_th + D_P*r_psi^2 + D_r*r_S^2 + 2*D_A*r_phA^2 + 2*D_D*r_phD^2;

Teq = W*r_p;

wn   = sqrt(Keq/Jeq);
zeta = Deq/(2*sqrt(Keq*Jeq));

fprintf('K = %.1f  K3 = %.1f  K4 = %.2f  N/m\n', K, K3, K4);
fprintf('Keq  = %.4e N*m/rad\n', Keq);
fprintf('Jeq  = %.4e kg*m^2\n', Jeq);
fprintf('Deq  = %.4e N*m*s/rad\n', Deq);
fprintf('Teq  = %.4e N*m\n', Teq);
fprintf('wn   = %.2f rad/s   zeta = %.3f\n', wn, zeta);

%% Simulation  Jeq*th'' + Deq*th' + Keq*th = Teq
th0  = [0; 0];
tEnd = 6;
f = @(t,x) [x(2); (Teq - Deq*x(2) - Keq*x(1))/Jeq];
[t, x] = ode45(f, [0 tEnd], th0);

deg_per_lb = 360/300;          % 300 lb per dial revolution
th_deg  = rad2deg(x(:,1));
reading = th_deg/deg_per_lb;

th_ss = Teq/Keq;
fprintf('Steady state  %.2f deg  ->  %.2f lb (expected %.1f lb)\n', ...
    rad2deg(th_ss), rad2deg(th_ss)/deg_per_lb, W_lb);

%% Plot
figure('Color','w','Position',[100 100 800 450]);
yyaxis left
plot(t, th_deg, 'LineWidth', 1.8); hold on
ylabel('Dial angle \theta (deg)');
yyaxis right
plot(t, reading, '--', 'LineWidth', 1.2);
yline(W_lb, ':', 'Actual load', 'LineWidth', 1.2);
ylabel('Dial reading (lb)');
xlabel('Time (s)');
title(sprintf('Scale dial response to a %g lb step load', W_lb));
legend('Model \theta(t)', 'Model reading', 'Actual load', 'Location', 'southeast');
grid on
