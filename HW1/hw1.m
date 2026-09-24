%% ME/RBE 4322 Homework 1 - Six-bar pick-and-place linkage
% Position, velocity, acceleration, static equilibrium and Newton's second
% law over a full 360 deg rotation of the input crank AB.
%
% Topology
%   Link 1  AB      crank, grounded at A, driven by the motor
%   Link 2  BC      coupler
%   Link 3  D-C-E   ternary link grounded at D (C and E are the other joints)
%   Link 4  EF      coupler
%   Link 5  GF      output rocker, grounded at G, carries the gripper
%
%   Loop 1  A -> B -> C -> D -> A
%   Loop 2  D -> E -> F -> G -> D

clc
clear
close all

tstart = tic;

%% ---------------- Assumptions ----------------
% 1. Planar mechanism, rigid links, frictionless revolute pins.
% 2. Links are 4340 steel bars of the given cross section (width 10 cm,
%    thickness 5 cm, bore 6 cm) with semicircular ends of radius w/2.
% 3. Link 3 is one straight bar from D to E with an extra bore at C
%    (the three given coordinates are collinear to within 5.5 mm).
% 4. COM of every link is at its geometric centre. Link 5 is the full
%    bar G-F-H, so its COM is the midpoint of G and H.
% 5. The artifact is a 200 N point mass rigidly gripped by link 5.
% 6. One pick-and-place per revolution of the input crank.
% 7. The artifact is carried for the whole revolution (worst case).
% 8. Input crank turns at constant speed, so alpha_AB = 0.
% 9. Gripper mass is ignored.

g     = 9.81;
rho   = 7850;     % steel density
wSec  = 0.10;     % link width
tSec  = 0.05;     % link thickness
dBore = 0.06;     % joint bore diameter

%% ---------------- First position geometry ----------------
A0 = [ 1.400 0.485 0];
B0 = [ 1.670 0.990 0];
C0 = [ 0.255 1.035 0];
D0 = [ 0.285 0.055 0];
E0 = [ 0.195 2.540 0];
F0 = [-0.980 2.570 0];
G0 = [ 0.050 0.200 0];

A = A0;  D = D0;  G = G0;          % grounded joints never move

L_AB = norm(B0-A0);
L_BC = norm(C0-B0);
L_DC = norm(C0-D0);
L_DE = norm(E0-D0);
L_EF = norm(F0-E0);
L_GF = norm(F0-G0);

% Load point H on link 5, 1.843 m from F, taken along GF extended
rGrip = 1.843;
gripOffsetDeg = 0;
H0 = F0 + rGrip*rotVec((F0-G0)/L_GF, deg2rad(gripOffsetDeg));
L_GH = norm(H0-G0);            % link 5 is the straight bar G-F-H

%% ---------------- Throughput requirement ----------------
partsRequired = 12500;
shiftSeconds  = 9*3600;
cycleTime = shiftSeconds/partsRequired;
omegaIn   = 2*pi/cycleTime;                    % rad/s, CCW positive
alphaIn   = 0;

fprintf('Cycle time            : %.4f s per part\n', cycleTime);
fprintf('Input angular velocity: %.4f rad/s (%.2f rpm)\n\n', omegaIn, omegaIn*60/(2*pi));


%% ---------------- Mass and mass moment of inertia ----------------
linkIner = @(m,L) m*((L+wSec)^2 + wSec^2)/12;


m1 = 23.3396; 
m2 = 56.5595; 
m3 = 97.3536; 
m4 = 46.9969; 
m5 = 173.5188;


J1 = linkIner(m1,L_AB);
J2 = linkIner(m2,L_BC);
J3 = linkIner(m3,L_DE);
J4 = linkIner(m4,L_EF);
J5 = linkIner(m5,L_GH);

Wart = 200;
mArt = Wart/g;

fprintf('Link masses  [kg]     : %.3f %.3f %.3f %.3f %.3f\n', m1,m2,m3,m4,m5);
fprintf('Link inertia [kg m^2] : %.3f %.3f %.3f %.3f %.3f\n', J1,J2,J3,J4,J5);
fprintf('Artifact mass [kg]    : %.3f\n\n', mArt);

mass = [m1 m2 m3 m4 m5];
Jcm  = [J1 J2 J3 J4 J5];

%% ---------------- Sweep the input crank ----------------
nStep  = 361;
sweep  = linspace(0,2*pi,nStep);
th1_0  = atan2(B0(2)-A0(2), B0(1)-A0(1));
phi3_0 = atan2(C0(2)-D0(2), C0(1)-D0(1));
phi5_0 = atan2(F0(2)-G0(2), F0(1)-G0(1));

Bp = zeros(nStep,3); Cp = Bp; Ep = Bp; Fp = Bp; Hp = Bp;
omegaLink = zeros(nStep,5);      % [AB BC DCE EF GF]
alphaLink = zeros(nStep,5);
vJoint = zeros(nStep,3,5);    % B C E F H
aJoint = zeros(nStep,3,5);
vCom   = zeros(nStep,3,5);
aCom   = zeros(nStep,3,5);
Fstat  = zeros(nStep,15);
Fdyn   = zeros(nStep,15);

Cprev = C0;  Fprev = F0;

for k = 1:nStep
    th1 = th1_0 + sweep(k);

    % --- position ---
    B = A + L_AB*[cos(th1) sin(th1) 0];
    C = circleIntersect(B, L_BC, D, L_DC, Cprev);
    dphi3 = wrapPi( atan2(C(2)-D(2), C(1)-D(1)) - phi3_0 );
    E = rotAbout(E0, D, dphi3);
    F = circleIntersect(E, L_EF, G, L_GF, Fprev);
    dphi5 = wrapPi( atan2(F(2)-G(2), F(1)-G(1)) - phi5_0 );
    H = rotAbout(H0, G, dphi5);
    Cprev = C;  Fprev = F;

    % --- centres of mass ---
    S1 = (A+B)/2;  S2 = (B+C)/2;  S3 = (D+E)/2;  S4 = (E+F)/2;  S5 = (G+H)/2;

    % --- angular velocity, vector loop closure ---
    % Loop 1: w1 x (B-A) + w2 x (C-B) + w3 x (D-C) = 0
    M1 = [ perp(C-B).' perp(D-C).' ];
    r1 = -cross([0 0 omegaIn], B-A).';
    s1 = M1(1:2,:)\r1(1:2);
    w2 = s1(1);  w3 = s1(2);

    % Loop 2: w3 x (E-D) + w4 x (F-E) + w5 x (G-F) = 0
    M2 = [ perp(F-E).' perp(G-F).' ];
    r2 = -cross([0 0 w3], E-D).';
    s2 = M2(1:2,:)\r2(1:2);
    w4 = s2(1);  w5 = s2(2);

    W1v=[0 0 omegaIn]; W2v=[0 0 w2]; W3v=[0 0 w3]; W4v=[0 0 w4]; W5v=[0 0 w5];

    % --- angular acceleration, same loops differentiated ---
    b1 = -( cross([0 0 alphaIn],B-A) - omegaIn^2*(B-A) - w2^2*(C-B) - w3^2*(D-C) );
    s3 = M1(1:2,:)\b1(1:2).';
    a2 = s3(1);  a3 = s3(2);

    b2 = -( cross([0 0 a3],E-D) - w3^2*(E-D) - w4^2*(F-E) - w5^2*(G-F) );
    s4 = M2(1:2,:)\b2(1:2).';
    a4 = s4(1);  a5 = s4(2);

    A1v=[0 0 alphaIn]; A2v=[0 0 a2]; A3v=[0 0 a3]; A4v=[0 0 a4]; A5v=[0 0 a5];

    % --- joint velocity and acceleration ---
    vB = cross(W1v,B-A);
    vC = cross(W3v,C-D);
    vE = cross(W3v,E-D);
    vF = cross(W5v,F-G);
    vH = cross(W5v,H-G);

    aB = cross(A1v,B-A) - omegaIn^2*(B-A);
    aC = cross(A3v,C-D) - w3^2*(C-D);
    aE = cross(A3v,E-D) - w3^2*(E-D);
    aF = cross(A5v,F-G) - w5^2*(F-G);
    aH = cross(A5v,H-G) - w5^2*(H-G);

    % --- centre of mass velocity and acceleration ---
    vS1 = cross(W1v,S1-A);
    vS2 = vB + cross(W2v,S2-B);
    vS3 = cross(W3v,S3-D);
    vS4 = vE + cross(W4v,S4-E);
    vS5 = cross(W5v,S5-G);

    aS1 = cross(A1v,S1-A) - omegaIn^2*(S1-A);
    aS2 = aB + cross(A2v,S2-B) - w2^2*(S2-B);
    aS3 = cross(A3v,S3-D) - w3^2*(S3-D);
    aS4 = aE + cross(A4v,S4-E) - w4^2*(S4-E);
    aS5 = cross(A5v,S5-G) - w5^2*(S5-G);

    pts  = [A;B;C;D;E;F;G;H];
    Scm  = [S1;S2;S3;S4;S5];
    aCm  = [aS1;aS2;aS3;aS4;aS5];
    alp  = [alphaIn a2 a3 a4 a5];

    % --- static equilibrium (all accelerations set to zero) ---
    Fart_s = [0 -Wart 0];
    Fstat(k,:) = solveLinkForces(pts,Scm,zeros(5,3),zeros(1,5),mass,Jcm,g,Fart_s).';

    % --- Newton's second law ---
    % Reaction of the gripped artifact on link 5, weight less its inertia
    Fart_d = [0 -Wart 0] - mArt*aH;
    Fdyn(k,:) = solveLinkForces(pts,Scm,aCm,alp,mass,Jcm,g,Fart_d).';

    % --- store ---
    Bp(k,:)=B; Cp(k,:)=C; Ep(k,:)=E; Fp(k,:)=F; Hp(k,:)=H;
    omegaLink(k,:) = [omegaIn w2 w3 w4 w5];
    alphaLink(k,:) = alp;
    vJoint(k,:,:) = [vB.' vC.' vE.' vF.' vH.'];
    aJoint(k,:,:) = [aB.' aC.' aE.' aF.' aH.'];
    vCom(k,:,:)   = [vS1.' vS2.' vS3.' vS4.' vS5.'];
    aCom(k,:,:)   = [aS1.' aS2.' aS3.' aS4.' aS5.'];
end


%% ---------------- First position results ----------------
lbl = {'AB','BC','DCE','EF','GF'};
fprintf('--- FIRST POSITION KINEMATICS ---\n');
fprintf('Angular velocity [rad/s]\n');
for i=1:5, fprintf('  omega %-4s %10.5f\n', lbl{i}, omegaLink(1,i)); end
fprintf('Angular acceleration [rad/s^2]\n');
for i=1:5, fprintf('  alpha %-4s %10.5f\n', lbl{i}, alphaLink(1,i)); end

jl = {'B','C','E','F','H'};
fprintf('Joint velocity [m/s]  and acceleration [m/s^2]\n');
for i=1:5
    fprintf('  %-2s v = [%9.5f %9.5f]   a = [%9.5f %9.5f]\n', jl{i}, ...
        vJoint(1,1,i), vJoint(1,2,i), aJoint(1,1,i), aJoint(1,2,i));
end
fprintf('Centre of mass velocity and acceleration\n');
for i=1:5
    fprintf('  S%d %-4s v = [%9.5f %9.5f]   a = [%9.5f %9.5f]\n', i, lbl{i}, ...
        vCom(1,1,i), vCom(1,2,i), aCom(1,1,i), aCom(1,2,i));
end

fn = {'FAx','FAy','FBx','FBy','FCx','FCy','FDx','FDy', ...
      'FEx','FEy','FFx','FFy','FGx','FGy','Tin'};
fprintf('\n--- FIRST POSITION FORCES [N] AND TORQUE [N m] ---\n');
fprintf('%-6s %14s %14s\n','','Static','Newton');
for i=1:15
    fprintf('%-6s %14.4f %14.4f\n', fn{i}, Fstat(1,i), Fdyn(1,i));
end

%% ---------------- Symbolic check at the first position ----------------
% Same five free bodies assembled with the symbolic toolbox
runSymbolic = true;
if runSymbolic
    B=Bp(1,:); C=Cp(1,:); E=Ep(1,:); F=Fp(1,:); H=Hp(1,:);
    S1=(A+B)/2; S2=(B+C)/2; S3=(D+E)/2; S4=(E+F)/2; S5=(G+H)/2;

    WAB=[0 -m1*g 0]; WBC=[0 -m2*g 0]; WDCE=[0 -m3*g 0];
    WEF=[0 -m4*g 0]; WGF=[0 -m5*g 0];

    syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin real
    ForceA=[FAx FAy 0]; ForceB=[FBx FBy 0]; ForceC=[FCx FCy 0];
    ForceD=[FDx FDy 0]; ForceE=[FEx FEy 0]; ForceF=[FFx FFy 0];
    ForceG=[FGx FGy 0]; InputTorque=[0 0 Tin];

    % ForceX is what the downstream link applies to the upstream link
    AppliedForce = [0 -Wart 0];

    eqn1 = ForceA + ForceB + WAB == 0;
    eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + InputTorque == 0;
    eqn3 = -ForceB + ForceC + WBC == 0;
    eqn4 = cross(B-S2,-ForceB) + cross(C-S2,ForceC) == 0;
    eqn5 = -ForceC + ForceD + ForceE + WDCE == 0;
    eqn6 = cross(C-S3,-ForceC) + cross(D-S3,ForceD) + cross(E-S3,ForceE) == 0;
    eqn7 = -ForceE + ForceF + WEF == 0;
    eqn8 = cross(E-S4,-ForceE) + cross(F-S4,ForceF) == 0;
    eqn9 = -ForceF + ForceG + WGF + AppliedForce == 0;
    eqn10= cross(F-S5,-ForceF) + cross(G-S5,ForceG) + cross(H-S5,AppliedForce) == 0;

    unk = [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin];
    symStatic = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10], unk);
    symVec = double([symStatic.FAx symStatic.FAy symStatic.FBx symStatic.FBy ...
        symStatic.FCx symStatic.FCy symStatic.FDx symStatic.FDy ...
        symStatic.FEx symStatic.FEy symStatic.FFx symStatic.FFy ...
        symStatic.FGx symStatic.FGy symStatic.Tin]);
    fprintf('\nSymbolic vs matrix statics, max difference: %.3e\n', ...
        max(abs(symVec - Fstat(1,:))));
end

%% ---------------- Motion range and load path ----------------
sw3 = rad2deg(max(unwrap(atan2(Cp(:,2)-D(2),Cp(:,1)-D(1)))) - ...
              min(unwrap(atan2(Cp(:,2)-D(2),Cp(:,1)-D(1)))));
sw5 = rad2deg(max(unwrap(atan2(Fp(:,2)-G(2),Fp(:,1)-G(1)))) - ...
              min(unwrap(atan2(Fp(:,2)-G(2),Fp(:,1)-G(1)))));
fprintf('\nAngular sweep of link 3 (DCE): %.2f deg\n', sw3);
fprintf('Angular sweep of link 5 (GF) : %.2f deg\n', sw5);

%% ---------------- Extremes over the revolution ----------------
jn = {'A','B','C','D','E','F','G'};
fprintf('\n--- JOINT FORCE MAGNITUDES OVER THE FULL REVOLUTION [N] ---\n');
fprintf('%-4s %12s %12s %12s %12s\n','','stat max','stat min','dyn max','dyn min');
for i=1:7
    ms = hypot(Fstat(:,2*i-1),Fstat(:,2*i));
    md = hypot(Fdyn(:,2*i-1), Fdyn(:,2*i));
    fprintf('%-4s %12.2f %12.2f %12.2f %12.2f\n', jn{i}, max(ms),min(ms),max(md),min(md));
end
[tsMax,iS] = max(abs(Fstat(:,15)));
[tdMax,iD] = max(abs(Fdyn(:,15)));
fprintf('\nPeak static torque : %.2f N m at input angle %.0f deg\n', tsMax, rad2deg(sweep(iS)));
fprintf('Peak dynamic torque: %.2f N m at input angle %.0f deg\n', tdMax, rad2deg(sweep(iD)));
fprintf('RMS dynamic torque : %.2f N m\n', sqrt(mean(Fdyn(:,15).^2)));

%% ---------------- Export kinematics ----------------
degIn = rad2deg(sweep).';
Tk = table(degIn, omegaLink(:,2),omegaLink(:,3),omegaLink(:,4),omegaLink(:,5), ...
                  alphaLink(:,2),alphaLink(:,3),alphaLink(:,4),alphaLink(:,5), ...
                  Bp(:,1),Bp(:,2),Cp(:,1),Cp(:,2),Ep(:,1),Ep(:,2),Fp(:,1),Fp(:,2), ...
    'VariableNames',{'inputDeg','w_BC','w_DCE','w_EF','w_GF', ...
                     'a_BC','a_DCE','a_EF','a_GF', ...
                     'Bx','By','Cx','Cy','Ex','Ey','Fx','Fy'});
writetable(Tk,'kinematics_matlab.csv');
fprintf('\nWrote kinematics_matlab.csv\n');
outFile = fullfile(pwd,'kinematics_matlab.csv');
writetable(Tk,outFile);
fprintf('Wrote %s\n', outFile);

%% ---------------- Plots ----------------
deg = rad2deg(sweep);

figure('Name','Kinematic outline, first position');
lk = [A0;B0]; plot(lk(:,1),lk(:,2),'-o','LineWidth',2); hold on
plot([B0(1) C0(1)],[B0(2) C0(2)],'-o','LineWidth',2);
plot([D0(1) C0(1) E0(1)],[D0(2) C0(2) E0(2)],'-o','LineWidth',2);
plot([E0(1) F0(1)],[E0(2) F0(2)],'-o','LineWidth',2);
plot([G0(1) F0(1) H0(1)],[G0(2) F0(2) H0(2)],'-o','LineWidth',2);
nm = {'A','B','C','D','E','F','G','H'};
pt = [A0;B0;C0;D0;E0;F0;G0;H0];
for i=1:8
    text(pt(i,1)+0.05,pt(i,2)+0.05,sprintf('%s (%.3f, %.3f)',nm{i},pt(i,1),pt(i,2)));
end
plot(pt([1 4 7],1),pt([1 4 7],2),'k^','MarkerSize',10,'MarkerFaceColor','k');
axis equal; grid on; xlabel('x (m)'); ylabel('y (m)');
title('Kinematic outline at the first position');
legend({'Link 1 AB','Link 2 BC','Link 3 DCE','Link 4 EF','Link 5 GFH'}, ...
    'Location','best');

figure('Name','Joint paths');
plot(Bp(:,1),Bp(:,2),Cp(:,1),Cp(:,2),Ep(:,1),Ep(:,2),Fp(:,1),Fp(:,2),Hp(:,1),Hp(:,2),'LineWidth',1.5);
axis equal; grid on; xlabel('x (m)'); ylabel('y (m)');
legend('B','C','E','F','H'); title('Joint paths over one revolution');

figure('Name','Joint positions');
subplot(2,1,1)
plot(deg,[Bp(:,1) Cp(:,1) Ep(:,1) Fp(:,1) Hp(:,1)],'LineWidth',1.2); grid on
ylabel('x (m)'); legend('B','C','E','F','H'); title('Joint positions');
subplot(2,1,2)
plot(deg,[Bp(:,2) Cp(:,2) Ep(:,2) Fp(:,2) Hp(:,2)],'LineWidth',1.2); grid on
xlabel('input crank angle (deg)'); ylabel('y (m)');

figure('Name','Angular velocity');
plot(deg,omegaLink,'LineWidth',1.4); grid on
xlabel('input crank angle (deg)'); ylabel('\omega (rad/s)');
legend(lbl); title('Link angular velocities');

figure('Name','Angular acceleration');
plot(deg,alphaLink,'LineWidth',1.4); grid on
xlabel('input crank angle (deg)'); ylabel('\alpha (rad/s^2)');
legend(lbl); title('Link angular accelerations');

figure('Name','Joint speed');
plot(deg,squeeze(hypot(vJoint(:,1,:),vJoint(:,2,:))),'LineWidth',1.4); grid on
xlabel('input crank angle (deg)'); ylabel('|v| (m/s)');
legend(jl); title('Joint speed magnitudes');

figure('Name','Joint acceleration');
plot(deg,squeeze(hypot(aJoint(:,1,:),aJoint(:,2,:))),'LineWidth',1.4); grid on
xlabel('input crank angle (deg)'); ylabel('|a| (m/s^2)');
legend(jl); title('Joint acceleration magnitudes');

figure('Name','Static joint forces');
plot(deg,hypot(Fstat(:,1:2:13),Fstat(:,2:2:14)),'LineWidth',1.4); grid on
xlabel('input crank angle (deg)'); ylabel('|F| (N)');
legend(jn); title('Static joint force magnitudes');

figure('Name','Dynamic joint forces');
plot(deg,hypot(Fdyn(:,1:2:13),Fdyn(:,2:2:14)),'LineWidth',1.4); grid on
xlabel('input crank angle (deg)'); ylabel('|F| (N)');
legend(jn); title('Dynamic joint force magnitudes');

figure('Name','Input torque');
plot(deg,Fstat(:,15),'LineWidth',1.5); hold on
plot(deg,Fdyn(:,15),'LineWidth',1.5); grid on
xlabel('input crank angle (deg)'); ylabel('T_{in} (N m)');
legend('static equilibrium','Newton''s second law');
title('Motor torque required');

%% ---------------- Local functions ----------------
function x = solveLinkForces(pts,Scm,aCm,alp,mass,Jcm,g,Fart)
% 15 unknowns [FAx FAy ... FGx FGy Tin], 15 equations from 5 free bodies.
% Zero aCm and alp gives static equilibrium.
A=pts(1,:); B=pts(2,:); C=pts(3,:); D=pts(4,:); E=pts(5,:);
F=pts(6,:); G=pts(7,:); H=pts(8,:);
S1=Scm(1,:); S2=Scm(2,:); S3=Scm(3,:); S4=Scm(4,:); S5=Scm(5,:);
m1=mass(1); m2=mass(2); m3=mass(3); m4=mass(4); m5=mass(5);

mc = @(Q,S) [-(Q(2)-S(2)), (Q(1)-S(1))];   % moment coefficients on [Fx Fy]
mz = @(r,Fv) r(1)*Fv(2) - r(2)*Fv(1);

iA=1:2; iB=3:4; iC=5:6; iD=7:8; iE=9:10; iF=11:12; iG=13:14; iT=15;
M = zeros(15); b = zeros(15,1);

% Link 1 AB
M(1,iA(1))=1; M(1,iB(1))=1;                  b(1)=m1*aCm(1,1);
M(2,iA(2))=1; M(2,iB(2))=1;                  b(2)=m1*aCm(1,2)+m1*g;
M(3,iA)=mc(A,S1); M(3,iB)=mc(B,S1); M(3,iT)=1;
                                             b(3)=Jcm(1)*alp(1);
% Link 2 BC
M(4,iB(1))=-1; M(4,iC(1))=1;                 b(4)=m2*aCm(2,1);
M(5,iB(2))=-1; M(5,iC(2))=1;                 b(5)=m2*aCm(2,2)+m2*g;
M(6,iB)=-mc(B,S2); M(6,iC)=mc(C,S2);         b(6)=Jcm(2)*alp(2);
% Link 3 DCE
M(7,iC(1))=-1; M(7,iD(1))=1; M(7,iE(1))=1;   b(7)=m3*aCm(3,1);
M(8,iC(2))=-1; M(8,iD(2))=1; M(8,iE(2))=1;   b(8)=m3*aCm(3,2)+m3*g;
M(9,iC)=-mc(C,S3); M(9,iD)=mc(D,S3); M(9,iE)=mc(E,S3);
                                             b(9)=Jcm(3)*alp(3);
% Link 4 EF
M(10,iE(1))=-1; M(10,iF(1))=1;               b(10)=m4*aCm(4,1);
M(11,iE(2))=-1; M(11,iF(2))=1;               b(11)=m4*aCm(4,2)+m4*g;
M(12,iE)=-mc(E,S4); M(12,iF)=mc(F,S4);       b(12)=Jcm(4)*alp(4);
% Link 5 GFH, artifact reaction Fart applied at H
M(13,iF(1))=-1; M(13,iG(1))=1;               b(13)=m5*aCm(5,1)-Fart(1);
M(14,iF(2))=-1; M(14,iG(2))=1;               b(14)=m5*aCm(5,2)+m5*g-Fart(2);
M(15,iF)=-mc(F,S5); M(15,iG)=mc(G,S5);       b(15)=Jcm(5)*alp(5)-mz(H-S5,Fart);

x = M\b;
end

function Pt = circleIntersect(C1,r1,C2,r2,Pref)
% Branch nearest Pref, which keeps the assembly continuous through the sweep.
d = C2 - C1; L = norm(d);
a = (r1^2 - r2^2 + L^2)/(2*L);
h2 = r1^2 - a^2;
if h2 < 0, h2 = 0; end
h = sqrt(h2);
u = d/L; n = [-u(2) u(1) 0];
Pm = C1 + a*u;
P1 = Pm + h*n; P2 = Pm - h*n;
if norm(P1-Pref) <= norm(P2-Pref), Pt = P1; else, Pt = P2; end
end

function q = rotAbout(p,ctr,ang)
c = cos(ang); s = sin(ang); v = p - ctr;
q = [ctr(1)+c*v(1)-s*v(2), ctr(2)+s*v(1)+c*v(2), 0];
end

function v = rotVec(u,ang)
c = cos(ang); s = sin(ang);
v = [c*u(1)-s*u(2), s*u(1)+c*u(2), 0];
end

function p = perp(r)
% cross([0 0 1], r), the coefficient of an unknown omega or alpha
p = [-r(2) r(1) 0];
end

function a = wrapPi(x)
a = mod(x+pi,2*pi) - pi;
end

tstop = toc(tstart)