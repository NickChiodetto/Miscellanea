% Cascade Speed and Current controller for NGHLS

% by NC @08/06/2020

close all; clear; clc;

%% MACHINE AND ELECTRONICS PARAMETERS

Power = 1e3;

% electro-magnetic
Ld = 21.4e-3;       % d-axis inductance 
Lq = 21.4e-3;       % q-axis inductance
R = 0.6;           % phase resistance (this might change real time due to transient)
phi_PM = 0;  % PM flux (this might change during transient as well)
p = 5;              % pole pairs
rpm = 12000;        % rated speed
ff = rpm*p/60;      % fundamental frequency
wm = 2*pi*50;       % ang elec speed 
fpwm = 20e3;        % switching freq Hz
Ts = 1/fpwm;        % switching period (for now assumed equal to sampling)
Vdc = 30;          % DC-Link voltage

% Fault parameters and power electronics parameters
Load_R = 0.1; %Short circuit test
OC_R = 100000; %Open  circuit test
Cbus = 100e-6; % bus capacitance
Rbus = 0.001; % bus cap ESR

%% CURRENT ERROR SENSOR MODELING

imean = 0;
iaccuracy = 12.8* (40+40)/100;
istdDev = iaccuracy/3;
ivar = istdDev^2; 

%% DC LINK FILTER and SUPPLY MODELING

% circuit parameters
% Power supply has ranges:
% L: 0.6 - 32 microH
% R: 2 - 31 mOhm
supply.V = 540;
supply.R = 3e-2;
supply.L = 30e-6;
% Filter
filter.R = 1e-3;
filter.L = 5e-3;
DCbus.RParallel = 10;
DCbus.CParallel = 40e-6;
% DC bus
DCbus.R = 1e-3;
DCbus.C = 20e-6;

% Extra old computed things for setting natural frequency around point of
% interest
% DC_Link.Req = supply.R + filter.R;
% DC_Link.Leq = supply.L + filter.L;
% 
% natural_f = 1/sqrt(DC_Link.Leq*DCbus.C)/2/pi;
natural_f = 500;

% transient parameters
DC_Link.Trip = 900;
Overspeed_Trip = 1400;
DC_Link.UV_trip = 400;

%% MODELLING 

% general matrices
I = eye(2);         % indentity matrix
J = skewdec(2,0);   % Orthogonal rotation matrix: this is a skew-symmetric matrix
% Resistance Matrices
Rs = R;%*eye(2);
b = [Rs/Ld; 0];
% Inductance Matrices
L = [Ld 0 ;0 Lq];
Linv = inv(L); 
% Permanent magnet flux vector
d = [-1/Ld; 0];
% stator fluxes as state variables (4 state space feedback)
A = [-Rs/Ld wm;
    -wm -Rs/Lq];
C = [1/Ld 0;
    0 1/Lq];

%% CONTROLLERS
% Simmetric optimum, Naslin polynomial
tau_c = 1/fpwm; % considering delay equal to one computational period
t2 = tau_c;
t1 = Ld/R;
ko = 1/R;
alpha = 3;
kr = t1/(alpha^3*ko*t2^2);
tr = alpha^2*t2;

% standard settings for a second order system of PI controllers -----------
ts = 1e-3;
xi = 0.7;
wn = 3/ts/xi;
kp = 2*xi*wn*Ld-R;
ki = wn^2*Ld;

% Defining transfer function (it is the same for both axis as it is jsut passive load)
PI = tf([kp ki],[1 0]);
Plant = tf([1],[Ld R]);
OL = PI*Plant;
CL = feedback(OL,1);
figure
step(CL);

%% plotting some stability stuff
% Bode plot
figure
bode(OL); hold on
bode(CL); hold on
bode(1/(1+OL)); grid on
legend('OL-Plant','CL','sensitivity')
% Nyquist plot
figure 
nyquistplot(OL,'-*'); hold on;

% unit circle
ang = 0:0.01:2*pi; 
xp = 1*cos(ang);    yp = 1*sin(ang);
plot(-1+xp,0+yp,'r');
xlim([-1.5 0.1]);   ylim([-1.5 1.5])
axis equal; hold on; title('Nyquist stability')
legend('Plant+PI'); 

%% SENSORLESS KALMAN FILTER

% noise generator seeds
seed_1 = round(rand*(2^32-1));
seed_2 = round(rand*(2^32-1));
seed_3 = round(rand*(2^32-1));
seed_4 = round(rand*(2^32-1));
seed_5 = round(rand*(2^32-1));
seed_6 = round(rand*(2^32-1));
seed_7 = round(rand*(2^32-1));
seed_8 = round(rand*(2^32-1));
seed_9 = round(rand*(2^32-1));
seed_10 = round(rand*(2^32-1));
seed_11 = round(rand*(2^32-1));

% system noise definition
sys_stdDev_i = 0.2; % currents standard deviation
sys_stdDev_w = 0.2; % speed standard deviation
sys_stdDev_theta = 0.2; % angle standard deviation

% system measurement definition
meas_stdDev_i = sqrt(0.5);

% meas noise covariance matrix
R_cov = (meas_stdDev_i^2)*eye(2); 
% System noise covariance matrix
Q = diag([sys_stdDev_i^2 sys_stdDev_i^2 sys_stdDev_w^2 sys_stdDev_theta^2]);
% Discretization: check the symbolic file in order to see the situation:
% no time dependency because of two reasons:
% 1) linearized state matrix does not depend on angle
% 2) inifinte inertia condition due to the very small sampling period, the
% speed does not really change
Qd = Q*Ts^2;


% %% Simulation
% out = sim('RL_analog_2017.slx');
% 
% figure
% plot(out.id_ref.Time,out.id_ref.Data);
% hold on
% plot(out.id_out.Time,out.id_out.Data);
% xlabel('time [s]'); ylabel('i_q [A]')
% legend('i_d^{ref}','i_d'); grid on

% %% nice plotting
% % Code to create the pdf nice pic from plot -------------------------------
% close(1)
% close(3)
% 
% addpath D:\MATLAB_Important_Scripts\export_fig
% label_size = 30;
% set(gca,'Fontsize',label_size)
% %Printing settings------------------------------
% set(gcf,'PaperUnits','centimeters')
% xSize = 41; ySize = 20;
% xLeft = (41-xSize)/2; yTop = (30-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% 
% % Printing commands for eps and pdf conversion
% print('-depsc','bode_plot'); 
% eps2pdf bode_plot.eps bode_plot.pdf


