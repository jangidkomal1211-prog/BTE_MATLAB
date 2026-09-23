tic
clc;
close all;

%% Physical Parameters
% Grid and geometric parameters
N = 100;                      % Number of grid points
L1 = 110e-9;                % Total length in meters
l = L1/2;                   % Half length
Temp = 500;                 % Temperature in Kelvin

% Material properties
Kn1 = 1;                    % Knudsen number for material 1 (Silicon)
Kn2 = 0.8;                  % Knudsen number for material 2 (Germanium)
thetaD1 = 636;              % Debye temperature for Silicon (K)
thetaD2 = 371;              % Debye temperature for Germanium (K)

% Temperature differences
del_T1 = 30;                % Temperature difference for material 1
del_T2 = -30;               % Temperature difference for material 2

% Lattice constants
a01 = 5.431e-10;            % Silicon lattice constant (m)
a02 = 5.6536e-10;           % Germanium lattice constant (m)

%% Density Calculations
% Silicon properties
density1 = 8/((5.431^3)*1e-30);     % Number of atoms per m³ for Silicon
ndof1 = 3*density1;                 % Total DOFs per m³

% Germanium properties  
density2 = 8/((5.6536^3)*1e-30);    % Number of atoms per m³ for Germanium
ndof2 = 3*density2;                 % Total DOFs per m³

%% Relaxation Time Functions
% Silicon relaxation time components (inverse)
tau1_wLinv = @(w) 7.10*(10^-20)*(10^24)*(w.^2)*Temp*(1-exp(-3*Temp/thetaD1)); %in s
tau1_wTinv = @(w) 9.51*(10^-47)*(10^48)*(w.^4)*Temp*(1-exp(-3*Temp/thetaD1)); %in s
tau1_wL1inv = @(w) 10.9*(10^-20)*(10^24)*(w.^2)*Temp*(1-exp(-3*Temp/thetaD1)); %in s
tau1_wT1inv = @(w) 37.8*(10^-47)*(10^48)*(w.^4)*Temp*(1-exp(-3*Temp/thetaD1)); %in s

% Combined relaxation time for Silicon
tau1_winv = @(w) tau1_wLinv(w) + 2*(tau1_wTinv(w)) + tau1_wL1inv(w) + 2*(tau1_wT1inv(w));

% Germanium relaxation time components (inverse)
tau2_wLinv = @(w) 17.42*(10^-20)*(10^24)*(w.^2)*Temp*(1-exp(-3*Temp/thetaD2)); %in s
tau2_wTinv = @(w) 37.35*(10^-47)*(10^48)*(w.^4)*Temp*(1-exp(-3*Temp/thetaD2)); %in s
tau2_wL1inv = @(w) 19.94*(10^-20)*(10^24)*(w.^2)*Temp*(1-exp(-3*Temp/thetaD2)); %in s
tau2_wT1inv = @(w) 162.06*(10^-47)*(10^48)*(w.^4)*Temp*(1-exp(-3*Temp/thetaD2)); %in s

% Combined relaxation time for Germanium
tau2_winv = @(w) tau2_wLinv(w) + 2*(tau2_wTinv(w)) + tau2_wL1inv(w) + 2*(tau2_wT1inv(w));

% Final relaxation times
tau_w1 = @(w) 1./tau1_winv(w);      % Silicon relaxation time
tau_w2 = @(w) 1./tau2_winv(w);      % Germanium relaxation time

%% Group Velocity Calculations
% GULP gives nu in cm^{-1} 
% GULP gives k in 0:100

% Conversion factors for frequency-k relationship
fac_omega_k1 = (2*pi)*3*1e10*a01*99/pi; % Silicon conversion factor % so that nu is in THz and k is in m^{-1}
fac_omega_k2 = (2*pi)*3*1e10*a02*99/pi;  % Germanium conversion factor

% Silicon dispersion parameters
% Transverse acoustic branch
bT  = 0.570232; % when nu in cm^{-1} is plotted against k in 0:99
cT  = -0.00370167;       
dT  = 2.41779e-05;

% For longitudinal acoustic: Silicon
bL  = 0.229788;
cL  = -0.000172436;        
dL  = 5.26561e-07;

% Average parameters for Silicon
b = (1/fac_omega_k1)*(2*bT+bL)/3;
c = (1/fac_omega_k1)*(2*cT+cL)/3;
d = (1/fac_omega_k1)*(2*dT+dL)/3;

% Silicon group velocity functions
vg1inv = @(w) b + 2*c*w+ 3*d*(w.*w); % in s/m
vg1 = @(w) 1./vg1inv(w);

% Germanium dispersion parameters
% Transverse acoustic branch
bT  = 0.838168; % when nu in cm^{-1} is plotted against k in 0:100
cT  = -0.00755758;       
dT  = 7.18578e-05;

% For longitudinal acoustic: Ge
bL  = 0.361361;
cL  = -0.000108137;        
dL  = 9.64882e-07;

% Average parameters for Germanium
b = (1/fac_omega_k2)*(2*bT+bL)/3;
c = (1/fac_omega_k2)*(2*cT+cL)/3;
d = (1/fac_omega_k2)*(2*dT+dL)/3;

% Germanium group velocity functions
vg2inv = @(w) b + 2*c*w+ 3*d*(w.*w); % in s/m 
vg2 = @(w) 1./vg2inv(w);

% Mean free path calculations
lamda1= @(w) tau_w1(w).*vg1(w); % in m
lamda2= @(w) tau_w2(w).*vg2(w); % in m

%% Density of States Data Loading
fileID = fopen('Si.dens','r');
tline = fgets(fileID);
tline = fgets(fileID);
nlines = 64;
A1 = fscanf(fileID, '%15f %25f' ,[2 nlines]);
x1 = A1(1,:)*0.03*2*pi;
y1 = A1(2,:)/0.03;

fileID = fopen('Ge.dens','r');
tline = fgets(fileID);
tline = fgets(fileID);
nlines = 64;
A2 = fscanf(fileID, '%15f %25f' ,[2 nlines]);
x2 = A2(1,:)*0.03*2*pi;
y2 = A2(2,:)/0.03;

%% Physical Constants
hbar=1.05457e-22; % Reduced Planck constant in J.ps
kB=1.380649e-23; % Boltzmann constant in J/K

%% Maximum Frequencies
wm1 = 235.44*0.03*2*pi;            % Maximum frequency for Silicon
wm2 = 154.73*0.03*2*pi;            % Maximum frequency for Germanium

%% DOS Normalization
% Calculate normalization factors
integral_dos1 = sum(y1)*(x1(2)-x1(1));
integral_dos2 = sum(y2)*(x2(2)-x2(1));

smallF = 1e-28;                     % Small factor for numerical stability
factor_dos1 = (ndof1/integral_dos1)*smallF;
factor_dos2 = (ndof2/integral_dos2)*smallF;

% Apply normalization
y1 = y1*factor_dos1;
y2 = y2*factor_dos2;

% Create interpolation functions
pp1 = csapi(x1,y1); % pp form where x is in THz
pp2 = csapi(x2,y2); % pp form where x is in THz

dos1 = @(w) ppval(pp1,w);          % Silicon DOS function
dos2 = @(w) ppval(pp2,w);          % Germanium DOS function

% Verify normalization
integral_dos1 = integral(dos1,0,x1(end));
integral_dos2 = integral(dos2,0,x2(end));

%% Heat Capacity Functions
% Silicon heat capacity per frequency
C_w1=@(w) (kB*(((hbar*w)/(kB*Temp)).^2)).*ppval(pp1,w).*exp(hbar*w/(kB*Temp)).*((exp(hbar*w/(kB*Temp))-1).^(-2)); %in J/K

% Germanium heat capacity per frequency
C_w2=@(w) (kB*(((hbar*w)/(kB*Temp)).^2)).*ppval(pp2,w).*exp(hbar*w/(kB*Temp)).*((exp(hbar*w/(kB*Temp))-1).^(-2)); %in J/K

%% Transmission and Reflection Coefficients
% Transmission coefficient
t12= @(w) (C_w2(w).*vg2(w)./(C_w1(w).*vg1(w)+C_w2(w).*vg2(w))).*(1-heaviside(w-wm2));

% Reflection coefficient
r12= @(w) 1-t12(w);

fun1 = @(w) C_w1(w)./tau_w1(w);
Dm1 = 1/integral(fun1,0,wm1);
fun2 = @(w) C_w2(w)./tau_w2(w);
Dm2 = 1/integral(fun2,0,wm2);

syms mu;

%% Exponential Integral Functions
% Define exponential integral functions for both materials
E1=@(n,mu,x) (mu.^(n-2)).*exp(-x./(Kn1*mu));
E11= @(x) integral(@(mu) E1(1,mu,x),0,1,'ArrayValued', true);
E12= @(x) integral(@(mu) E1(2,mu,x),0,1,'ArrayValued', true);
E13= @(x) integral(@(mu) E1(3,mu,x),0,1,'ArrayValued', true);

E2=@(n,mu,x) (mu.^(n-2)).*exp(-x./(Kn2*mu));
E21= @(x) integral(@(mu) E2(1,mu,x),0,1);
E22= @(x) integral(@(mu) E2(2,mu,x),0,1);
E23= @(x) integral(@(mu) E2(3,mu,x),0,1);

%% MAIN CALCULATION LOOP - INTERFACE REGIME (Silicon-Ge)

% Set interface emissivity values
e1 = 1; % Emissivity at left interface (Si)
e2 = 1; % Emissivity at right interface (Ge)
a  = l / L1; % Normalized interface position

% Constants of gneralsolutions for each layer
A_1w=  @(w) (1-e1)*vg1(w).*E12(a).*E13(a)-(vg1(w)./r12(w))+(1-e2)*(t12(w)./r12(w)).*vg2(w).*E22(1-a).*E23(1-a);
B_1w= @(w) e1*del_T1-(e1*del_T1*(1-e1)*vg1(w).*E13(a).*E12(a))./(A_1w(w));
B_2w= @(w) -(e2*del_T2*(1-e1)*vg2(w).*E23(1-a).*E12(a))./(A_1w(w));
B_3w= @(w) (1-e1)*(1-((1-e1)*vg1(w).*E13(a).*E12(a))./(A_1w(w)));
B_4w= @(w) -((1-e1)*vg1(w).*E12(a))./(A_1w(w));
B_5w= @(w) -((1-e2)*(1-e1)*vg2(w).*E12(a).*E23(1-a))./(A_1w(w));
B_6w= @(w) -e1*del_T1*vg1(w).*E13(a)./(A_1w(w));
B_7w= @(w) -e2*del_T2*vg2(w).*E23(1-a)./(A_1w(w));
B_8w= @(w) -(1-e1)*vg1(w).*E13(a)./(A_1w(w));
B_9w= @(w) -vg1(w)./A_1w(w);
B_10w= @(w) -(1-e2)*vg2(w).*E23(1-a)./A_1w(w);
B_11w= @(w) -e1*del_T1*((vg1(w).*vg1(w))./vg2(w)).*t12(w).*E13(a)./(A_1w(w).*r12(w));
B_12w= @(w) -e2*del_T2*t12(w).*vg1(w).*E23(1-a)./(A_1w(w).*r12(w));
B_13w= @(w) -(1-e1)*t12(w).*((vg1(w).*vg1(w))./vg2(w)).*E13(a)./(A_1w(w).*r12(w));
B_14w= @(w) -((vg1(w).*vg1(w))./vg2(w)).*t12(w)./(A_1w(w).*r12(w));
B_15w= @(w) -(1-e2)*vg1(w).*t12(w).*E23(1-a)./(A_1w(w).*r12(w));
B_16w= @(w) -e1*del_T1*vg1(w).*t12(w).*E23(1-a).*E13(a)./(A_1w(w).*r12(w));
B_17w= @(w) e2*del_T2*(1-(1-e2)*vg2(w).*t12(w).*E22(1-a).*E23(1-a)./(A_1w(w).*r12(w)));
B_18w= @(w) -(1-e2)*(1-e1)*vg1(w).*t12(w).*E22(1-a).*E13(a)./(A_1w(w).*r12(w));
B_19w= @(w) -(1-e2)*vg1(w).*t12(w).*E23(1-a)./(A_1w(w).*r12(w));
B_20w= @(w) (1-e2)*(1-(1-e2)*vg2(w).*t12(w).*E22(1-a).*E23(1-a)./(A_1w(w).*r12(w)));
B_44w= @(w) -((1-e1)*vg2(w).*E12(a))./(A_1w(w));
B_99w= @(w) -vg2(w)./A_1w(w);
B_144w= @(w) -vg1(w).*t12(w)./(A_1w(w).*r12(w));
B_199w= @(w) -(1-e2)*vg2(w).*t12(w).*E23(1-a)./(A_1w(w).*r12(w));

% --------- Allocate Main Matrices ----------
K=sym(zeros(N+1,N+1));
L=sym(zeros(N+1,N+1));
M=sym(zeros(N+1,N+1));
O=sym(zeros(N+1,N+1));

%Compute Matrix Elements (for [1,1]; loop for others)
fun3= @(w) (C_w1(w).*Kn1./tau_w1(w)).*(2*a/Kn1-1+2*E13(a)+(B_3w(w)+B_8w(w)).*((1/2)-E13(a))*(1-E12(a))+(B_4w(w)+B_9w(w)).*((1/2)-E13(a))^2);
K(1,1)=2*Dm1*integral(fun3,0,wm1);

fun4= @(w) (C_w2(w).*Kn1./tau_w1(w)).*((B_5w(w)+B_10w(w)).*((1/2)-E13(a))*(1-E22(1-a))+(B_44w(w)+B_99w(w)).*((1/2)-E13(a))*((1/2)-E23(1-a)));
L(1,1)=2*Dm1*integral(fun4,0,wm1);

fun5= @(w) (C_w1(w).*Kn2./tau_w2(w)).*((B_13w(w)+B_18w(w)).*((1/2)-E23(1-a))*(1-E12(a))+(B_14w(w)+B_19w(w)).*((1/2)-E13(a))*((1/2)-E23(1-a)));
M(1,1)=2*Dm2*integral(fun5,0,wm2);

fun6= @(w) (C_w2(w).*Kn2./tau_w2(w)).*(2*(1-a)/Kn2-1+2*E23(1-a)+(B_15w(w)+B_20w(w)).*((1/2)-E23(1-a))*(1-E22(1-a))+(B_144w(w)+B_199w(w)).*((1/2)-E23(1-a))^2);
O(1,1)=2*Dm2*integral(fun6,0,wm2);

% Continue for m=2:N+1 and n=2:N+1 using loops
for m=2:N+1

    % % First Column of First Matrix K % %
    A3= @(mu) Kn1*(2*sin((m-1)*a*pi)+(1-(exp(-a./(Kn1*mu)))).*(m-1)*pi.*mu*Kn1.*(-1-cos((m-1)*a*pi)+(m-1)*pi.*mu*Kn1*sin((m-1)*a*pi)))./((m-1)*pi*(1+((m-1)*pi*Kn1*mu).^2));
    A4= integral(A3,0,1);
    A5= @(w) C_w1(w).*A4./(Kn1*tau_w1(w));
    
    A6= @(mu) Kn1*mu.*(1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
    A7= integral(A6,0,1);
    A8= @(w) (C_w1(w)./tau_w1(w)).*(B_3w(w)).*(1-E12(a))*A7;
    
    A9= @(w) (C_w1(w)./tau_w1(w)).*(B_4w(w)).*((1/2)-E13(a))*A7;
    
    A10= @(mu) Kn1*mu.*(-exp(-a./(Kn1*mu))+(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
    A11= integral(A10,0,1);
    A12= @(w) (C_w1(w)./tau_w1(w)).*(B_8w(w)).*(1-E12(a))*A11;
    
    A13= @(w) (C_w1(w)./tau_w1(w)).*(B_9w(w)).*((1/2)-E13(a))*A11;
    
    K(m,1)=2*Dm1*integral(A5,0,wm1)+2*Dm1*integral(A8,0,wm1)+2*Dm1*integral(A9,0,wm1)+2*Dm1*integral(A12,0,wm1)+2*Dm1*integral(A13,0,wm1);
    
    % % First Column of Second Matrix L % %
    A14= @(mu) Kn1*mu.*(1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
    A15= integral(A14,0,1);
    A16= @(w) (C_w2(w)./tau_w1(w)).*(B_5w(w))*(1-E22(1-a))*A15;
    
    A17= @(w) (C_w2(w)./tau_w1(w)).*(B_44w(w))*((1/2)-E23(1-a))*A15;
    
    A18= @(mu) Kn1*mu.*(-exp(-a./(mu*Kn1))+(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
    A19= integral(A18,0,1);
    A20= @(w) (C_w2(w)./tau_w1(w)).*(B_10w(w))*(1-E22(1-a))*A19;
    
    A21= @(w) (C_w2(w)./tau_w1(w)).*(B_99w(w))*((1/2)-E23(1-a))*A19;
    
    L(m,1)=2*Dm1*integral(A16,0,wm1)+2*Dm1*integral(A17,0,wm1)+2*Dm1*integral(A20,0,wm1)+2*Dm1*integral(A21,0,wm1);
    
    % % First Column Third Matrix M % %
    A22= @(mu) Kn2*mu.*(((-1)^m)*exp(-(1-a)./(mu*Kn2))+cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
    A23= integral(A22,0,1);
    A24= @(w) (C_w1(w)./tau_w2(w)).*(B_13w(w)).*(1-E12(a))*A23;
    
    A25= @(w) (C_w1(w)./tau_w2(w)).*(B_14w(w)).*((1/2)-E13(a))*A23;
    
    A26= @(mu) Kn2*mu.*(((-1)^(m-1))-exp(-(1-a)./(mu*Kn2)).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
    A27= integral(A26,0,1);
    A28= @(w) (C_w1(w)./tau_w2(w)).*(B_18w(w)).*(1-E12(a))*A27;
    
    A29= @(w) (C_w1(w)./tau_w2(w)).*(B_19w(w)).*((1/2)-E13(a))*A27;
    
    M(m,1)=2*Dm2*integral(A24,0,wm2)+2*Dm2*integral(A25,0,wm2)+2*Dm2*integral(A28,0,wm2)+2*Dm2*integral(A29,0,wm2);
    
    % % First Column of Forth Matrix O % %
    A32= @(mu) Kn2.*(2*(1./((m-1)*pi)).*(sin((m-1)*pi)-sin((m-1)*a*pi))+(1-exp((a-1)./(mu*Kn2))).*mu*Kn2.*(-cos((m-1)*pi)-cos((m-1)*a*pi)+(m-1)*pi.*mu.*(sin((m-1)*pi)-sin((m-1)*a*pi))*Kn2))./(1+((m-1)*pi.*mu.*Kn2).^2);
    A33= integral(A32,0,1);
    A34= @(w) C_w2(w).*A33./(Kn2*tau_w2(w));
    
    A35= @(mu) Kn2*mu.*(((-1)^(m))*(exp(-(1-a)./(Kn2*mu)))+cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
    A36= integral(A35,0,1);
    A37= @(w) (C_w2(w)./tau_w2(w)).*(B_15w(w)).*(1-E22(1-a))*A36;
    
    A38= @(w) (C_w2(w)./tau_w2(w)).*(B_144w(w)).*((1/2)-E23(1-a))*A36;
    
    A39= @(mu) Kn2*mu.*(((-1)^(m-1))-(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
    A40= integral(A39,0,1);
    A41= @(w) (C_w2(w)./tau_w2(w)).*(B_20w(w)).*(1-E22(1-a))*A40;
    
    A42= @(w) (C_w2(w)./tau_w2(w)).*(B_199w(w)).*((1/2)-E23(1-a))*A40;
    
    O(m,1)=2*Dm2*integral(A34,0,wm2)+2*Dm2*integral(A37,0,wm2)+2*Dm2*integral(A38,0,wm2)+2*Dm2*integral(A41,0,wm2)+2*Dm2*integral(A42,0,wm2);
    
    for n=2:N+1
    % % First Row of First Matrix K % %
    B3= @(mu) 2*cos(a*pi*(n-1)/2)*Kn1.*(2*sin((n-1)*a*pi/2)+(1-(exp(-a./(Kn1*mu)))).*(n-1)*pi.*mu*Kn1.*(-cos((n-1)*a*pi/2)+(n-1)*pi.*mu*Kn1*sin((n-1)*a*pi/2)))./((n-1)*pi*(1+((n-1)*pi*Kn1.*mu).^2));
    B4= integral(B3,0,1);
    B5= @(w) C_w1(w).*B4./(Kn1*tau_w1(w));
    
    B6= @(mu) (1-(exp(-a./(Kn1*mu))).*(cos((n-1)*a*pi)-(n-1)*pi*mu*Kn1*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn1*mu).^2);
    B7= integral(B6,0,1);
    B8= @(w) (C_w1(w)*Kn1./tau_w1(w)).*(B_3w(w)+B_8w(w)).*((1/2)-E13(a))*B7;
    
    B9= @(mu) mu.*(-(exp(-a./(Kn1*mu)))+cos((n-1)*a*pi)+(n-1)*pi*mu*Kn1*sin((n-1)*a*pi))./(1+(((n-1)*pi)^2).*(Kn1*mu).^2);
    B10= integral(B9,0,1);
    B11= @(w) (C_w1(w)*Kn1./tau_w1(w)).*(B_4w(w)+B_9w(w)).*((1/2)-E13(a))*B10;
    
    K(1,n)=2*Dm1*integral(B5,0,wm1)+2*Dm1*integral(B8,0,wm1)+2*Dm1*integral(B11,0,wm1);
    
    % % First Row of Second Matrix L % %
    B12= @(mu) ((-1)^(n-1)-exp(-(1-a)./(mu*Kn2)).*(cos((n-1)*a*pi)+(n-1)*pi*mu*Kn2*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
    B13= integral(B12,0,1);
    B14= @(w) (C_w2(w)*Kn1./tau_w1(w)).*(B_5w(w)+B_10w(w)).*((1/2)-E13(a))*B13;
    
    B15= @(mu) mu.*(((-1)^n)*exp(-(1-a)./(mu*Kn2))+cos((n-1)*a*pi)-(n-1)*pi*mu*Kn2*sin((n-1)*a*pi))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
    B16= integral(B15,0,1);
    B17= @(w) (C_w2(w)*Kn1./tau_w1(w)).*(B_44w(w)+B_99w(w)).*((1/2)-E13(a))*B16;
    
    L(1,n)=2*Dm1*integral(B14,0,wm1)+2*Dm1*integral(B17,0,wm1);
    
    % % First Row of Third Matrix M % %
    B18= @(mu) (1-exp(-a./(mu*Kn1)).*(cos((n-1)*a*pi)-(n-1)*pi*mu*Kn1*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn1*mu).^2);
    B19= integral(B18,0,1);
    B20= @(w) (C_w1(w)*Kn2./tau_w2(w)).*(B_13w(w)+B_18w(w)).*((1/2)-E23(1-a))*B19;
    
    B21= @(mu) mu.*(-exp(-a./(mu*Kn1))+cos((n-1)*a*pi)+(n-1)*pi*mu*Kn1*sin((n-1)*a*pi))./(1+(((n-1)*pi)^2).*(Kn1*mu).^2);
    B22= integral(B21,0,1);
    B23= @(w) (C_w1(w)*Kn2./tau_w2(w)).*(B_14w(w)+B_19w(w)).*((1/2)-E23(1-a))*B22;
    
    M(1,n)=2*Dm2*integral(B20,0,wm2)+2*Dm2*integral(B23,0,wm2);
    
    % % First Row of Forth Matrix O % %
    B36= @(mu) 2*cos((1+a)*pi*(n-1)/2)*Kn2.*(-2*sin((n-1).*(-1+a)*pi/2)+(-1+(exp((-1+a)./(Kn2*mu)))).*(n-1)*pi.*mu*Kn2.*(cos((n-1).*(-1+a)*pi/2)+(n-1)*pi.*mu*Kn2*sin((n-1).*(-1+a)*pi/2)))./((n-1)*pi.*(1+((n-1)*pi*Kn2.*mu).^2));
    B37= integral(B36,0,1);
    B38= @(w) C_w2(w).*B37./(Kn2*tau_w2(w));
    
    B39= @(mu) ((-1)^(n-1)-(exp(-(1-a)./(Kn2*mu))).*(cos((n-1)*a*pi)+(n-1)*pi*mu*Kn2*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
    B40= integral(B39,0,1);
    B41= @(w) (C_w2(w)*Kn2./tau_w2(w)).*(B_15w(w)+B_20w(w)).*((1/2)-E23(1-a))*B40;
    
    B42= @(mu) mu.*(((-1)^n)*(exp(-(1-a)./(Kn2*mu)))+cos((n-1)*a*pi)-(n-1)*pi*mu*Kn2*sin((n-1)*a*pi))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
    B43= integral(B42,0,1);
    B44= @(w) (C_w2(w)*Kn2./tau_w2(w)).*(B_144w(w)+B_199w(w)).*((1/2)-E23(1-a))*B43;
    
    O(1,n)=2*Dm2*integral(B38,0,wm2)+2*Dm2*integral(B41,0,wm2)+2*Dm2*integral(B44,0,wm2);
     
    if m==n
        % % Rest Elements of First Matrix K % %
        C3= @(mu) Kn1.*(4*(m-1)*pi.*mu*Kn1.*(exp(-a./(mu*Kn1))).*(cos((m-1)*a*pi)-(m-1)*pi.*mu.*sin((m-1)*a*pi).*Kn1)+(2*a*(m-1)*pi+sin(2*(m-1)*a*pi)+(m-1)*pi.*mu*Kn1.*(-3-cos(2*(m-1)*a*pi)+(m-1)*pi.*mu*Kn1.*(2*(m-1)*a*pi+sin(2*(m-1)*a*pi)+2*(m-1)*pi.*mu*Kn1.*(sin((m-1)*a*pi)).^2))))./(2*(m-1)*pi.*(1+((m-1)*pi.*mu*Kn1).^2).^2);
        C4= integral(C3,0,1);
        C5= @(w) C_w1(w).*C4./(Kn1*tau_w1(w));
    
        C6= @(mu) (1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        C7= @(mu) Kn1*mu.*C6(mu);
        C8= integral(C6,0,1);
        C9= integral(C7,0,1);
        C10= @(w) (C_w1(w)./tau_w1(w)).*(B_3w(w))*C8*C9;
    
        C11= @(mu) (-(exp(-a./(Kn1*mu)))+cos((m-1)*a*pi)+(m-1)*pi*mu*Kn1*sin((m-1)*a*pi))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        C12= @(mu) mu.*C11(mu);
        C13= integral(C11,0,1);
        C14= integral(C12,0,1);
        C15= @(w) (C_w1(w)./tau_w1(w)).*(B_4w(w))*C9*C14;
    
        C16= @(w) (C_w1(w)./tau_w1(w)).*(B_8w(w))*C9*C13;
    
        C17= @(mu) Kn1*mu.*C11(mu);
        C18= integral(C17,0,1);
        C19= @(w) (C_w1(w)./tau_w1(w)).*(B_9w(w))*C18*C14;
    
        K(m,n)=2*Dm1*integral(C5,0,wm1)+2*Dm1*integral(C10,0,wm1)+2*Dm1*integral(C15,0,wm1)+2*Dm1*integral(C16,0,wm1)+2*Dm1*integral(C19,0,wm1);
    
        % % Rest Elements of Second Matrix L % %
        C20= @(mu) Kn1*mu.*(1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        C21= @(mu) (((-1)^(m-1))-(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C22= integral(C20,0,1);
        C23= integral(C21,0,1);
        C24= @(w) (C_w2(w)./tau_w1(w)).*(B_5w(w))*C22*C23;
    
        C26= @(mu) mu.*(((-1)^(m))*(exp(-(1-a)./(Kn2*mu)))+cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C28= integral(C26,0,1);
        C29= @(w) (C_w2(w)./tau_w1(w)).*(B_44w(w))*C22*C28;
    
        C30= @(mu) Kn1*mu.*(-(exp(-(a)./(Kn1*mu)))+(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        C31= integral(C30,0,1);
        C32= @(w) (C_w2(w)./tau_w1(w)).*(B_10w(w))*C31*C23;
    
        C33= @(mu) Kn1*mu.*(-exp(-a./(Kn1*mu))+(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        C34= integral(C33,0,1);
        C35= @(w) (C_w2(w)./tau_w1(w)).*(B_99w(w))*C34*C28;
    
        L(m,n)=2*Dm1*integral(C24,0,wm1)+2*Dm1*integral(C29,0,wm1)+2*Dm1*integral(C32,0,wm1)+2*Dm1*integral(C35,0,wm1);
    
        % % Rest Elements of Third Matrix M % %
        C36= @(mu) Kn2*mu.*(((-1)^(m)).*(exp(-(1-a)./(Kn2*mu)))+cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C37= @(mu) (1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        C38= integral(C36,0,1);
        C39= integral(C37,0,1);
        C40= @(w) (C_w1(w)./tau_w2(w)).*(B_13w(w))*C38*C39;
    
        C41= @(mu) mu.*(-(exp(-a./(Kn1*mu)))+cos((m-1)*a*pi)+(m-1)*pi*mu*Kn1*sin((m-1)*a*pi))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        C42= integral(C41,0,1);
        C43= @(w) (C_w1(w)./tau_w2(w)).*(B_14w(w))*C42*C38;
    
        C44= @(mu) Kn2*mu.*(((-1)^(m-1))-(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C45= integral(C44,0,1);
        C46= @(w) (C_w2(w)./tau_w1(w)).*(B_18w(w))*C45*C39;
    
        C47= @(w) (C_w2(w)./tau_w1(w)).*(B_19w(w))*C42*C45;
    
        M(m,n)=2*Dm2*integral(C40,0,wm2)+2*Dm2*integral(C43,0,wm2)+2*Dm2*integral(C46,0,wm2)+2*Dm2*integral(C47,0,wm2);
    
        % % Rest Elements of Forth Matrix O % %
        C50= @(mu) Kn2.*(-4*exp((-1+a)./(mu*Kn2)).*(m-1)*pi.*mu*Kn2.*(-cos((m-1)*pi)+(m-1)*pi*mu.*sin((m-1)*pi).*Kn2).*(cos((m-1)*a*pi)+(m-1)*pi.*mu.*sin((m-1)*a*pi)*Kn2)+(-2*(-1+a)*(m-1)*pi+sin(2*(m-1)*pi)-sin(2*a*(m-1)*pi)+(m-1)*pi.*mu*Kn2.*(-2-cos(2*(m-1)*pi)-cos(2*a*(m-1)*pi)+(m-1)*pi.*mu*Kn2.*(-2*(-1+a)*(m-1)*pi+sin(2*(m-1)*pi)-sin(2*a*(m-1)*pi)+2*(m-1)*pi.*mu.*((sin((m-1)*pi).^2)+(sin(a*(m-1)*pi)).^2)*Kn2))))./(2*(m-1)*pi*(1+(pi.*mu*Kn2*(m-1)).^2).^2);
        C51= integral(C50,0,1);
        C52= @(w) C_w2(w).*C51./(Kn2*tau_w2(w));
    
        C53= @(mu) Kn2*mu.*(-(exp(-(1-a)./(Kn2*mu)))+((-1)^(m-1))*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C54= @(mu) (1+((-1)^m)*(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C55= integral(C53,0,1);
        C56= integral(C54,0,1);
        C57= @(w) (C_w2(w)./tau_w2(w)).*(B_15w(w))*C55*C56;
    
        C58= @(mu) mu.*(((-1)^(m-1))*(exp(-(1-a)./(Kn2*mu)))-(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C59= @(mu) Kn2*C58(mu);
        C60= integral(C58,0,1);
        C61= integral(C59,0,1);
        C62= @(w) (C_w2(w)./tau_w2(w)).*(B_144w(w))*C60*C61;
    
        C63= @(mu) (((-1)^(m-1))-(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C64= @(mu) Kn2*mu.*C63(mu);
        C65= integral(C63,0,1);
        C66= integral(C64,0,1);
        C67= @(w) (C_w2(w)./tau_w2(w)).*(B_20w(w))*C65*C66;
    
        C68= @(mu) mu.*(((-1)^m)+(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        C69= integral(C68,0,1);
        C70= @(w) (C_w2(w)./tau_w2(w)).*(B_199w(w))*C69*C61;
    
        O(m,n)=2*Dm2*integral(C52,0,wm2)+2*Dm2*integral(C57,0,wm2)+2*Dm2*integral(C62,0,wm2)+2*Dm2*integral(C67,0,wm2)+2*Dm2*integral(C70,0,wm2);
    
    else
        % % First Matrix % %
        D3= @(mu) Kn1.*(-(exp(-a./(mu*Kn1))).*(m-n).*(m+n-2)*pi.*mu*Kn1.*(-cos((m-1)*a*pi)-cos((n-1)*a*pi)+pi.*mu.*((m-1)*sin((m-1)*a*pi)+(n-1).*sin((n-1)*a*pi))*Kn1)+(2*(m-1)*cos((n-1)*a*pi).*sin((m-1)*a*pi)-2*(n-1)*cos((m-1)*a*pi).*sin((n-1)*a*pi)+pi.*mu*Kn1.*(-(m-n)*(m+n-2).*(1+cos((m-1)*a*pi).*cos((n-1)*a*pi))+pi.*mu*Kn1.*(((m-1)^2+(n-1)^2).*((m-1)*cos((n-1)*a*pi).*sin((m-1)*a*pi)-(n-1)*cos((m-1)*a*pi).*sin((n-1)*a*pi))+(m-1)*(m-n)*(n-1)*(m+n-2)*pi.*mu*sin((m-1)*a*pi).*sin((n-1)*a*pi)*Kn1))))./((m-n)*(m+n-2)*pi*(1+((m-1)^2).*mu.*mu*Kn1*Kn1*pi*pi).*(1+((n-1)*pi.*mu*Kn1).^2));
        D4= integral(D3,0,1);
        D5= @(w) C_w1(w).*D4./(Kn1*tau_w1(w));
    
        D6= @(mu) Kn1*mu.*(1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        D7= @(mu) (1-(exp(-a./(Kn1*mu))).*(cos((n-1)*a*pi)-(n-1)*pi*mu*Kn1*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn1*mu).^2);
        D8= integral(D6,0,1);
        D9= integral(D7,0,1);
        D10= @(w) (C_w1(w)./tau_w1(w)).*(B_3w(w))*D8*D9;
    
        D11= @(mu) mu.*(-(exp(-a./(Kn1*mu)))+cos((n-1)*a*pi)+(n-1)*pi*mu*Kn1*sin((n-1)*a*pi))./(1+(((n-1)*pi)^2).*(Kn1*mu).^2);
        D12= integral(D11,0,1);
        D13= @(w) (C_w1(w)./tau_w1(w)).*(B_4w(w))*D12*D8;
    
        D14= @(mu) Kn1*mu.*(-(exp(-a./(Kn1*mu)))+cos((m-1)*a*pi)+(m-1)*pi*mu*Kn1*sin((m-1)*a*pi))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        D15= integral(D14,0,1);
        D16= @(w) (C_w1(w)./tau_w1(w)).*(B_8w(w))*D15*D9;
    
        D17= @(w) (C_w1(w)./tau_w1(w)).*(B_9w(w))*D15*D12;
    
        K(m,n)=2*Dm1*integral(D5,0,wm1)+2*Dm1*integral(D10,0,wm1)+2*Dm1*integral(D13,0,wm1)+2*Dm1*integral(D16,0,wm1)+2*Dm1*integral(D17,0,wm1);
    
        % % Second Matrix % %
        D18= @(mu) Kn1*mu.*(1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        D19= @(mu) (((-1)^(n-1))-(exp(-(1-a)./(Kn2*mu))).*(cos((n-1)*a*pi)+(n-1)*pi*mu*Kn2*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
        D20= integral(D18,0,1);
        D21= integral(D19,0,1);
        D22= @(w) (C_w2(w)./tau_w1(w)).*(B_5w(w))*D20*D21;
    
        D24= @(mu) mu.*(((-1)^(n))*(exp(-(1-a)./(Kn2*mu)))+(cos((n-1)*a*pi)-(n-1)*pi*mu*Kn2*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
        D26= integral(D24,0,1);
        D27= @(w) (C_w2(w)./tau_w1(w)).*(B_44w(w))*D20*D26;
    
        D28= @(mu) Kn1*mu.*(-(exp(-(a)./(Kn1*mu)))+(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn1*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
        D30= integral(D28,0,1);
        D32= @(w) (C_w2(w)./tau_w1(w)).*(B_10w(w))*D30*D21;
    
        D35= @(w) (C_w2(w)./tau_w1(w)).*(B_99w(w))*D30*D26;
    
        L(m,n)=2*Dm1*integral(D22,0,wm1)+2*Dm1*integral(D27,0,wm1)+2*Dm1*integral(D32,0,wm1)+2*Dm1*integral(D35,0,wm1);
    
        % % Third Matrix % %
        D36= @(mu) Kn2*mu.*(((-1)^(m))*(exp(-(1-a)./(Kn2*mu)))+cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        D37= @(mu) (1-(exp(-a./(Kn1*mu))).*(cos((n-1)*a*pi)-(n-1)*pi*mu*Kn1*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn1*mu).^2);
        D38= integral(D36,0,1);
        D39= integral(D37,0,1);
        D40= @(w) (C_w1(w)./tau_w2(w)).*(B_13w(w))*D38*D39;
    
        D41= @(mu) mu.*(-(exp(-a./(Kn1*mu)))+cos((n-1)*a*pi)+(n-1)*pi*mu*Kn1*sin((n-1)*a*pi))./(1+(((n-1)*pi)^2).*(Kn1*mu).^2);
        D42= integral(D41,0,1);
        D43= @(w) (C_w1(w)./tau_w2(w)).*(B_14w(w))*D42*D38;
    
        D44= @(mu) Kn2*mu.*(((-1)^(m-1))-(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        D45= integral(D44,0,1);
        D46= @(w) (C_w1(w)./tau_w2(w)).*(B_18w(w))*D45*D39;
    
        D47= @(w) (C_w2(w)./tau_w1(w)).*(B_19w(w))*D42*D45;
    
        M(m,n)=2*Dm2*integral(D40,0,wm2)+2*Dm2*integral(D43,0,wm2)+2*Dm2*integral(D46,0,wm2)+2*Dm2*integral(D47,0,wm2);
    
        % % Forth Matrix % %
        D50= @(mu) Kn2.*(2*((m-1).*cos((n-1)*pi).*sin((m-1)*pi)-(m-1).*cos((n-1)*a*pi).*sin((m-1)*a*pi)-(n-1)*cos((m-1)*pi).*sin((n-1)*pi)+(n-1)*cos((m-1)*a*pi).*sin((n-1)*a*pi))+pi.*mu*Kn2.*((m-n)*(m+n-2)*(exp((-1+a)./(mu*Kn2)).*(cos((m-1)*a*pi).*cos((n-1)*pi)+cos((m-1)*pi).*cos((n-1)*a*pi))-(cos((m-1)*pi).*cos((n-1)*pi)+cos((m-1)*a*pi).*cos((n-1)*a*pi)))+pi.*mu*Kn2.*(-exp((-1+a)./(mu*Kn2)).*(m-n)*(m+n-2).*((m-1)*cos((n-1)*a*pi).*sin((m-1)*pi)-(m-1)*cos((n-1)*pi).*sin((m-1)*a*pi)+(n-1)*cos((m-1)*a*pi).*sin((n-1)*pi)-(n-1)*cos((m-1)*pi).*sin((n-1)*a*pi)+(m-1)*(n-1)*pi.*mu.*(sin((m-1)*a*pi).*sin((n-1)*pi)+sin((m-1)*pi).*sin((n-1)*a*pi))*Kn2)+(((m-1)^2+(n-1)^2)*((m-1)*cos((n-1)*pi).*sin((m-1)*pi)-(m-1)*cos((n-1)*a*pi).*sin((m-1)*a*pi)-(n-1)*cos((m-1)*pi).*sin((n-1)*pi)+(n-1)*cos((m-1)*a*pi).*sin((n-1)*a*pi))+(m-1)*(m-n)*(n-1)*(m+n-2)*pi.*mu*(sin((m-1)*pi).*sin((n-1)*pi)+sin((m-1)*a*pi).*sin((n-1)*a*pi))*Kn2))))./((m-n)*(m+n-2)*pi*(1+((m-1)*pi.*mu*Kn2).^2).*(1+((n-1)*pi.*mu*Kn2).^2));
        D51= integral(D50,0,1);
        D52= @(w) C_w2(w).*D51./(Kn2*tau_w2(w));
    
        D53= @(mu) Kn2*mu.*(((-1)^(m)).*(exp(-(1-a)./(Kn2*mu)))+(cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        D54= @(mu) (((-1)^(n-1))-(exp(-(1-a)./(Kn2*mu))).*(cos((n-1)*a*pi)+(n-1)*pi*mu*Kn2*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
        D55= integral(D53,0,1);
        D56= integral(D54,0,1);
        D57= @(w) (C_w2(w)./tau_w2(w)).*(B_15w(w))*D55*D56;
    
        D58= @(mu) mu.*(((-1)^(n))*(exp(-(1-a)./(Kn2*mu)))+(cos((n-1)*a*pi)-(n-1)*pi*mu*Kn2*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
        D59= integral(D58,0,1);
        D60= @(w) (C_w2(w)./tau_w2(w)).*(B_144w(w))*D59*D55;
    
        D61= @(mu) Kn2*mu.*(((-1)^(m))+(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        D62= @(mu) (((-1)^(n))+(exp(-(1-a)./(Kn2*mu))).*(cos((n-1)*a*pi)+(n-1)*pi*mu*Kn2*sin((n-1)*a*pi)))./(1+(((n-1)*pi)^2).*(Kn2*mu).^2);
        D63= integral(D61,0,1);
        D64= integral(D62,0,1);
        D65= @(w) (C_w2(w)./tau_w2(w)).*(B_20w(w))*D63*D64;
    
        D66= @(mu) Kn2*mu.*(((-1)^(m-1))-(exp(-(1-a)./(Kn2*mu))).*(cos((m-1)*a*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
        D67= integral(D66,0,1);
        D68= @(w) (C_w2(w)./tau_w2(w)).*(B_199w(w))*D67*D59;
    
        O(m,n)=2*Dm2*integral(D52,0,wm2)+2*Dm2*integral(D57,0,wm2)+2*Dm2*integral(D60,0,wm2)+2*Dm2*integral(D65,0,wm2)+2*Dm2*integral(D68,0,wm2);
    
    end
    end
end

% Global Matrix Containing matrices K, L, M, and O
A=sym(zeros(2*N+2,2*N+2));

A(1,1)= 1-K(1,1)/4;
A(1,2)=-L(1,1)/4;
A(2,1)=-M(1,1)/4;
A(2,2)=1-O(1,1)/4;

for m=2:N+1
    A(2*m-1,1)=-K(m,1)/2;
    A(2*m-1,2)=-L(m,1)/2;
    A(2*m,1)=-M(m,1)/2;
    A(2*m,2)=-O(m,1)/2;
    for n=2:N+1
        A(1,2*n-1)=-K(1,n)/4;
        A(1,2*n)=-L(1,n)/4;
        A(2,2*n-1)=-M(1,n)/4;
        A(2,2*n)=-O(1,n)/4;

        if m==n
            A(2*m-1,2*n-1)=1-K(m,n)/2;
            A(2*m-1,2*n)=-L(m,n)/2;
            A(2*n,2*m-1)=-M(m,n)/2;
            A(2*m,2*n)=1-O(m,n)/2;
        else
            A(2*m-1,2*n-1)=-K(m,n)/2;
            A(2*m-1,2*n)=-L(m,n)/2;
            A(2*m,2*n-1)=-M(m,n)/2;
            A(2*m,2*n)=-O(m,n)/2;
        end
     end
end

% --------- Build Forcing Vectors ----------
F=sym(zeros(2*N+2,1));
F1=sym(zeros(N+1,1));
F2=sym(zeros(N+1,1));

fun7= @(w) (C_w1(w)*Kn1./tau_w1(w)).*((B_1w(w)+B_6w(w)).*((1/2)-E13(a)))+(C_w2(w)*Kn1./tau_w1(w)).*((B_2w(w)+B_7w(w)).*((1/2)-E13(a)));
F1(1)=(Dm1/2)*integral(fun7,0,wm1);

fun8= @(w) (C_w1(w)*Kn2./tau_w2(w)).*((B_11w(w)+B_16w(w)).*((1/2)-E23(1-a)))+(C_w2(w)*Kn2./tau_w2(w)).*((B_12w(w)+B_17w(w)).*((1/2)-E23(1-a)));
F2(1)=(Dm2/2)*integral(fun8,0,wm2);

for m=2:N+1
    % % Elements of First Matrix F1 (First Material) % %
    fun11= @(w,mu) (C_w1(w)./tau_w1(w)).*Kn1.*mu.*(B_1w(w).*(1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi.*mu*Kn1.*sin((m-1)*a*pi)))+B_6w(w).*(-(exp(-a./(Kn1*mu)))+cos((m-1)*a*pi)+(m-1)*pi.*mu*Kn1.*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
    fun12= @(w,mu) (C_w2(w)./tau_w1(w)).*Kn1.*mu.*(B_2w(w).*(1-(exp(-a./(Kn1*mu))).*(cos((m-1)*a*pi)-(m-1)*pi.*mu*Kn1.*sin((m-1)*a*pi)))+B_7w(w).*(-(exp(-a./(Kn1*mu)))+cos((m-1)*a*pi)+(m-1)*pi.*mu*Kn1.*sin((m-1)*a*pi)))./(1+(((m-1)*pi)^2).*(Kn1*mu).^2);
    fun13= @(w,mu) fun11(w,mu)+fun12(w,mu);
    fun14= @(w) integral(@(mu) fun13(w,mu),0,1,'ArrayValued', true);
    F1(m)=Dm1*integral(fun14,0,wm1);
    
    % % Elements of Second Matrix F1 (Second Material) % %
    fun15= @(w,mu) (C_w1(w)./tau_w2(w))*Kn2*mu.*(B_11w(w).*(cos((m-1)*a*pi)-(m-1)*pi*mu*sin((m-1)*a*pi)*Kn2+(exp(-(1-a)./(Kn2*mu)))*(-cos((m-1)*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*pi)))+B_16w(w).*(cos((m-1)*pi)+(m-1)*pi*mu*sin((m-1)*pi)*Kn2+(exp(-(1-a)./(Kn2*mu)))*(-cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi))))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
    fun16= @(w,mu) (C_w2(w)./tau_w2(w))*Kn2*mu.*(B_12w(w).*(cos((m-1)*a*pi)-(m-1)*pi*mu*sin((m-1)*a*pi)*Kn2+(exp(-(1-a)./(Kn2*mu)))*(-cos((m-1)*pi)+(m-1)*pi*mu*Kn2*sin((m-1)*pi)))+B_17w(w).*(cos((m-1)*pi)+(m-1)*pi*mu*sin((m-1)*pi)*Kn2+(exp(-(1-a)./(Kn2*mu)))*(-cos((m-1)*a*pi)-(m-1)*pi*mu*Kn2*sin((m-1)*a*pi))))./(1+(((m-1)*pi)^2).*(Kn2*mu).^2);
    fun17= @(w,mu) fun15(w,mu)+fun16(w,mu);
    fun18= @(w) integral(@(mu) fun17(w,mu),0,1,'ArrayValued', true);
    F2(m)=Dm2*integral(fun18,0,wm2);
end

% Sorting Elements of Global Force Vector F from Vectors F1 and F2
for j=1:N+1
    F(2*j-1)=F1(j);
    F(2*j)=F2(j);
end

% Solve SYSTEM for Temperature Coefficients which gives Fourier coefficients of Temperature field
syms x y;
X=A\F;

C=sym(zeros(N+1,1));
X1=zeros(1,N+1);
X2=zeros(1,N+1);

for m=1:N+1
    C(m)=cos((m-1)*pi*x);
    X1(m)=X(2*m-1);
    X2(m)=X(2*m);
end

% Final Temperature Values after substituting Fourier coefficients in Temperature Fourier Series
del_TT1=X1*C;
del_TT2=X2*C;
 

% Space domain discretization
p=0:0.01:1;

final_del_TT1=piecewise((x>=0) & (x<=a),del_TT1,(x>=a) & (x<=1),del_TT2);
final_delT1=vpa(subs(final_del_TT1,x,p))+Temp*ones(1,size(final_del_TT1,2));

data3 = [p(:), final_delT1(:)];  % Ensure both are column vectors
% Write to file
fid = fopen('Temperature', 'w');
for i = 1:length(p)
    fprintf(fid, '%.6f\t%.6f\n', L1*1e9*data3(i, 1), data3(i, 2));  % Tab-separated
end
fclose(fid);

figure;
plot(p*L1*1e9, final_delT1, '-', 'MarkerFaceColor','k', 'Color', 'k', 'MarkerSize', 7, 'LineWidth', 3);
xlabel('x(nm)', 'Interpreter', 'latex');
ylabel('T(K)', 'Interpreter', 'latex');
xlim([0 110]);
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
% Get current figure handle
f = gcf;
% Export the figure to PNG and PDF formats
exportgraphics(f, 'Temperature.png', 'Resolution', 150);
exportgraphics(f, 'Temperature.pdf', 'ContentType', 'vector');
hold off;

toc
