clc
clear all
close all
tic

%% user inputs (non dim)

sigma = 5643;              % dimensionless scan rate
dTheta = 0.194;            % dimensionless potential
r1 = 950;               % dimensionless rate constant for Step 1
r2 = 120;               % dimensionless rate constant for Step 2

epsilon = 1E-01;           % Electrode radius (cm)
D = 6.9E-06;               % Diffusivity (cm^2 sec^-1)
R = 8.314;                 % Universal gas constant (J mol^-1 K^-1)
F = 96485.332;             % Faraday's constant (C mol^-1)
T = 298;                   % Standard temperature (K)

v = (sigma*R*T*D)/(epsilon^2 * F);      %scan rate (V s^-1)
dE = (dTheta*R*T)/F;                    %Potential (V)

k01 = (r1*D)/(epsilon^2);               %Apparent rate constant for Step A
k02 = (r2*D)/(epsilon);                 %Apparent rate constant for Step B


%% dimensional solver

Cba = 1E-06;              %Bulk concentration of FT (mol cm^-3)
alfa1 = 0.95;              %Cathodic transfer coeff for Step A
alfaP1 = 1-alfa1;
alfa2 = 0.775;              %Cathodic transfer coeff for Step B
alfaP2 = 1-alfa2;        
Area = 1;                 %Electode area (cm^2)   

EL = -0.8;                %Voltametry starting voltage
ER = 0.5;                 %Voltametry ending voltage
EM = -0.2;
Ef01 = -0.61;            %Formal electrode potential for Step A
E0rds = 0.005;            %Formal potential for second electron trasnfer of step B 
E0nonrds = 0.001;         %Formal potential for first electron trasnfer of step B 

tmax = 2*(abs(EL-ER)/v);  %Maximum duration of simulation (s)
dt = dE/v;                %Duration of each time step (s)   

N = 500;                  %Number of grid points
L = 0.2*sqrt(D*tmax);     %Span of simulation domain (cm)
dx = L/N;                 %Space steps

t = dt:dt:tmax;           %Time vector
C0A = Cba*ones(1,N);      %Initial concentraion field

aP0 = dx/dt;              %Mass diffusion coeff for time neighbour
aE = D/dx;                %Mass diffusion for east face
aW = D/dx;                %Mass diffusion for west face

E(1) = EM;                %Initialization of potential vector
CA = C0A';                %Initialization of concentration vector (species A)
CB = (zeros(1,N))';       %Initialization of concentration vector (species B)
CC = (zeros(1,N))';       %Initialization of concentration vector (species B)

%% Voltage vector evaluation for transient boundary condition
for m=1:length(t)
    if m < (0.25*length(t))
    E(m+1) = E(m) - dE;
    elseif m < (0.75*length(t))
    E(m+1) = E(m) + dE;
    else
    E(m) = E(m-1) - dE;
    end
end

%% Main time loop
for m=1:length(t)
    
    % coefficient matrix
    %for specie A
    MatA(1,1) = aE + aP0 + k01*(exp(((-(1+alfa1)*F)/(R*T))*(E(m)-Ef01)));
    MatA(1,2) = -aE;
    for i = 2:N-1
    MatA(i,i) = aE + aW + aP0;
    MatA(i,i-1) = -aW;
    MatA(i,i+1) = -aE;
    end
    MatA(N,N-1) = -aW;
    MatA(N,N) = aW + aP0;
    
    %for specie B
    MatB(1,1) = aE + aP0 + k02*exp((1-alfa2)*((F)/(R*T))*(E(m)-E0rds));
    MatB(1,2) = -aE;
    for i = 2:N-1
    MatB(i,i) = aE + aW + aP0;
    MatB(i,i-1) = -aW;
    MatB(i,i+1) = -aE;
    end
    MatB(N,N-1) = -aW;
    MatB(N,N) = aW + aP0;
    
    %for specie C
    MatC(1,1) = aE + aP0 + k02*exp(((-F)/(R*T))*(E0rds-E0nonrds))*exp((-(1+alfa2))*(((F)/(R*T))*(E(m)-E0rds)));
    MatC(1,2) = -aE;
    for i = 2:N-1
    MatC(i,i) = aE + aW + aP0;
    MatC(i,i-1) = -aW;
    MatC(i,i+1) = -aE;
    end
    MatC(N,N-1) = -aW;
    MatC(N,N) = aW + aP0;
    
    %% source terms (spatial and temporal sources)
    %for specie A
    SourceA(1) = aP0*CA(1); 
    for i = 2:N-1
       SourceA(i) = aP0*CA(i);
    end
    SourceA(N) = aP0*CA(N);
    
    %for specie B
   SourceB(1) = aP0*CB(1) + CA(1)*k01*exp(-(((1+alfa1)*F)/(R*T))*(E(m)-Ef01)) + CC(1)*k02*exp(((-F)/(R*T))*(E0rds-E0nonrds))*exp((-(1+alfa2))*(((F)/(R*T))*(E(m)-E0rds))); 
    for i = 2:N-1
       SourceB(i) = aP0*CB(i);
    end
    SourceB(N) = aP0*CB(N);
     
    %for specie C
      SourceC(1) = aP0*CC(1) + CB(1)*k02*exp((1-alfa2)*((F)/(R*T))*(E(m)-E0rds)); 
     for i = 2:N-1
        SourceC(i) = aP0*CC(i);
     end
     SourceC(N) = aP0*CC(N);

     % use matrix inversion to find the concentration of species at current time     
     invMA = inv(MatA);
     invMB = inv(MatB);
     invMC = inv(MatC);
     CA = invMA*(SourceA');
     CB = invMB*(SourceB');
     CC = invMC*(SourceC');
     
     % plot concentration profiles in realtime
     figure(1)
     plot(CA,'r');
     hold on
     plot(CB,'b');
     plot(CC,'k');
     hold off
     grid on
     legend('FT (sp A)','{FT_{red}} (sp B)','{FT_{ox}} (sp C)')
     % Calculation of current repsonse at the electrode surface
     I(m) = (F*D*(1/dx)*(4*(-CA(2)+CA(1))+(2*(-CC(2)+CC(1))))*Area); % dimensional current
     pause(0.01)
     
end
%%
toc
Istar = (I*epsilon)/(F*Area*D*Cba); %dimensionless current
%Theta = E*(F/(R*T));

figure(2)
plot(E,Istar)
Istar = Istar';
E = E';
