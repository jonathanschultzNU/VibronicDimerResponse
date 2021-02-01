%Calculation of first- and third-order response of a molecular dimer using a Frenkel-Exciton,
%time-independent Hamiltonian

%	Author: Jonathan D. Schultz
%	Email: jonathanschultz2022@u.northwestern.edu
%	Last revision date: February 1st, 2021
%
%	Copyright: Jonathan D. Schultz, 2021

%   Please see readme file for information about this package

%   Notes: 
%   This code uses several structures to to minimize workspace clutter. The
%   structure are as follows:
%         Par = parameters
%         FEH = Frenkel-exciton Hamiltonian; anything related to the Hamiltonian generation
%         Res = Response; anything related to calculating the time-dependent system response function

close all;
clear;
clear var;

%Define constants and parameters

Par.c = 2.9979 * 10^-5;      % speed of light [cm/fs]
Par.e1 = 17000;              % center energy of local excited electronic state for molecule 1 [cm-1]
Par.e2 = 17000;              % center energy of local excited electronic state for molecule 2 [cm-1]
Par.w0 = 1300;               % vibration coupled to electronic state [cm-1]
Par.J = 200;                 % dipolar coupling [cm-1]
Par.lambda = sqrt(0.6);      % linear vibronic coupling
Par.v_max = 3;               % maximum vibrational quanta included in Hamitonian

%Calculate Hamiltonian
FEH = DimHamGen(Par);        %external function

%% Prepare for first-order response calculation

Res.tmax = 128;         % max time for coherence and rephasing time vectors (fs)
Res.dt = 2;             % step size (fs)

Res.nsteps = Res.tmax/Res.dt+1;   
Res.Tdeph = 20;         % dephasing time constant for coherence and rephasing time dimensions

Res.t1 = (0:Res.dt:Res.tmax);        %coherence time vector

Res.R_t1 = zeros(1,length(Res.t1));   %pre-allocate R1 response function matrix

Res.Use = @(t) expm(1i*FEH.Hse*t*2*pi);  %singly-excited-state propagator for numerical integration

%assuming no thermal distribution
Res.ket = FEH.fockgg(:,1);    
Res.bra = Res.ket';

%% Calculate the linear response

for i = 1:length(Res.t1)
    Res.R_t1(i) = Res.bra*(FEH.MUse_g')*Res.Use(Res.t1(i))*(FEH.MUg_se')*Res.ket.*exp(-Res.t1(i)/Res.Tdeph);
end

%% Convert linear response to the frequency domain and plot

Res.n1=2^8;
Res.dw1 = 1 / (Par.c*Res.dt*Res.n1);                       %defining the frequency step [cm-1]
Res.w1 = ((-Res.n1*Res.dw1)/2:Res.dw1:((Res.n1*Res.dw1)/2-Res.dw1));   %frequency axis [cm-1%]
Res.w1=Res.w1+Par.e1;                                  %shift energies back by how much was removed from diagonal

Res.R1_w1 = real(fftshift(fft(Res.R_t1,Res.n1)));

figure
plot(Res.w1,normdim(Res.R1_w1),'k','LineWidth',2)
xlim([15000 20000])
ylim([0 1.1])
xlabel('\omega/2\pic (cm^{-1})');
ylabel('linear response (a.u.)');
set(gca,'FontSize',20,'YTickLabel',[]);

%% Prepare for third-order response calculation

Res.tmax = 128;         % max time for coherence and rephasing time vectors (fs)
Res.dt = 2;             % step size (fs)
%PLEASE NOTE: for 65 coherence and rephasing timepoints and 3 vibrational quanta, the model takes ~16
%seconds per spectrum to compute

Res.dt2 = Res.dt;           % waiting time step size (fs)
Res.nsteps = Res.tmax/Res.dt+1;   
Res.t2max = 0;          % maximum waiting time value
Res.Tdeph = 20;         % dephasing time constant for coherence and rephasing time dimensions
Res.T2deph = 80;        % dephasing time for waiting time

Res.t1 = (0:Res.dt:Res.tmax);        %coherence time vector
Res.t2 = (0:Res.dt2:Res.t2max);      %waiting time vector
Res.t3 = (0:Res.dt:Res.tmax);        %rephasing time vector

Res.R1_t1t3 = zeros(Res.nsteps,Res.nsteps,length(Res.t2));   %pre-allocate R1 response function matrix
Res.R2_t1t3 = zeros(Res.nsteps,Res.nsteps,length(Res.t2));   %pre-allocate R2 response function matrix
Res.R3_t1t3 = zeros(Res.nsteps,Res.nsteps,length(Res.t2));   %pre-allocate R3 response function matrix
Res.R4_t1t3 = zeros(Res.nsteps,Res.nsteps,length(Res.t2));   %pre-allocate R4 response function matrix

Res.Ugg = @(t) expm(-1i*FEH.Hgg*t*2*pi);  %ground-state propagator for numerical integration 
Res.Use = @(t) expm(-1i*FEH.Hse*t*2*pi);  %singly-excited-state propagator for numerical integration

%assuming no thermal distribution
Res.ket = FEH.fockgg(:,1);    
Res.bra = Res.ket';

%% Calculate third-order response

%time bottleneck for the code - loop through every index of all response pathways

for k = 1:length(Res.t2)
    for j = 1:length(Res.t1)    %column index is coherence time
        for i= 1:length(Res.t3) %row index is rephasing time

            %R1: Non-rephasing ground state bleach
            Res.R1_t1t3(i,j,k) = Res.bra*FEH.MUse_g'*Res.Use(Res.t3(i))*FEH.MUg_se'*Res.Ugg(Res.t2(k))*FEH.MUse_g'*Res.Use(Res.t1(j))*FEH.MUg_se'*Res.ket*exp(-(Res.t3(i)+Res.t1(j))/(Res.Tdeph))*exp(-Res.t2(k)/Res.T2deph); 

            %R2: Rephasing stimulated emission
            Res.R2_t1t3(i,j,k) = Res.bra*FEH.MUg_se*conj(Res.Use(Res.t1(j)))*conj(Res.Use(Res.t2(k)))*FEH.MUse_g*(FEH.MUse_g')*Res.Use(Res.t3(i))*Res.Use(Res.t2(k))*FEH.MUg_se'*Res.ket*exp(-(Res.t3(i)+Res.t1(j))/(Res.Tdeph))*exp(-Res.t2(k)/Res.T2deph); 

            %R3: Rephasing ground state bleach
            Res.R3_t1t3(i,j,k) = Res.bra*FEH.MUg_se*conj(Res.Use(Res.t1(j)))*FEH.MUse_g*conj(Res.Ugg(Res.t2(k)))*FEH.MUse_g'*Res.Use(Res.t3(i))*FEH.MUg_se'*Res.ket*exp(-(Res.t3(i)+Res.t1(j))/(Res.Tdeph))*exp(-Res.t2(k)/Res.T2deph);
 
            %R4: Non-rephasing stimulated emission
            Res.R4_t1t3(i,j,k) = Res.bra*FEH.MUg_se*conj(Res.Use(Res.t2(k)))*FEH.MUse_g*conj(Res.Ugg(Res.t3(i)))*FEH.MUse_g'*Res.Use(Res.t3(i))*Res.Use(Res.t2(k))*Res.Use(Res.t1(j))*FEH.MUse_g*Res.ket*exp(-(Res.t3(i)+Res.t1(j))/(Res.Tdeph))*exp(-Res.t2(k)/Res.T2deph);            
            
        end
    end
end


%% Convert third-order response to the frequency domain

% Trap rule correction (see Hamm and Zanni book, chapter 9)

for k = 1:length(Res.t2)
    Res.R1_t1t3(:,1,k) = Res.R1_t1t3(:,1,k)./2;   Res.R1_t1t3(1,:,k) = Res.R1_t1t3(1,:,k)./2;
    Res.R2_t1t3(:,1,k) = Res.R2_t1t3(:,1,k)./2;   Res.R2_t1t3(1,:,k) = Res.R2_t1t3(1,:,k)./2;
    Res.R3_t1t3(:,1,k) = Res.R3_t1t3(:,1,k)./2;   Res.R3_t1t3(1,:,k) = Res.R3_t1t3(1,:,k)./2;
    Res.R4_t1t3(:,1,k) = Res.R4_t1t3(:,1,k)./2;   Res.R4_t1t3(1,:,k) = Res.R4_t1t3(1,:,k)./2;
end

Res.n3 = 2^8; Res.n1 = 2^8; 

%preallocate frequency-domain matrices
Res.R1_w1w3 = zeros(Res.n1,Res.n3,length(Res.t2));
Res.R2_w1w3 = zeros(size(Res.R1_w1w3));
Res.R3_w1w3 = zeros(size(Res.R1_w1w3));
Res.R4_w1w3 = zeros(size(Res.R1_w1w3));
Res.R_GSB_abs = zeros(size(Res.R1_w1w3));       %ground state bleach
Res.R_SE_abs = zeros(size(Res.R1_w1w3));        %stimulated emission

temp.R1 = zeros(Res.n1,Res.n3);
temp.R2 = zeros(size(temp.R1));
temp.R3 = zeros(size(temp.R1));
temp.R4 = zeros(size(temp.R1));

for k = 1:length(Res.t2)
    temp.R1 = fft2(Res.R1_t1t3(:,:,k),Res.n3,Res.n1);   
    temp.R1 = rot90(temp.R1,2);             %not particularly sure why I have to do this step, but it works
    Res.R1_w1w3(:,:,k) = -(fftshift(temp.R1));

    temp.R2 = fft2(Res.R2_t1t3(:,:,k),Res.n3,Res.n1); 
    temp.R2 = fliplr(circshift(temp.R2,-1,2)); 
    temp.R2 = rot90(temp.R2,2);
    Res.R2_w1w3(:,:,k) = -fftshift(temp.R2);

    temp.R3 = fft2(Res.R3_t1t3(:,:,k),Res.n3,Res.n1); 
    temp.R3 = fliplr(circshift(temp.R3,-1,2));
    temp.R3 = rot90(temp.R3,2);
    Res.R3_w1w3(:,:,k) = -fftshift(temp.R3);

    temp.R4 = fft2(Res.R4_t1t3(:,:,k),Res.n3,Res.n1);  
    temp.R4 = rot90(temp.R4,2);
    Res.R4_w1w3(:,:,k) = -(fftshift(temp.R4));
    
    Res.R_GSB_abs(:,:,k) = real(Res.R1_w1w3(:,:,k))+real(Res.R3_w1w3(:,:,k));
    Res.R_SE_abs(:,:,k) = real(Res.R2_w1w3(:,:,k))+real(Res.R4_w1w3(:,:,k));
    
end

% Define frequency axes

Res.dw3 = 1 / (Par.c*Res.dt*Res.n3);                       %defining the frequency step [cm-1]
Res.w3 = ((-Res.n3*Res.dw3)/2:Res.dw3:((Res.n3*Res.dw3)/2-Res.dw3));   %frequency axis [cm-1%]
Res.w3=Res.w3+Par.e1;                                  %shift energies back by how much was removed from diagonal

Res.dw1 = 1 / (Par.c*Res.dt*Res.n1);                       %defining the frequency step [cm-1]
Res.w1 = ((-Res.n1*Res.dw1)/2:Res.dw1:((Res.n1*Res.dw1)/2-Res.dw1));   %frequency axis [cm-1%]
Res.w1=Res.w1+Par.e1;                                  %shift energies back by how much was removed from diagonal
 

%% Plot frequency-domain third-order response

Res.indw1 = Res.w1>15000 & Res.w1<20000;  % wide view
Res.indw3 = Res.w3>15000 & Res.w3<20000;

Res.tsel = [1];     %select the desired timepoints to plot

fig = figure;
for i=1:length(Res.tsel)
subplot1 = subplot(1,length(Res.tsel),i,'Parent',fig);
hold(subplot1,'on');
MDplot(Res.w1(Res.indw1),Res.w1(Res.indw3),normdim(real(Res.R1_w1w3(Res.indw3,Res.indw1,Res.tsel(i)))),Res.t2(i));
end