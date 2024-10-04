%% EXAMPLE WORKFLOW
% First Calculate the Recombination paramters of the CT and Exciton
% here we also calculate the absorption profile of the device and the
% emission spectra based on the properties of the states 

%% Setup
% add folders to path
addpath(genpath(pwd));  

% open figure
% handle = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', 1);
% if isempty(handle)
%     open(pwd + "\Example_Data.fig")
% end

%%
for i = 1:length(figure(fignumber).Children)
    set(figure(fignumber).Children(i).Children, {'color'}, num2cell(bone2,2));
end

%%
for i = 1:length(gcf().Children)
    set(gcf().Children(i).Children, {'color'}, num2cell(flipud(jet(11)),2));
end

%% SET PARAMETERS and make Device
Prec                        = paramsRec;                    % initiliase the recombination parameters (default values)
offset                      = 0.01;                 % eV    % energy difference between the excited state and the CT state
Prec.params.tickness        = 100 * 1e-9;           % m     % thickness of the active layer
Prec.params.Ex.DG0          = 1.35;                 
Prec.params.CT.DG0          = Prec.params.Ex.DG0 - offset;
Prec.params.Ex.f            = 2;
Prec.params.CT.f            = 1e-10;
Prec.params.Ex.sigma        = 0.1;
Prec.params.CT.sigma        = 0.01;
Prec.params.Ex.numbrestate  = 5;
Prec.params.CT.numbrestate  = 5;
Prec.params.Ex.L0           = 0.10;  %0.10
Prec.params.Ex.Li           = 0.05;   %0.15   %0.04-0.150
Prec.params.CT.L0           = 0.07;   %0.18   %CT smoothing 
Prec.params.CT.Li           = 0.1;    %0.15
Prec.params.RCTE            = 1;  
Prec.params.Excitondesnity  = 3e27;
Prec.params.Vstar           = 0.0; % hybridisation: on: 0.01 /0.02 coupling between CT and ex
Prec.const.T                = 100;
Prec                        = paramsRec.calcall(Prec); % Update the Recombination Parameters

krecCT    = Prec.params.CT.results.knr  +  Prec.params.CT.results.krTot;
krecex    = Prec.params.Ex.results.knr  +  Prec.params.Ex.results.krTot;
Voc_Prec  = Prec.results.Vocrad - Prec.results.Dvnr;

figure(5)
DP.simulate_PL(Prec,5);

% Generate a device with the defined parameters
%% Parameters are from Prec which is defined above and from the PINDevice file, which is loaded below

activelayer = 2;        % Active Layer Index                % integer
NC          = 1e18;     % Number of Charge Carriers         % cm^-3
Bfor        = 5e-11;    % Rate Constant CS to CT            % cm^3 / s
kdisEx      = 1e10;      % Rate Constatn Ex dissociation     % 1 / s
kdisCT      = 5e7;    % Rate Constant CT dissociation     % 1 / s
mobility    = 9e-5;     % Charge Carrier Mobility           % cm^2 / V / s
Rshunt      = 10^5;     % ohm cm2 (?)

delta_kEx = kdisEx - krecex;
delta_kCT = kdisCT - krecCT;
deviceParameterFile = 'DeviceParameters_Default.xlsx';
DP = deviceparams(['parameters\',deviceParameterFile]);

DP.light_properties.OM      = 0;
DP.Time_properties.tpoints  = 100;
DP.Layers{activelayer}.tp   = Prec.params.tickness * 100; % [cm] = [m] * 100
DP.Layers{2}.r0_CT          = 0;
DP.Layers{2}.r0_Ex          = 0;

DP.External_prop.Rseries    = 0;
DP.Experiment_prop.BC       = 3;
sn                          = 1e8;
DP.External_prop.sn_l       = sn;  %left Electron extraction/surface recombination coefficient in cm s^-1
DP.External_prop.sn_r       = sn;  %right Electron extraction/surface recombination coefficient in cm s^-1
DP.External_prop.sp_l       = sn;  %left hole extraction/surface recombination coefficient in cm s^-1
DP.External_prop.sp_r       = sn;  %right hole extraction/surface recombination coefficient in cm s^-1 
DP.physical_const.T         = Prec.const.T;

DP = DP.generateDeviceparams(NC, activelayer, mobility, kdisCT, kdisEx, Prec, Bfor, 0);
clear activelayer Tq1exp mobility kdis kdisex


%% Run the JV scans here
Vstart  = 0;
Vend    = 1.5;
tic
DP.Layers{2}.r0 = 0; % R0 for field dependence is 1 nm
DV2 = device(DP);
DV2.Prec = Prec;
toc

% for different suns
suns = 1;
for Gen = suns
    tic
    DV2 = device.runsolJsc(DV2,Gen);
    %DV2 = device(DP);
    toc
    
    tic
    DV2 = device.runsolJV(DV2,Gen,Vstart,Vend);
    toc
end
    
%% PLOTTING
fignumber = 1;

% Plot JV
figure(fignumber)
subplot(2,3,1)
[Jsc,Voc,FF,JJ,VV] = dfplot.JV_new(DV2.sol_JV(end),1,Rshunt);
xlim([-0, 1]);
box on

% plot EQE from Prec
figure(fignumber)
subplot(2,3,4)
semilogy(Prec.const.Edistribution, Prec.results.AbsLJ,'LineWidth',2,'Color',[1,0,0]); hold on
xlabel('Energy (eV)')
ylabel('EQE')
ylim([1*1e-7, 1.5])
xlim([0.5,2])
%sgtitle(['Vocrad= ' num2str(DP.results.Vocrad) ' DVnr= ' num2str(DP.results.DVnr) ' Voc= ' num2str(DP.results.Voc)])
box on

% Plot Photoluminescence (from full simulation)
figure(fignumber)
subplot(2,3,2)
hold on
dfplot.photoluminescence(DV2,2,0,1);
xlim([0.5,2]);
ylim([1e-1,1e3]);
set(gca,'YScale','log')
legend off
box on

%DP.simulate_PL(Prec,fignumber);

% Plot EL
figure(fignumber)
subplot(2,3,5)
DP.simulate_EL(Prec,fignumber);

% PLOT TAS
subplot(2,3,3)
DP.simulateTAS(2e10,1e27,fignumber);
xlabel("Time")
ylabel("TAS")
legend off


%% Output Values
kbT = DP.physical_const.kB*300;
CT0 = DP.results.CT0;
ni  = sqrt(CT0*kdisCT/Bfor);

E1_Ex = Prec.params.Ex.DG0;
E2_CT = Prec.params.CT.DG0;
E3_CS = 2*kbT*log(NC/ni);
%ECS_IP = DP.Layers{2}.EA - DP.Layers{2}.IP;

D1_ex_ct = offset;
D2_ct_cs = E2_CT - E3_CS;
D3_ex_cs = D1_ex_ct + D2_ct_cs;

disp("--------")
DVnr_new = DV2.sol_JV.params.results.Vocrad - Voc;
disp("Voc (full simulation) = " + Voc + " V");
disp("Non radiative voltage loss (full simulation) = " + DVnr_new + " V");
disp("--------")

