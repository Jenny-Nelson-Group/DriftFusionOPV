%% EXAMPLE WORKFLOW
% First Calculate the Recombination paramters of the CT and Exciton
% here we also calculate the absorption profile of the device and the
% emission spectra based on the properties of the states 

%% SET PARAMETERS and make Device
% Set and adjust the Recombination Parameters
%open("C:\Users\jolan\DATA_Local\01 Documents\03 Work\04 Imperial PhD\10 TO SORT\Exp_CN-Tp_3.fig")
handle = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', 1);
if isempty(handle)
    %open("C:\Users\jsm121\LOCAL DATA\Imperial PhD\10 TO SORT\Exp_all_3_nice.fig")
    %open("C:\Users\jolan\DATA_Local\01 Documents\03 Work\04 Imperial PhD\10 TO SORT\Exp_all_3_nice.fig")
    open("C:\Users\"+getenv('username')+"\OneDrive - Imperial College London\Imperial PhD\10 TO SORT\Exp_all_3_nice.fig")
end
addpath(genpath(pwd)); % add folders to path

%%
Prec                        = paramsRec;                    % initiliase the recombination parameters (default values)
offset                      = 0.14;                  % eV    % energy difference between the excited state and the CT state
Prec.params.tickness        = 100 * 1e-9;           % m     % thickness of the active layer
Prec.params.Ex.DG0          = 1.35;                 
Prec.params.CT.DG0          = Prec.params.Ex.DG0 - offset;
Prec.params.Ex.f            = 2.56;
Prec.params.CT.f            = 1e-4;
Prec.params.Ex.sigma        = 0.01;
Prec.params.CT.sigma        = 0.01;
Prec.params.Ex.numbrestate  = 1;
Prec.params.CT.numbrestate  = 1;
Prec.params.Ex.L0           = 0.075;  %0.10
Prec.params.Ex.Li           = 0.05; %0.15   %0.04-0.150
Prec.params.CT.L0           = 0.07;  %0.18  %CT smoothing 
Prec.params.CT.Li           = 0.1;   %0.15
Prec.params.RCTE            = 0.1;  
Prec.params.Excitondesnity  = 3e27;
Prec.params.Vstar           = 0.0; % hybridisation: on: 0.01 /0.02 coupling between CT and ex
Prec.const.T                = 300;
Prec                        = paramsRec.calcall(Prec); % Update the Recombination Parameters

krecCT  = Prec.params.CT.results.knr  +  Prec.params.CT.results.krTot;
krecex  = Prec.params.Ex.results.knr  +  Prec.params.Ex.results.krTot;
Voc     = Prec.results.Vocrad - Prec.results.Dvnr;

% Generate a device with the defined parameters
%% Parameters are from Prec which is defined above and from the PINDevice file, which is loaded below

activelayer = 2;        % Active Layer Index                % integer
NC          = 1e18;     % Number of Charge Carriers         % cm^-3
Bfor        = 5e-11;    % Rate Constant CS to CT            % cm^3 / s
kdisEx      = 6e8;      % Rate Constatn Ex dissociation     % 1 / s
kdisCT      = 3e7;      % Rate Constant CT dissociation     % 1 / s
mobility    = 5e-5;     % Charge Carrier Mobility           % cm^2 / V / s
Rshunt      = 10^5;     % ohm cm2 (?)

delta_kEx = kdisEx - krecex;
delta_kCT = kdisCT - krecCT;
% deviceParameterFile = uigetfile('parameters\DeviceParameters_Default.xlsx');
% if isequal(deviceParameterFile,0)
%     deviceParameterFile = 'DeviceParameters_Default.xlsx';
% end
deviceParameterFile = 'DeviceParameters_Qatar.xlsx';
DP = deviceparams(['parameters\',deviceParameterFile]);

DP.light_properties.OM      = 0; %to consider the transfer matrix generation profile
DP.Time_properties.tpoints  = 100;
DP.Layers{activelayer}.tp   = Prec.params.tickness * 100; % [cm] = [m] * 100
DP.Layers{2}.r0_CT= 0;
DP.Layers{2}.r0_Ex = 0;

DP.External_prop.Rseries=0;
DP.Experiment_prop.BC =3;
sn=1e8;
DP.External_prop.sn_l=sn;%left Electron extraction/surface recombination coefficient in cm s^-1
DP.External_prop.sn_r=sn;%right Electron extraction/surface recombination coefficient in cm s^-1
DP.External_prop.sp_l=sn;%left hole extraction/surface recombination coefficient in cm s^-1
DP.External_prop.sp_r=sn;%right hole extraction/surface recombination coefficient in cm s^-1

DP = DP.generateDeviceparams(NC, activelayer, mobility, kdisCT, kdisEx, Prec, Bfor, 0);
clear activelayer Tq1exp mobility kdis kdisex

% SIMULATION 0D model
% AAA
% kforr=DP.Layers{2}.kfor;
fignumber = 1;
%DP.simulate_PL(Prec,fignumber);
DP.simulate_EL(Prec,fignumber);
%DP.simulateTAS(2e10,1e27,fignumber);

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


% plot EQE from Prec
figure(fignumber)
subplot(2,2,1)
semilogy(Prec.const.Edistribution, Prec.results.AbsLJ,'LineWidth',2,'Color',[1,0,0]); hold on
xlabel('Energy [eV]')
ylabel('EQE [a.u]')
ylim([1*1e-7, 1.5])
xlim([1,2])
sgtitle(['Vocrad= ' num2str(DP.results.Vocrad) ' DVnr= ' num2str(DP.results.DVnr) ' Voc= ' num2str(DP.results.Voc)])

pause(0.1) %give the figure time to finish

%% Run the JV scans here
    Vstart  = 0;
    Vend    = 1.5;
    tic
    DP.Layers{2}.r0=0; %R0 for field dependence is 1 nm

    DV2=device(DP);
    DV2.Prec=Prec;
    toc

    % for different suns
    suns = [1];
    for Gen=suns
        tic
        DV2=device.runsolJsc(DV2,Gen);
        %DV2=device(DP);
        toc

        tic
        DV2=device.runsolJV(DV2,Gen,Vstart,Vend);
        toc
    end
    
    %%
    % plot JV
    figure(fignumber)
    subplot(2,2,4)
    %dfplot.JV_new(DV2.sol_JV(2),1)
    hold on
    [Jsc,Voc,FF,JJ,VV] = dfplot.JV_new(DV2.sol_JV(end),1,Rshunt);
    %title("Jsc ="+Jsc+", FF = "+FF+", Voc = "+Voc)
    hold on
    %plot([0,1],[0,0]);
    xlim([-0., 1]);
%%
    figure(fignumber)
    subplot(2,2,2)
    hold on
    dfplot.photoluminescence(DV2,2,0,1);
    xlim([800,1400]);
    ylim([1e-1,1e3]);

%%
DVnr_new = DV2.sol_JV.params.results.Vocrad - Voc;