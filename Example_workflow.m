%% EXAMPLE WORKFLOW
% First Calculate the Recombination paramters of the CT and Exciton
% here we also calculate the absorption profile of the device and the
% emission spectra based on the properties of the states 

%% SET PARAMETERS and make Device
startup
% Set and adjust the Recombination Parameters
Prec                        = paramsRec;                    % initiliase the recombination parameters (default values)
offset                      = 0.3;                  % eV    % energy difference between the excited state and the CT state
Prec.params.tickness        = 100 * 1e-9;           % m     % thickness of the active layer
Prec.params.Ex.DG0          = 1.5;                % shift pl and band eqe by the same amount 
Prec.params.CT.DG0          = Prec.params.Ex.DG0 - offset;
Prec.params.Ex.f            = 5;        
Prec.params.CT.f            = 2e-2; %(smaller > down)
Prec.params.Ex.sigma        = 0.01;  % difference between band edge and pl peak position (move PL(ex) right)
Prec.params.CT.sigma        = 0.08;  %move el & CT in PL left/right
Prec.params.Ex.numbrestate  = 1;
Prec.params.CT.numbrestate  = 1;
Prec.params.Ex.L0           = 0.065;  %sharpness pl (does not affect el)   rounder eqe egde
Prec.params.Ex.Li           = 0.03;   %lifts the second bend in PL (does not affect el)
Prec.params.CT.L0           = 0.07;  % sharpness el (smaller = sharper)
Prec.params.CT.Li           = 0.03;  % lifts the secnod bend in el, and pl excited state
Prec.params.RCTE            = 5e-2;
Prec.params.Excitondesnity  = 3e30;  % eqe roudness
Prec.params.Vstar           = 0.001;
Prec.const.T                = 300;
Prec                        = paramsRec.calcall(Prec); % Update the Recombination Parameters

krecCT  = Prec.params.CT.results.knr;
krecex  = Prec.params.Ex.results.knr;
Voc     = Prec.results.Vocrad - Prec.results.Dvnr;
% SIMULATION 0D model

% kforr=DP.Layers{2}.kfor;
fignumber = 1;
DP.simulate_PL(Prec,fignumber);
DP.simulate_EL(Prec,fignumber);
DP.simulateTAS(2e10,1e27,fignumber);

% plot EQE from Prec
figure(fignumber)
subplot(2,2,1)
plot(Prec.const.Edistribution, Prec.results.AbsLJ); hold on
xlabel('Energy [eV]')
ylabel('EQE [a.u]')
ylim([1*1e-7, 1])
xlim([1,2])
sgtitle(['Vocrad= ' num2str(DP.results.Vocrad) ' DVnr= ' num2str(DP.results.DVnr) ' Voc= ' num2str(DP.results.Voc)])

pause(0.1) %give the figure time to finish
%% Generate a device with the defined parameters
% Parameters are from Prec which is defined above and from the PINDevice file, which is loaded below

activelayer = 2;        % Active Layer Index                % integer
NC          = 2e19;     % Number of Charge Carriers         % cm^-3
Kfor        = 1e-10;    % Rate Constant CS to CT            % cm^3 / s
kdis        = 1e12;     % Rate Constant CT dissociation     % cm^3 / s
kdisex      = 1e14;     % Rate Constatn Ex dissociation     % cm^3 / s
mobility    = 5e-4;     % Charge Carrier Mobility           % cm^2 / V / s

deviceParameterFile = uigetfile('parameters\DeviceParameters_Default.xlsx');
if isequal(deviceParameterFile,0)
    deviceParameterFile = 'DeviceParameters_Default.xlsx';
end
DP = deviceparams(['parameters\',deviceParameterFile]);

DP.light_properties.OM      = 0; %to consider the transfer matrix generation profile
DP.Time_properties.tpoints  = 100;
DP.Layers{activelayer}.tp   = Prec.params.tickness * 100; % [cm] = [m] * 100

DP = DP.generateDeviceparams(NC, activelayer, mobility, kdis, kdisex, Prec, Kfor, 0);
clear NC activelayer Tq1exp mobility kdis kdisex


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
    figure(2)
    %dfplot.JV_new(DV2.sol_JV(2),1)
    hold on
    dfplot.JV_new(DV2.sol_JV(end),1)
    
    