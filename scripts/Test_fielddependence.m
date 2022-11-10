%% first set the parameters for the recombination model

Prec                        = DV_D4A2_nofield.Prec;                    % initiliase the recombination parameters (default values)
% fixed model parameters for the efficiency limit exploration
Prec.params.tickness        = 100 * 1e-9;
Prec.params.Ex.sigma        = 0.0001;Prec.params.CT.sigma        = 0.0001;
Prec.params.Ex.numbrestate  = 1;Prec.params.CT.numbrestate  = 1;
Prec.const.Edistribution=0.5:0.005:4;%the energy distribution for the absorption spectra%need to be linearly spaced
  Prec.params.Ex.Number_Vibronic_Mode_initial=5;   Prec.params.Ex.Number_Vibronic_Mode_final=15;%number of vibronic mode considered for the ground state
  Prec.params.CT.Number_Vibronic_Mode_initial=5;  Prec. params.CT.Number_Vibronic_Mode_final=15;%number of vibronic mode considered for the ground state
Prec.params.Excitondesnity  = 8e27; %in m-3
%energetics that will be changed along the exploration of the highest
%efficiency
Prec.const.T                = 300;
Prec                        = paramsRec.calcall(Prec); % Update the Recombination Parameters
%% initialise the device parameters


deviceParameterFile = 'DeviceParameters_Default.xlsx';
DP = deviceparams(['parameters\',deviceParameterFile]);
activelayer = 2;
DP.light_properties.OM      = 0; %to consider the transfer matrix generation profile
DP.Time_properties.tpoints  = 100;
DP.Layers{activelayer}.tp   = Prec.params.tickness * 100; % [cm] = [m] * 100
% Active Layer Index                % integer
NC          = 2e19;    %set to the same degeneracy as the CT state
mobility    = DV_D4A2_nofield.DP.Layers{2}.mue;     % Charge Carrier Mobility           % cm^2 / V / s
kdis=DV_D4A2_nofield.DP.Layers{2}.kdis;
kdisex=DV_D4A2_nofield.DP.Layers{2}.kdisexc;
kfor=DV_D4A2_nofield.DP.Layers{2}.kfor;
DP=DP.generateDeviceparams(NC,activelayer,mobility,kdis,kdisex,Prec,kfor,0);


    %% 0D model
    %     kforr=DP.Layers{2}.kfor;
    fignumber=1;
    DP.simulate_PL(Prec,fignumber);
    DP.simulate_EL(Prec,fignumber);
    DP.simulateTAS(2e10,1e27,fignumber);
    figure(fignumber)
    subplot(2,2,1)
    semilogy(Prec.const.Edistribution,Prec.results.AbsLJ)
    xlabel('Energy [eV]')
    ylabel('EQE [a.u]')
    ylim([1*1e-7, 1])
    xlim([1,2])
    sgtitle(['Vocrad= ' num2str(DP.results.Vocrad) ' DVnr= ' num2str(DP.results.DVnr) ' Voc= ' num2str(DP.results.Voc)])
    %             legend("ToT","From CT","FromEx")
    pause(0.1)
    %% Run the JV scans here
    Vstart=0;Vend=1.5;
    counter=3;
    tic
    
    DP.Layers{2}.r0_CT=0; %R0 for field dependence is 1 nm
    DP.Layers{2}.r0_Ex=1e-7; 
    DV2=device(DP);
    DV2.Prec=Prec;
    toc
    %%
    suns=[0.01,0.05,0.1,0.5,1];
    for Gen=suns
        tic
        DV2=device.runsolJsc(DV2,Gen);
        toc
        
        tic
        DV2=device.runsolJV(DV2,Gen,Vstart,Vend);
        toc
        tic
        DV2=device.runsolJV(DV2,Gen,Vstart,-Vend);
        toc
    end
    assignin('base',"DV_D4A2_nofield"+num2str(counter),DV2)
    %%
    figure(2)
    dfplot.JV_new(DV2.sol_JV(9),1);
    dfplot.JV_new(DV2.sol_JV(10),1);
    
    %%
    figure(3)
      counter=0;
    for Gen=suns
        counter=counter+1;
        [tableres.Jsc(counter),tableres.Voc(counter),tableres.FF(counter)]=dfplot.JV_new(DV2.sol_JV(2*(counter-1)+1),1);
    end
    tableres.Jsc=tableres.Jsc';
    tableres.Voc=tableres.Voc';
    tableres.FF=tableres.FF';
    tableres.suns=suns';
    tableres=struct2table(tableres);
    %%
    figure(4)
    loglog(tableres.suns,tableres.suns*max(tableres.Jsc))
    hold on
    loglog(tableres.suns,tableres.Jsc,'*')
    %%
    [Jsc,Voc,FF]=dfplot.JV_new(DV2.sol_JV(2),1);
    figure(3)
    dfplot.photoluminescence_mult(DV2,2,0,1," K");
    figure(4)
    dfplot.Electroluminescence_multi(DV2,2,0,2," K")

