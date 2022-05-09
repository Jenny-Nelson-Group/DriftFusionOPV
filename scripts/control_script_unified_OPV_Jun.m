%% this is mohammed's code for running an unified model in OPV: EES2022
%% this model includes a molecular model, based on paper PRX2018, JACS2019, NC2021, a Kinetic model, and a drift diffusion model
%% using this code, we can run molecular model to calculate the rates of transition first, then kinetic model,
%% then drift diffusion model to get populations of CT and LE, as well as charge transport
%% modified by Jun Yan

clc
clear all
close all

%% First Calculate the Recombination paramters of the CT and Exciton
%here we also calculate the absorption profile of the device and the
%emission spectra based on the properties of the states, 
Prec=paramsRec;%first initiliase the recombination parameters ( the parameters are set here by default and can be changed under or in the paramRec function
%below are a set of parameters that can be changed

%% input for LE and CT properties in the molecular model
Prec.params.tickness=100e-9;% active layer thickness in m
Prec.const.T=300;% temperature of the cell
Prec.params.RCTE=1e0;%ratio CT to S1
Prec.params.Excitondesnity=1e28;% exciton density 
%% EXciton properties          
Prec.params.Ex.hW=0.15;%main vibronic energy considered in eV
%Prec.params.Ex.hW2=0.012;%low frequency vibronic energy considered in eV
Prec.params.Ex.f=5;%oscillator strength of the CT state
Prec.params.Ex.L0=0.1;%outer reorganisation energy(High frequency) in eV
Prec.params.Ex.Li=0.1;%inner reorganisation energy in eV
Prec.params.Ex.DG0=1.5;%energy of the CT state in  eV
Prec.params.Ex.sigma=0.001;%gaussian distribution for CT state
Prec.params.Ex.numbrestate=1;
Prec.params.Ex.Dmus=3*3.33e-30/1.6e-19;%difference in static dipole moment (10 in DEbye ) then th
%% CT properties
Prec.params.CT.hW=0.15;%main vibronic energy considered in eV
%Prec.params.CT.hW2=0.012;%low frequency vibronic energy considered in eV
Prec.params.CT.f=1e-2;%oscillator strength of the CT state
Prec.params.CT.L0=0.1;%outer reorganisation energy(High frequency) in eV
Prec.params.CT.Li=0.1;%inner reorganisation energy in eV
Prec.params.CT.DG0=1.4;%energy of the CT state in  eV
Prec.params.CT.sigma=0.001;%gaussian distribution for CT state
Prec.params.CT.numbrestate=1;
Prec.params.CT.Dmus=3*3.33e-30/1.6e-19;%difference in static dipole moment (10 in DEbye ) then th
%%%%%%% add this to account for the effect of Hybredisation
Prec.params.Vstar=0.000; % Coupling between S1 and CT in eV

%% this line calculated the different properties of the states.
Prec=paramsRec.calcall(Prec);
krecCT=Prec.params.CT.results.knr+Prec.params.CT.results.krTot;
krecex=Prec.params.Ex.results.knr+Prec.params.Ex.results.krTot;
%Voc=Prec.results.Vocrad-Prec.results.Dvnr;

%% in this part we generate a device with a certain number of properties
%the default properties are set in pnParamsHCT or in the excel file p3hTPCBM.xlsx
runjv=0;
NC=1e20;activelayer=2;Kfor=1e-11;%V in V, K in S-1, NC in Cm-3, Jsc in mA cm-2,
mobility=1e-3;kdis=1e11;kdisex=1e11;% Tq1 in s,mobility in Cm2V-1s-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deviceParameterFile='DeviceParameters_Default.xlsx';
DP = deviceparams(['parameters\',deviceParameterFile]);

DP.light_properties.OM=0;%to consider the transfer matrix generation profile
DP.Time_properties.tpoints=100;
DP.Layers{activelayer}.tp=Prec.params.tickness*100;% in cm
DP=DP.generateDeviceparams(NC,activelayer,mobility,kdis,kdisex,Prec,1e-11,0);






if runjv==1
%% JV simulations
clear NC activelayer Tq1exp mobility kdis kdisex
%% Run the JV scans here
Vstart=0;Vend=2;
tic
DP.Layers{2}.r0=0; %R0 for field dependence is 1 nm

DV2=device(DP);

DV2=device.runsolJsc(DV2,1);
    toc
% light JV
    tic
    DV2=device.runsolJV(DV2,1,Vstart,Vend);
    toc

    assignin('base',"DV_light",DV2)
    figure(2)
    
dfplot.JV_new(DV2.sol_JV(1),1,0)
hold on
lg=legend();
lg.String{end}='light';

%%
clear EL_maxint PL_maxint
%close all
kk=0;
variablelist=[1,5,10,50,100];
for var=variablelist
%     DP.External_prop.Rseries=0;
%     DP.Experiment_prop.BC =3;
%     sn=1e5;
%     DP.External_prop.sn_l=sn;%left Electron extraction/surface recombination coefficient in cm s^-1
%     DP.External_prop.sn_r=sn;%right Electron extraction/surface recombination coefficient in cm s^-1
%     DP.External_prop.sp_l=sn;%left hole extraction/surface recombination coefficient in cm s^-1
%     DP.External_prop.sp_r=sn;%right hole extraction/surface recombination coefficient in cm s^-1
    tic
NC=1e20;activelayer=2;Kfor=1e-11;%V in V, K in S-1, NC in Cm-3, Jsc in mA cm-2,
mobility=1e-3;kdis=1e11;kdisex=1e11;% Tq1 in s,mobility in Cm2V-1s-1
    DP=DP.generateDeviceparams(NC,activelayer,mobility/var,kdis,kdisex,Prec,1e-11,0);

    DV2=device(DP,DV_light.sol_eq);
    DV2.Prec=Prec;
    toc
    
    G=1e-6;
    tic
    DV2=device.runsolJsc(DV2,G);
    toc
    
    tic
    DV2=device.runsolJV(DV2,G,Vstart,Vend);
    toc

    assignin('base',"DV_MT"+num2str(var),DV2)
    
    %%
    kk=kk+1;
    figure(11)
    subplot(2,2,1)
dfplot.JV_new(DV2.sol_JV(1),1,0)
ax=gca;
ax.YScale='log';
ax.YLim=[0,1e-1];
hold on
lg=legend();
lg.String{end}=num2str(var);%num2str(log(Rseries)/log(10));
subplot(2,2,2)
[~,ELplot]=dfplot.Electroluminescence_multi(DV2,2,0,0.1,num2str(var),0);
EL_maxint(kk)=max(ELplot);
hold on 
xlim([0.9,1.6])

subplot(2,2,3)
[~,PLplot,~]=dfplot.photoluminescence_mult(DV2,2,1,1e-6,num2str(var));
xlim([0.9,1.6])
PL_maxint(kk)=max(PLplot);
hold on 

subplot(2,2,4)
plot(variablelist(1:kk),PL_maxint/max(PL_maxint))
hold on
plot(variablelist(1:kk),EL_maxint/max(EL_maxint))
hold off
legend('maxPL','maxEL')
end

%% Plot the results for the series resitance

kk=0;
clear EL_maxint PL_maxint
variablelist=[1,100,500,1000];

for var=variablelist
    kk=kk+1;
%    eval("DV2=DV_R"+num2str(var));
   figure(22)
    subplot(2,2,1)
dfplot.JV_new(DV2.sol_JV(1),1,var)
ax=gca;
ax.YScale='log';
ax.YLim=[1e-7,1e-1];
ax.XLim=[0,2];
hold on
lg=legend();
lg.String{end}=num2str(var);%num2str(log(var)/log(10));
lg.Title.String='series resitance ohm cm-2';
title('dark JV')
subplot(2,2,2)
[~,ELplot]=dfplot.Electroluminescence_multi(DV2,2,0,0.1,num2str(var),var);
EL_maxint(kk)=max(ELplot);
hold on 
xlim([0.9,2])
lg=legend;
lg.Title.String='series resitance ';
subplot(2,2,3)
[~,PLplot,~]=dfplot.photoluminescence_mult(DV2,2,0.5,1e-6,num2str(var));
xlim([0.9,2])
PL_maxint(kk)=max(PLplot);
hold on 
lg=legend;
lg.Title.String='series resitance';
subplot(2,2,4)
plot(variablelist(1:kk),PL_maxint/max(PL_maxint),'*-')
hold on
plot(variablelist(1:kk),EL_maxint/max(EL_maxint),'o-')
hold off
ylim([0,1.1])
legend('maxPL','maxEL')
 xlabel('series resitance (ohm cm-2)')
 ax=gca;
 ax.XScale='log';
end
%% Plot the results for the shunt resitance
kk=0;
clear EL_maxint PL_maxint
variablelist=[100,500,1000,5000,10000];
%eval("DV2=DV_1;");
for var=variablelist
    kk=kk+1;
%     eval("DV2=DV_R"+num2str(Rseries));
   figure(33)
    subplot(2,2,1)
dfplot.JV_new(DV2.sol_JV(1),1,0,var)
ax=gca;
ax.YScale='log';
ax.YLim=[1e-7,1e-1];
ax.XLim=[0,2];
hold on
lg=legend();
lg.String{end}=num2str(var);%num2str(log(Rseries)/log(10));
lg.Title.String='shunt resitance ohm cm-2';
title('dark JV')
subplot(2,2,2)
[~,ELplot]=dfplot.Electroluminescence_multi(DV2,2,0,2,num2str(var),0,var);
EL_maxint(kk)=max(ELplot);
hold on 
xlim([0.9,2])
lg=legend;
lg.Title.String='shunt resitance ';
subplot(2,2,3)
[~,PLplot,~]=dfplot.photoluminescence_mult(DV2,2,0.5,1e-6,num2str(var));
xlim([0.9,2])
PL_maxint(kk)=max(PLplot);
hold on 
lg=legend;
lg.Title.String='shunt resitance';
subplot(2,2,4)
plot(variablelist(1:kk),PL_maxint/max(PL_maxint),'*-')
hold on
plot(variablelist(1:kk),EL_maxint/max(EL_maxint),'o-')
hold off
ylim([0,1.1])
legend('maxPL','maxEL')
 xlabel('shunt resitance (ohm cm-2)')
 ax=gca;
 ax.XScale='log';
end

%% Plot the results for the mobility
kk=0;
clear EL_maxint PL_maxint EL_maxint_05 EL_maxint_1 EL_maxint_2
variablelist=[1,5,10,50,100];
%eval("DV2=DV_1;");
 
for var=variablelist
    kk=kk+1;
    eval("DV2=DV_MT"+num2str(var));
  figure(44)
    subplot(2,2,1)
dfplot.JV_new(DV2.sol_JV(1),1,0)
ax=gca;
ax.YScale='log';
ax.YLim=[1e-7,1e-1];
ax.XLim=[0,2];
hold on
lg=legend();
lg.String{end}=num2str(mobility/var,'%1.0e');%num2str(log(Rseries)/log(10));
lg.Title.String='mobility  ohm cm2 V-1 s-1';
title('dark JV')
subplot(2,2,2)
[~,ELplot]=dfplot.Electroluminescence_multi(DV2,2,0,2,num2str(mobility/var,'%1.0e'),0);
EL_maxint_2(kk)=max(ELplot);
hold on 

xlim([0.9,2])
lg=legend;
lg.Title.String='mobility ';
subplot(2,2,3)
[~,PLplot,~]=dfplot.photoluminescence_mult(DV2,2,0.5,1e-6,num2str(mobility/var,'%1.0e'));
xlim([0.9,2])
PL_maxint(kk)=max(PLplot);
hold on 
lg=legend;
lg.Title.String='mobility';
subplot(2,2,4)
plot(mobility./variablelist(1:kk),PL_maxint/max(PL_maxint),'*-')
hold on
plot(mobility./variablelist(1:kk),EL_maxint_2/max(EL_maxint_2),'o-')
hold off
ylim([0,1.1])
figure(3)
[~,ELplot]=dfplot.Electroluminescence_multi(DV2,2,0,1,num2str(mobility/var,'%1.0e'),0);
EL_maxint_1(kk)=max(ELplot);
[~,ELplot]=dfplot.Electroluminescence_multi(DV2,2,0,0.5,num2str(mobility/var,'%1.0e'),0);
EL_maxint_05(kk)=max(ELplot);
% [~,ELplot]=dfplot.Electroluminescence_multi(DV2,2,0,5,num2str(mobility/Rseries,'%1.0e'));
% EL_maxint_5(kk)=max(ELplot);
figure(44)
subplot(2,2,4)
plot(mobility./variablelist(1:kk),PL_maxint/max(PL_maxint),'*-')
hold on
% plot(mobility./variablelist(1:kk),EL_maxint_5/max(EL_maxint_5),'o-')
plot(mobility./variablelist(1:kk),EL_maxint_2/max(EL_maxint_2),'o-')
plot(mobility./variablelist(1:kk),EL_maxint_1/max(EL_maxint_1),'o-')

plot(mobility./variablelist(1:kk),EL_maxint_05/max(EL_maxint_05),'o-')

hold off
ylim([0,1.1])
legend('maxPL','maxEL at 2mA','maxEL at 1mA','maxEL at 0.5mA')
 xlabel('mobility ( cm2 V-1s-1)')
 ax=gca;
 ax.XScale='log';
end
else
end