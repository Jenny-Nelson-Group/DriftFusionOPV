%% this is mohammed's code for running an unified model in OPV: EES2022
%% this model includes a molecular model, based on paper PRX2018, JACS2019, NC2021, a Kinetic model, and a drift diffusion model
%% using this code, we can run molecular model to calculate the rates of transition first, then kinetic model,
%% then drift diffusion model to get populations of CT and LE, as well as charge transport
%% modified by Jun Yan

clc
clear all
close all
%% defalut plot setting
%set(groot,'defaultfigureposition',[100 100 1200 600])
%set(0, 'DefaultLineLineWidth', 1.5);
set(gcf,'color','w')
set(gca,'FontSize',16)
set(gca,'linewidth',0.5)
%% 
Temp=linspace(80,300,10);
W=Temp';
fignumber=1;
EL_OUTPUT=zeros(501,1);
PL_OUTPUT=zeros(501,1);
i=0;
for var=Temp%
    
    i=i+1;

%% First Calculate the Recombination paramters of the CT and Exciton
%here we also calculate the absorption profile of the device and the
%emission spectra based on the properties of the states, 
Prec=paramsRec;%first initiliase the recombination parameters ( the parameters are set here by default and can be changed under or in the paramRec function
%below are a set of parameters that can be changed

%% input for LE and CT properties in the molecular model
Prec.params.tickness=100e-9;% active layer thickness in m
Prec.const.T=var;% temperature of the cell
Prec.params.RCTE=1e0;%ratio CT to S1
Prec.params.Excitondesnity=1e28;% exciton density 
Prec.params.CTdesnity=Prec.params.Excitondesnity*Prec.params.RCTE;% CT density 

%% EXciton properties          
Prec.params.Ex.hW=0.15;%main vibronic energy considered in eV
%Prec.params.Ex.hW2=0.012;%low frequency vibronic energy considered in eV
Prec.params.Ex.f=5;%oscillator strength of the CT state
Prec.params.Ex.L0=0.1;%outer reorganisation energy(High frequency) in eV
Prec.params.Ex.Li=0.1;%inner reorganisation energy in eV
Prec.params.Ex.DG0=1.5;%energy of the CT state in  eV
Prec.params.Ex.sigma=0.001;%gaussian distribution for CT state
Prec.params.Ex.numbrestate=3;
Prec.params.Ex.Dmus=3*3.33e-30/1.6e-19;%difference in static dipole moment (10 in DEbye ) then th
%% CT properties
Prec.params.CT.hW=0.15;%main vibronic energy considered in eV
%Prec.params.CT.hW2=0.012;%low frequency vibronic energy considered in eV
Prec.params.CT.f=1e-2;%oscillator strength of the CT state
Prec.params.CT.L0=0.1;%outer reorganisation energy(High frequency) in eV
Prec.params.CT.Li=0.1;%inner reorganisation energy in eV
Prec.params.CT.DG0=1.4;%energy of the CT state in  eV
Prec.params.CT.sigma=0.001;%gaussian distribution for CT state
Prec.params.CT.numbrestate=3;
Prec.params.CT.Dmus=3*3.33e-30/1.6e-19;%difference in static dipole moment (10 in DEbye ) then th
%%%%%%% add this to account for the effect of Hybredisation
Prec.params.Vstar=0.000; % Coupling between S1 and CT in eV

%% this line calculated the different properties of the states.
Prec=paramsRec.calcall(Prec);
krecCT=Prec.params.CT.results.knr+Prec.params.CT.results.krTot;
krecex=Prec.params.Ex.results.knr+Prec.params.Ex.results.krTot;
%Voc=Prec.results.Vocrad-Prec.results.Dvnr;

%% in this part we generate a device with a certain number of properties
% the default properties are set in pnParamsHCT or in the excel file p3hTPCBM.xlsx
NC=1e20;activelayer=2;Kfor=1e-11;%V in V, K in S-1, NC in Cm-3, Jsc in mA cm-2,
mobility=1e-3;kdis=1e11;kdisex=1e11;% Tq1 in s,mobility in Cm2V-1s-1
% to include thermal activated exciton diffusion
% d_domain=20;% nm
% Ld=(0.015*Prec.const.T+0.5);% nm
% ex_diff_eff=Ld/d_domain^2*(2*d_domain-Ld);
% kdisex=ex_diff_eff*kdisex;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deviceParameterFile='DeviceParameters_Default.xlsx';
DP = deviceparams(['parameters\',deviceParameterFile]);

DP.light_properties.OM=0;%to consider the transfer matrix generation profile
DP.Time_properties.tpoints=100;
DP.Layers{activelayer}.tp=Prec.params.tickness*100;% in cm
%DP=DP.generateDeviceparams(NC,activelayer,mobility,kdis,kdisex,Prec,1e-11,0);
%% temperature dependent exciton lifetime

%% 0D model
% PL
Gex=1e25;
DP=DP.generateDeviceparams(NC,activelayer,mobility,kdis,kdisex,Prec,Kfor,0);
[X_PL,Y_PL,Z_PL]=DP.simulate_PL(Prec,fignumber,Gex);
DP.simulate_PL_normalized(Prec,fignumber,Gex);
PL_PEAK(i)=max(Y_PL);
PL_OUTPUT=[PL_OUTPUT,Y_PL];
% EL
Ginj=1e25;
DP=DP.generateDeviceparams(NC,activelayer,mobility,kdis,kdisex,Prec,Kfor,0);
[X_EL,Y_EL,Z_EL]=DP.simulate_EL(Prec,fignumber,Ginj);
[X_EL_norm,Y_EL_norm,Z_EL_norm]=DP.simulate_EL_normalized(Prec,fignumber,Ginj);

EL_PEAK(i)=max(Y_EL);
%EL_PEAK_NEW(i)=max(NEW_EL);
EL_OUTPUT=[EL_OUTPUT,Y_EL];
%DP.simulateTAS(2e10,1e27,fignumber);
figure(fignumber)
subplot(3,4,1)
semilogy(Prec.const.Edistribution,Prec.results.AbsLJ)
xlabel('Energy [eV]')
ylabel('EQE [a.u]')
ylim([1*1e-7, 10])
xlim([0.5,2])
hold on
%% test reciprocity relation
[X,Y,Z]=simulate_EL(DP,Prec,fignumber,Ginj);% calculate EL
EQE=Prec.results.AbsLJ;
EQE_From_EL=Y./Prec.const.bb(:,2);
ind1=find(Prec.results.AbsLJ~=0);
ind=ind1(10);
Fit_EQE=EQE_From_EL/EQE_From_EL(ind)*EQE(ind);
figure(fignumber)
subplot(3,4,1)
% semilogy(X,Y/max(Y))
% hold on
semilogy(X,Fit_EQE,'k--','linewidth',1)
hold on

%% plot DVnr vs Voc,rad
figure (fignumber)
subplot(3,4,3)
plot(DP.results.Vocrad,DP.results.DVnr,'o','MarkerSize',10)
xlabel('V_{oc,rad} [V]')
ylabel('DV_{nr} (V)')
hold on 
% Prec.results.Dvnr
% DP.results.DVnr
% Prec.results.Vocrad
% DP.results.Vocrad
%% plot rate 
figure (fignumber)
subplot(3,4,2)
semilogy(Prec.const.T,Prec.params.Ex.results.krTot+Prec.params.CT.results.krTot,'o','MarkerSize',10)
hold on
semilogy(Prec.const.T,Prec.params.Ex.results.knr+Prec.params.CT.results.knr,'d','MarkerSize',10)
xlabel('Temp (T)')
ylabel('Kr or Knr (s^{-1})')
hold on 

%% to include FC dynamics
Vstart=0;Vend=5;
 DP.Layers{2}.r0=0;
    tic
    DP=DP.generateDeviceparams(NC,activelayer,mobility*exp(-(2/3*0.06/Prec.const.kb/Prec.const.T)),kdis,kdisex,Prec,Kfor,0);%*exp(-(2/3*0.06/Prec.const.kb/Prec.const.T))
    DV2=device(DP);
    DV2.Prec=Prec;
    toc
    
    G=1e-6;
    tic
    DV2=device.runsolJsc(DV2,G);
    toc
    
    tic
    DV2=device.runsolJV(DV2,G,Vstart,Vend);
    toc
    
%jv
figure (100)
dfplot.JV_new(DV2.sol_JV(1),1,0)
hold on
% PL
figure (fignumber)
subplot(3,4,11)
[~,PLplot,~]=dfplot.photoluminescence_mult(DV2,2,1,1e-6,num2str(var));
PL_max_FC(i)=max(PLplot);
hold on
figure (fignumber)
subplot(3,4,12)
semilogy(Prec.const.Edistribution,PLplot/max(PLplot))
            xlabel('Energy [eV]')
            ylabel('PL norm. [a.u]')
            ylim([1E-3, 10])
            xlim([0.5,2])
hold on 
% EL
figure (fignumber)
subplot(3,4,9)
[~,ELplot]=dfplot.Electroluminescence_multi(DV2,2,0,1,num2str(var),0);
EL_max_FC(i)=max(ELplot);
hold on 
figure (fignumber)
subplot(3,4,10)
semilogy(Prec.const.Edistribution,ELplot/max(ELplot))
            xlabel('Energy [eV]')
            ylabel('EL norm. [a.u]')
            ylim([1E-3, 10])
            xlim([0.5,2])
            hold on
end

%% plot state distribution
% CT
Prec.params.CT.DOS=1/Prec.params.CT.sigma/sqrt(2*pi).*exp(-1/2*((Prec.params.CT.Statedistribution-Prec.params.CT.DG0)/Prec.params.CT.sigma).^2);
Prec.params.CT.DOStrapz=trapz(Prec.params.CT.Statedistribution,Prec.params.CT.DOS);
% EX
Prec.params.Ex.DOS=1/Prec.params.Ex.sigma/sqrt(2*pi).*exp(-1/2*((Prec.params.Ex.Statedistribution-Prec.params.Ex.DG0)/Prec.params.Ex.sigma).^2);
Prec.params.Ex.DOStrapz=trapz(Prec.params.Ex.Statedistribution,Prec.params.Ex.DOS);
paramsRec.calcall(Prec);
figure (fignumber)
subplot(3,4,4)
semilogy(Prec.params.Ex.Statedistribution,Prec.params.Ex.DOS,'r')
hold on
semilogy(Prec.params.CT.Statedistribution,Prec.params.CT.DOS*Prec.params.RCTE,'k')
hold on
xlabel('Energy [eV]')
ylabel('DOS [a.u]')
% ylim([1e-3, 1000])
xlim([0.5,2.5])
legend('LE','CT')
%% show key variables
sgtitle(['offset = ' num2str(Prec.params.Ex.DG0-Prec.params.CT.DG0) '; ' 'fosc,CT = ' num2str(Prec.params.CT.f) ])


%% EL PL PEAK INTENSITY
figure (22)
semilogy(Temp,PL_PEAK/PL_PEAK(end),'o-')
hold on
semilogy(Temp,EL_PEAK/EL_PEAK(end),'o-')
hold on
semilogy(Temp,PL_max_FC/PL_max_FC(end),'ko-')
hold on
semilogy(Temp,EL_max_FC/EL_max_FC(end),'ro-')
hold on

%semilogy(Temp,EL_PEAK_NEW/EL_PEAK_NEW(end),'*-')
legend ('PL','EL','PL-FC','EL-FC')
%hold on
xlabel('Temperature (K)','fontsize',15)
ylabel('Peak intensity (norm.)','fontsize',15)
W=[W,PL_PEAK'/PL_PEAK(end),EL_PEAK'/EL_PEAK(end),PL_max_FC'/PL_max_FC(end),EL_max_FC'/EL_max_FC(end)];

EL_OUTPUT(:,1)=Prec.const.Edistribution';
PL_OUTPUT(:,1)=Prec.const.Edistribution';

%% normalize EL PL
for ii=2:length(EL_OUTPUT(1,:))
EL_OUTPUT_NORM(:,ii)=EL_OUTPUT(:,ii)/max(EL_OUTPUT(:,ii));
end
EL_OUTPUT_NORM(:,1)=EL_OUTPUT(:,1);
for jj=2:length(EL_OUTPUT(1,:))
PL_OUTPUT_NORM(:,jj)=PL_OUTPUT(:,jj)/max(PL_OUTPUT(:,jj));
end
PL_OUTPUT_NORM(:,1)=PL_OUTPUT(:,1);

