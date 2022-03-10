% AL_thickness=200*1e-7;%active layer tickness

%% First Calculate the Recombination paramters of the CT and Exciton
%here we also calculate the absorption profile of the device and the
%emission spectra based on the properties of the states, 
% 
%% do not forget to install the symbolic math package*
Temp=300;
AL_thickness=50*1e-7;
Prec=paramsRec;%first initiliase the recombination parameters ( the parameters are set here by default and can be changed under or in the paramRec function
%below are a set of parameters that can be changed
Prec.params.Ex.DG0=1.63;
Prec.params.CT.f=5e-2;
offset=0.3;
Prec.params.tickness=AL_thickness*1e-2;%in m
Prec.params.CT.DG0=Prec.params.Ex.DG0-offset;
Prec.params.RCTE=1e-1;
Prec.params.Ex.Li=0.15;
Prec.params.CT.Li=0.15;
Prec.params.CT.L0=0.18;
Prec.params.Ex.f=5;
Prec.params.Excitondesnity=1e28;
Prec.params.Vstar=0.001;
Prec.const.T=300;
Prec.params.CT.sigma=0.04;
Prec.params.Ex.sigma=0.02;
Prec.params.Ex.numbrestate=1;
Prec.params.CT.numbrestate=1;
Prec=paramsRec.calcall(Prec);%this line calculated the different properties of the states.
krecCT=Prec.params.CT.results.knr;krecex=Prec.params.Ex.results.knr;
Voc=Prec.results.Vocrad-Prec.results.Dvnr;
semilogy(Prec.const.Edistribution,Prec.results.AbsLJ)
hold on
xlabel('Energy [eV]')
ylabel('EQE [a.u]')
ylim([1*1e-7, 1])
xlim([1,2])
sgtitle(['Vocrad= ' num2str(Voc) ' Jscrad= ' num2str(Prec.results.Jscrad)])
%             legend("ToT","From CT","FromEx")
pause(0.1)
%% in this part we generate a device with a certain number of properties
%the default properties are set in pnParamsHCT or in the excel file p3hTPCBM.xlsx

NC=2e19;activelayer=2;Kfor=1e-11;%V in V, K in S-1, NC in Cm-3, Jsc in mA cm-2,
mobility=3e-4;kdis=1.3e10;kdisex=1e11;% Tq1 in s,mobility in Cm2V-1s-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DP=deviceparams('PINDevice.xlsx');

DP.light_properties.OM=0;%to consider the transfer matrix generation profile
DP.Time_properties.tpoints=100;
DP.Layers{activelayer}.tp=AL_thickness;
DP=DP.generateDeviceparams(NC,activelayer,mobility,kdis,kdisex,Prec,1e-11,0);
clear NC activelayer Tq1exp mobility kdis kdisex
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
tic
DP.Layers{2}.r0=0; %R0 for field dependence is 1 nm

DV2=device(DP);
DV2.Prec=Prec;
toc
%%
suns=[1];
for Gen=suns
    tic
    DV2=device.runsolJsc(DV2,Gen);
    DV2=device(DP);
    toc
    
    tic
    DV2=device.runsolJV(DV2,Gen,Vstart,Vend);
    toc
end
assignin('base',"DV_"+num2str(Temp)+"K",DV2)
%%
figure(2)
dfplot.JV_new(DV2.sol_JV(2),1)
hold on
dfplot.JV_new(DV2.sol_JV(end),1)
% figure(3)
% dfplot.photoluminescence_mult(DV2,2,0,1,num2str(Temp)+" K");
% figure(4)
% dfplot.Electroluminescence_multi(DV2,2,0,2,num2str(Temp)+" K")
