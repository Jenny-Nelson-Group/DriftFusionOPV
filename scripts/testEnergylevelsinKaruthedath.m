Hab_CT_dis=0.004;
L0_CT_dis=0.11;
Hab_LE_dis=0.03;
L0_LE_dis=0.54;
marcusrate =  @(Hab,L0,x) Hab^2*exp(-(-x+L0).^2/0.026/4/L0)/sqrt(4*pi*L0*0.026)/6.62e-34*1.6e-19;
NC=2e19;activelayer=2;mobility=3e-4;
devicenames={'D0A2','D2A2','D4A2'};
kk=2;
devicename=devicenames{kk};
% eval(['expresults=' devicename ';']);

eval(['Prec=DV_' devicename '_nofield.Prec;']);    % close all
Prec.params.CT.numbrestate=1;
Prec.params.Ex.numbrestate=1;

%%



offsets=[];

Jsc=[];
Voc=[];
FF=[];
index2=0;
for LiCT=0.07:0.03:0.16
    index=0;
    index2=index2+1;
    for offset=0:0.05:0.2
        %% this part for the energy levels all energies in EV
        EA_A=-3.85-offset;
        EA_D=-3.1;
        IP_A=-5.75;
        IP_D=-5.1;
        E_Coul_CT=0.28;
        E_Coul_Ex=0.27;
        B=0.37;%quadrupole moment effect
        E_Ex=(min(EA_A-IP_A-E_Coul_Ex,EA_D-IP_D-E_Coul_Ex));
        E_CT=(min(EA_A-IP_D-E_Coul_CT+B,EA_D-IP_A-E_Coul_CT+B));
        ECS=-(max(IP_A,IP_D)-min(EA_A,EA_D));
        %%
        Prec.params.CT.DG0=E_CT;
        Prec.params.Ex.DG0=E_Ex;
        Prec.params.CT.Li=0.07;
        Prec.params.Ex.Li=LiCT;
        index=index+1;
        kdis_LE=marcusrate(Hab_LE_dis,L0_LE_dis,E_Ex-E_CT);
        offset_CTCS=Prec.params.CT.DG0-ECS;
        kdis_CT=marcusrate(Hab_CT_dis,L0_CT_dis,offset_CTCS);
        Prec=paramsRec.calcall(Prec);
        krecCT=Prec.params.CT.results.knr;krecex=Prec.params.Ex.results.knr;
        DP=pnParamsHCT;
        DP.Time_properties.tpoints=100;
        DP=DP.generateDeviceparams(NC,activelayer,mobility,kdis_CT,kdis_LE,Prec,ECS,2);
        fignumber=kk;
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
        %%
        Vstart=0;Vend=1.5;
        tic
        DP.Layers{2}.r0=0; %R0 for field dependence is 1 nm
        
        DV2=device(DP);
        DV2.Prec=Prec;
        
        %%
        suns=[1];
        for Gen=suns
            tic
            DV2=device.runsolJsc(DV2,Gen);
            toc
            
            tic
            DV2=device.runsolJV(DV2,Gen,Vstart,Vend);
            toc
        end
        offsets(index)=offset;
        figure(3)
        [Jsc1,Voc1,FF1]=dfplot.JV_new(DV2.sol_JV(1),1);
        Jsc(index,index2)=Jsc1;
        Voc(index,index2)=Voc1;
        FF(index,index2)=FF1;
        figure(22)
        subplot(1,3,1)
        plot(offsets,Jsc(:,index2),'*')
        lg=legend;
        lg.String{end}="CT Li= "+num2str(LiCT)+" eV";
        subplot(1,3,2)
        plot(offsets,Voc(:,index2),'*')
        subplot(1,3,3)
        plot(offsets,FF(:,index2),'*')
        pause(0.1)
        devicename=devicenames{kk}+"_"+num2str(index)+"_"+num2str(index2);
        
        eval("DV_"+devicename+"=DV2;");
        
        
    end
end

%%
kk=1;
index=0;
eval(['ECS=-DV_' devicenames{kk} '_nofield.DP.Layers{2}.IP;']);
offsets=zeros(8,1);
Krrec_CT=zeros(8,1);Kdis_Ex=zeros(8,1);
Kdis_CT=zeros(8,1);Kfor=zeros(8,1);offset_CTCS=zeros(8,1);Knrec_CT=zeros(8,1);
Krrec_Ex=zeros(8,1);Knrec_Ex=zeros(8,1);Voc_est=zeros(8,1);
% Jsc=zeros(8,1);Voc=zeros(8,1);FF=zeros(8,1);

index2=0;
for LiCT=0.07:0.03:0.16
    index=0;
    index2=index2+1;
    for offset=0:0.05:0.2%[0.1:0.04:0.4 0.02 0.06]%0.1:0.04:0.4
    index=index+1;
    devicename=devicenames{kk}+"_"+num2str(index);
    
    eval("DV2=DV_"+devicename+";")
    offsets(index)=0:0.05:0.2;
    
    %     figure(kk+1)
    %     [Jsc1,Voc1,FF1]=dfplot.JV_new(DV2.sol_JV(1));
    %     Jsc(index)=Jsc1;
    %     Voc(index)=Voc1;
    %     FF(index)=FF1;
    figure(kk+2)
    subplot(2,3,1)
    plot(offsets,Jsc,'*-')
    xlabel('Energy offset (LE-CT) [eV]')
    ylabel('Jsc [mA cm-2]')
    subplot(2,3,2)
    Voc_est(index)=DV2.Prec.results.Vocrad-DV2.Prec.results.Dvnr;
    plot(offsets,Voc,'*-')
    hold on

    xlabel('Energy offset (LE-CT) [eV]')
    ylabel('Voc [V]')
    subplot(2,3,3)
    plot(offsets,FF,'*')
    xlabel('Energy offset (LE-CT) [eV]')
    ylabel('FF [%]')
    pause(0.1)
    
    Knrec_CT(index)=DV2.Prec.params.CT.results.knr;
    Krrec_CT(index)=DV2.Prec.params.CT.results.krTot;
    Knrec_Ex(index)=DV2.Prec.params.Ex.results.knr;
    Krrec_Ex(index)=DV2.Prec.params.Ex.results.krTot;
    Kdis_CT(index)=DV2.DP.Layers{2}.kdis;
    Kdis_Ex(index)=DV2.DP.Layers{2}.kdisexc;
    Kfor(index)=DV2.DP.Layers{2}.kfor;
    offset_CTCS(index)=DV2.Prec.params.CT.DG0-ECS;
    subplot(2,3,4)
    plot(offsets,Jsc.*Voc.*FF,'*')
    hold on
    plot(offsets,Jsc.*Voc_est.*FF,'*')
    hold off
    legend("Voc drift diffusion","Voc estimated from the nrad model")
    
    xlabel('Energy offset (LE-CT) [eV]')
    ylabel('PCE [%]')
    subplot(2,3,5)
    [~,b]=sort(offsets);
    semilogy(sort(offsets),Knrec_CT(b),'*-')
    hold on
    %       semilogy(sort(offsets),Krrec_CT(b),'*-')
    semilogy(sort(offsets),Kdis_Ex(b),'*-')
    semilogy(sort(offsets),Kdis_CT(b),'*-')
    semilogy(sort(offsets),Knrec_Ex(b),'*-')
    %       semilogy(sort(offsets),Krrec_Ex(b),'*-')
    hold off
    xlabel('Energy offset (LE-CT) [eV]')
    ylabel('Dissociation and recombination rate constant [s-1]')
    legend("CT nr recombination","LE dissociation","CT dissociation","LE nr recombination")
    
    subplot(2,3,6)
    semilogy(offsets,Kfor,'*')
    xlabel('Energy offset (LE-CT) [eV]')
    ylabel('Reformation rate constant [cm3s-1]')
    %     plot(offsets,offset_CTCS,'*')
    end
end
sgtitle("calculation based on "+devicenames{kk}+" where E_C_S="+num2str(ECS)+" eV and E_L_E="+num2str(DV2.Prec.params.Ex.DG0+"eV"));


%%
kk=2;
index2=0;
for LiCT=0.07:0.03:0.16
    index=0;
    index2=index2+1;
    for offset=0:0.05:0.2
         index=index+1;
        devicename=devicenames{kk}+"_"+num2str(index)+"_"+num2str(index2);
        devicename=char(devicename);
        eval(['results.ECS(index,index2)=-DV_' devicename '.DP.Layers{2}.IP;']);
        eval(['results.DGLE(index,index2)=DV_' devicename '.Prec.params.Ex.DG0;']);
        eval(['results.DGCT(index,index2)=DV_' devicename '.Prec.params.CT.DG0;']);
        results.offset_LE_CT=results.DGLE-results.DGCT;
        results.offset_CT_CS=results.DGCT-results.ECS;
        eval(['results.kdisCT(index,index2)=DV_' devicename '.DP.Layers{2}.kdis;']);
        eval(['results.kdisLE(index,index2)=DV_' devicename '.DP.Layers{2}.kdisexc;']);
        eval(['results.knrecCT(index,index2)=DV_' devicename '.Prec.params.CT.results.knr;']);
        eval(['results.knrecLE(index,index2)=DV_' devicename '.Prec.params.Ex.results.knr;']);
        eval(['results.Kfor(index,index2)=DV_' devicename '.DP.Layers{2}.kfor;']);
        results.EA_A(index,index2)=-3.85-offset;
        results.offsets=offsets;
        results.LiCT=0.07:0.03:0.16;
        results.FF=FF;
        results.Jsc=Jsc;
        results.Voc=Voc;
        results.PCE=FF.*Voc.*Jsc;
    end
end

%%
close all
offset=0:0.01:0.4;
fun=@(x,xdata) x(1)^2*exp(-(-xdata+x(2)).^2/0.026/4/x(2))/sqrt(4*pi*x(2)*0.026)/6.62e-34*1.6e-19;
x=lsqcurvefit(fun,[0.01,0.5],offset_LE_CT(2:3),kdisLE(2:3));
figure
plot(offset_LE_CT,kdisLE,'*',offset,fun(x,offset))


%%
figure(22)
subplot(1,3,1)
plot(offsets(2:4),Jsc(2:4),'*')
subplot(1,3,2)
plot(offsets(2:4),Voc(2:4),'*')
subplot(1,3,3)
plot(offsets(2:4),FF(2:4),'*')
pause(0.1)

%%
index2=0;
for LiCT=0.07:0.03:0.16
    index=0;
    index2=index2+1;
    figure(22)
        subplot(1,3,1)
        
        plot(offsets,Jsc(:,index2),'*-')
        lg=legend;
        lg.String{end}="CT Li= "+num2str(LiCT)+" eV";
        hold on
        subplot(1,3,2)
        plot(offsets,Voc(:,index2),'*-')
        hold on
        
        subplot(1,3,3)
        plot(offsets,FF(:,index2),'*-')
        hold on
        pause(0.1)
end