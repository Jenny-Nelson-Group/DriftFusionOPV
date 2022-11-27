% AL_thickness=200*1e-7;%active layer tickness
kdis0=DV_D4A2_nofield.DP.Layers{2}.kdis;
kfor0=DV_D4A2_nofield.DP.Layers{2}.kfor;
counter=1;
for kdis=logspace(9,11.5,8)
    counter=counter+1;
    Prec=DV_D4A2_nofield.Prec;
    %% in this part we generate a device with a certain number of properties
    %the default properties are set in pnParamsHCT or in the excel file p3hTPCBM.xlsx
    
    NC=2e19*sqrt(kdis/kdis0);activelayer=2;Kfor=kfor0;%V in V, K in S-1, NC in Cm-3, Jsc in mA cm-2,
    mobility=3e-4;kdisex= DV_D4A2_nofield.DP.Layers{2}.kdisexc;% Tq1 in s,mobility in Cm2V-1s-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DP=pnParamsHCT;
    
    DP.light_properties.OM=0;%to consider the transfer matrix generation profile
    DP.Time_properties.tpoints=100;
    DP=DP.generateDeviceparams(NC,activelayer,mobility,kdis,kdisex,Prec,Kfor,0);
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
    suns=[1e-6,1];
    for Gen=suns
        tic
        DV2=device.runsolJsc(DV2,Gen);
        toc
        
        tic
        DV2=device.runsolJV(DV2,Gen,Vstart,Vend);
        toc
    end
    assignin('base',"DV_D4A2_fixBfor"+num2str(counter),DV2)
    %%
    figure(2)
    dfplot.JV_new(DV2.sol_JV(1),1)
    [Jsc(counter),Voc(counter),FF(counter)]=dfplot.JV_new(DV2.sol_JV(2),1);
    figure(3)
    dfplot.photoluminescence_mult(DV2,2,0,1,num2str(counter)+" K");
    figure(4)
    dfplot.Electroluminescence_multi(DV2,2,0,2,num2str(counter)+" K")
end


%%
counter=2;
for kdis=1:6
        eval("DV2=DV_D4A2_fixBfor"+num2str(counter)+";")
        counter=counter+1; 
[GSBtableD4(:,1),GSBtableD4(:,counter)]=DV2.DP.simulateTAS(2e10,1e27,fignumber);
figure(3)
[JVcat(3,counter),JVcat(4,counter),JVcat(5,counter),JVmodtableD4(:,counter),JVmodtableD4(:,1)]=dfplot.JV_new(DV2.sol_JV(2),1);
JVcat(1,counter)=DV2.DP.Layers{2}.kdis;
JVcat(2,counter)=DV2.DP.Layers{2}.kfor;
JVcat(7,counter)=DV2.DP.Layers{2}.N0C;
end
JVcat=JVcat';
JVcat(:,6)=JVcat(:,3).*JVcat(:,4).*JVcat(:,5);