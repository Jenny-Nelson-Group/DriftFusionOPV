Hab_CT_dis=0.004;
L0_CT_dis=0.11;
Hab_LE_dis=0.03;
L0_LE_dis=0.54;
marcusrate =  @(Hab,L0,x) Hab^2*exp(-(-x+L0).^2/0.026/4/L0)/sqrt(4*pi*L0*0.026)/6.62e-34*1.6e-19;
NC=2e19;activelayer=2;mobility=3e-4;
devicenames={'D0A2','D2A2','D4A2'};


for kk=1:3
    devicename=devicenames{kk};
    eval("DV2=DV_"+devicename+"_nofield;");
%     Gen=1e-6;
%                 tic
%             DV2=device.runsolJsc(DV2,Gen);
%             toc
%             
%             tic
%             DV2=device.runsolJV(DV2,Gen,Vstart,Vend);
%             toc
    
    for V=2:2:4
       DV2=device.current_transient(DV2,1e-6,0.4,V,1);
    end
     eval("DV_"+devicename+"_nofield=DV2;");
end

%%
close all
dfplot.transient_EL(DV_D0A2_nofield,DV_D0A2_nofield.sol_Vpulse(2),2,0,"D0")
dfplot.PL_T(DV_D2A2_nofield,DV_D2A2_nofield.sol_Vpulse(2),2,1,"D2")
dfplot.PL_T(DV_D4A2_nofield,DV_D4A2_nofield.sol_Vpulse(2),2,1,"D4")
dfplot.PL_T(DV_D0A2_nofield,DV_D0A2_nofield.sol_Vpulse(end),2,0,"D0")
dfplot.PL_T(DV_D2A2_nofield,DV_D2A2_nofield.sol_Vpulse(end),2,1,"D2")
dfplot.PL_T(DV_D4A2_nofield,DV_D4A2_nofield.sol_Vpulse(end),2,1,"D4")


%%
close all
dfplot.transient_EL(DV_D0A2_nofield,DV_D0A2_nofield.sol_Vpulse(7),2,0,"D0")
dfplot.transient_EL(DV_D2A2_nofield,DV_D2A2_nofield.sol_Vpulse(end),2,1,"D2")
dfplot.transient_EL(DV_D4A2_nofield,DV_D4A2_nofield.sol_Vpulse(end),2,1,"D4")
dfplot.transient_EL(DV_D0A2_nofield,DV_D0A2_nofield.sol_Vpulse(4),2,0,"D0")
dfplot.transient_EL(DV_D2A2_nofield,DV_D2A2_nofield.sol_Vpulse(4),2,1,"D2")
dfplot.transient_EL(DV_D4A2_nofield,DV_D4A2_nofield.sol_Vpulse(4),2,1,"D4")
%%
close all
dfplot.Electrondistrbution(DV_D0A2_nofield,DV_D0A2_nofield.sol_Vpulse(end-1),[0,0.5e-7,0.8e-7,1e-7,1.5e-7,2e-7,3e-7,4e-7])
dfplot.Electrondistrbution(DV_D4A2_nofield,DV_D4A2_nofield.sol_Vpulse(end-1),[0,0.5e-7,0.8e-7,1e-7,1.5e-7,2e-7,3e-7,4e-7])