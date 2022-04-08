function fitresult=FitPL_TEMP_test(BE,Templist)
% code to fit the temp PL using a single state
% get the PLs in the temp list
% [1.4,0.12,0.08,0.15,0.022] data for Y6:PS
counter=0;
for kk=1:length(BE.PLs)
    if max(str2num((BE.PLs(kk).temperature))==Templist)
        counter=counter+1;
        Energy_exp(:,counter)=1240./BE.PLs(kk).spectra.wavelength_nm;
        Signal_exp(:,counter)=BE.PLs(kk).spectra.Correctedsignal_Ev;
        Temp(counter)=str2double(BE.PLs(kk).temperature);
    end
end

%% generating the data to fit to, and clean it
centerenergy=1.3;
Energy=centerenergy-0.5:0.01:centerenergy+0.2;

yData=0;
xData=0;
for kk=1:size(Signal_exp,2)
    signal0=interp1(Energy_exp(:,kk),Signal_exp(:,kk),Energy);
    %     signal0=signal0./signal0*max(signal0);
    if(xData==0)
        xData=Energy;
        yData=signal0;
    else
        xData=[xData,Energy];
        yData=[yData,signal0];
    end
end
yData=yData./max(max(yData));
yData(isnan(yData))=0;
%%
Prec=paramsRec;%first initiliase the recombination parameters ( the parameters are set here by default and can be changed under or in the paramRec function
%below are a set of parameters that can be changed
StartPoint(1)=0.1;StartPoint(2)=0.1;StartPoint(3)=0.02;%StartPoint(4)=1.34;StartPoint(5)=0.01;
Upper(1)=0.2;Upper(2)=0.2;Upper(3)=0.05;%Upper(4)=1.5;Upper(5)=0.03;
Lower(1)=0.05;Lower(2)=0.05;Lower(3)=0.005;%Lower(4)=1;Lower(5)=0.001;
function_text="@(x,xdata) simulateemission(Prec,Temp,xdata(1:length(Energy)),x(1),x(2),0.145,1.38,x(3));";
Prec.params.Ex.DG0=1.38;%StartPoint(1);
Prec.params.Ex.Li=StartPoint(1);
Prec.params.Ex.L0=StartPoint(2);
Prec.params.Ex.f=5;
Prec.params.Ex.hW=0.15;%StartPoint(3);
Prec.params.Ex.sigma=StartPoint(3);
Prec.params.Ex.numbrestate=5;
Prec = paramsRec.update(Prec);
Prec.params.Ex=Prec.FCWD(Prec.params.Ex,Prec.const); 
Prec.params.Ex=Prec.Calcrate(Prec.params.Ex,Prec.const);

%% manual fitting of the data and setting up the free parameters for the fit 
move_to_fit='N';
    function_text="@(x,xdata) simulateemission(Prec,Temp,xdata(1:length(Energy)),x(2),x(3),x(4),x(1),x(5));";
    eval("fun="+function_text);
while move_to_fit=='N'
    prompt = "input the 5 parameter of the state of interrest in the format [DG0,Li,L0,hw,sigma] ";
    data_input= input(prompt);

    results=fun(data_input,xData(1:length(Energy)));
    % plot the results
    figure (1)
    subplot(1,2,1)
    results=reshape(results,[length(Energy),size(Signal_exp,2)]);
    plot(Energy,results);
    title('simulated results')
    lgd = legend;
    lgd.String=num2str(Temp');
    
    subplot(1,2,2)
    yData=reshape(yData,[length(Energy),size(Signal_exp,2)]);
    plot(Energy,yData);
    title('experiemental results')
    lgd = legend;
    lgd.String=num2str(Temp');

    
        figure (2)
        count=0;
        for TT=Temp
            count=count+1;
            subplot(3,2,count)
            plot(Energy,results(:,count)./max(results(:,count)),'.-','DisplayName',"fit");
            hold on
            plot(Energy,yData(:,count)./max(yData(:,count)),'DisplayName',"experimental");
            title([num2str(TT) ' K'])
            legend()
            hold off
        end
        prompt = "Do you want to move to the fit ? Y/N [Y]: ";
        move_to_fit = input(prompt,"s");
end

%%
StartPoint(1)=0.1;StartPoint(2)=0.1;StartPoint(3)=0.02;%StartPoint(4)=1.34;StartPoint(5)=0.01;
Upper(1)=0.2;Upper(2)=0.2;Upper(3)=0.05;%Upper(4)=1.5;Upper(5)=0.03;
Lower(1)=0.05;Lower(2)=0.05;Lower(3)=0.005;%Lower(4)=1;Lower(5)=0.001;
function_text="@(x,xdata) simulateemission(Prec,Temp,xdata(1:length(Energy)),x(1),x(2),0.145,1.38,x(3));";
results=fun(StartPoint,xData(1:length(Energy)));

%% Set up fittype and options.
option = optimoptions('lsqcurvefit');
option.Display = 'iter';
option.MaxIterations=5;
[fitresult,resnorm,resid,exitflag,output,lambda,J]= lsqcurvefit(fun,StartPoint ,xData,yData,Lower,Upper,option );

results=fun(fitresult,xData(1:length(Energy)));
%% plotting the results
figure (1)
subplot(1,2,1)
results=reshape(results,[length(Energy),size(Signal_exp,2)]);
plot(Energy,results);
title('simulated results')
lgd = legend;
lgd.String=num2str(Temp');

subplot(1,2,2)
yData=reshape(yData,[length(Energy),size(Signal_exp,2)]);
plot(Energy,yData);
title('experiemental results')
lgd = legend;
lgd.String=num2str(Temp');
figure (2)
count=0;
for TT=Temp
    count=count+1;
    subplot(3,2,count)
    plot(Energy,results(:,count)./max(results(:,count)),'.-','DisplayName',"fit");
    hold on
    plot(Energy,yData(:,count)./max(yData(:,count)),'DisplayName',"experimental");
    title([num2str(TT) ' K'])
    legend()
    hold off
end

end

function results=simulateemission(Prec,Temp_list,Edistribution,Li_LE,L0_LE,hW_LE,E_LE,Sig_LE)
results=[];
for Temp=Temp_list
Prec.const.T=Temp;
% Prec.params.CT.DG0=Prec.params.Ex.DG0-offset;
Prec.params.Ex.Li=Li_LE;
Prec.params.Ex.L0=L0_LE;
Prec.params.Ex.hW=hW_LE;
Prec.params.Ex.DG0=E_LE;
Prec.params.Ex.sigma=Sig_LE;
Prec.const.Edistribution=Edistribution;
Prec = paramsRec.update(Prec);
Prec.params.Ex=Prec.FCWD(Prec.params.Ex,Prec.const); 
Prec.params.Ex=Prec.Calcrate(Prec.params.Ex,Prec.const);

krE=Prec.params.Ex.results.krE./(Prec.params.Ex.results.krTot+Prec.params.Ex.results.knr);%assuming a fixed generation rate of exciton
X=(Prec.const.Edistribution)';
Y=(krE)';

results=[results,interp1(X,Y,Edistribution)];
end
results=results/max(max(results));
end