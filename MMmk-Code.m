clear all;close all;clc;

% % % % % Initializations
noQueueVector = [];srvingState = [];arrTotalSysTime = [];arrWaitTime = [];noSysVector = [];
mu = 1;m = 1;k = 100;lambda = 0.9*m;
noArr = 500000; d = 1; iArr = 1; iState = 0; state = 0; firstArrived = 0;j=1;

% % % % % Arrival Times
step = 10^-d; tolerance = 10^-(d+1);rho = lambda/mu;
arrTimes = round((exprnd((1/lambda),1,noArr)+(10^-d)),d);
arrTimeline = cumsum(arrTimes);
timeLine = 0:step:ceil(arrTimeline(end));
server(1:m) = 0;check = [];

disp('Simulation Running............');
% % % % % Simulation Starts

for i = 0:step:ceil(arrTimeline(end))
    iArr
    srvingState = srvingState - step;
    server = server - step;
    
    for iServer=1:m
        if(server(iServer)<0)
            server(iServer) = 0;
        end
    end
    
    monitor(iState+1,:) = server;

    if(iArr>length(arrTimeline))
        iArr = length(arrTimeline);
    end
    
    if(abs(i-arrTimeline(iArr))<tolerance)
        state = state + 1;
        [minVal, minIndex] = min(server);
        
        if(state>k)
            state = state - 1;
            arrWaitTime(iArr) = inf;
            arrTotalSysTime(iArr) = inf;
            srvingState(iArr) = -1;
        else
            iArrSrvTime = round((exprnd((1/mu))+(10^-d)),d);
            arrSrvTime(iArr) = iArrSrvTime;
            server(minIndex) = server(minIndex) + iArrSrvTime;
            arrWaitTime(iArr) = minVal;
            arrTotalSysTime(iArr) = iArrSrvTime + arrWaitTime(iArr);
            srvingState(iArr) = arrTotalSysTime(iArr); 
        end
        iArr = iArr + 1;
    end
        
    if(~isempty(find(srvingState<tolerance & srvingState>-tolerance, 1)))
        state = state - sum(srvingState<tolerance & srvingState>-tolerance);
    end
    
    iState = iState + 1;
    noSysVector(iState) = state;
    check(iState) = ~isempty(find(srvingState<tolerance & srvingState>-tolerance, 1));

    
    if(state>m)
        noQueueVector(iState) = state - m;
    else
        noQueueVector(iState) = 0;
    end
    
end
disp('Simulation Ended.....');
% % % % % Simulation Ends

% plot(arrSrvTime,'--r');hold on;plot(arrTimes);hold off;

figure;stairs(timeLine,noSysVector);
title(['m = ' num2str(m) ' and k = ' num2str(k)]);
xlabel('t - time axis');ylabel('N(t) - no. Packets in system');

% % % % % Time Averages
avgExpNoSystem = trapz(timeLine,noSysVector)/timeLine(end);
avgExpNoQueue = trapz(timeLine,noQueueVector)/timeLine(end);
avgExpWaitTime = mean(arrWaitTime(arrWaitTime<inf));
avgExpSysTime = mean(arrTotalSysTime(arrTotalSysTime<inf));


% % % % % Statistical Averages

% M/M/m/k system
p0 = 1./(sum((rho.^(0:(m-1)))./factorial(0:(m-1)))+((rho^m/factorial(m))*((1-((rho/m)^(k-m+1)))/(1-(rho/m)))));
pk = (rho^k)*p0/(factorial(m)*(m^(k-m)));
avgThryNoQueue = (rho^(m+1)*p0/factorial(m-1))*((1-(rho/m)^(k-m+1)-((1-rho/m)*(k-m+1)*((rho/m)^(k-m))))/(m-rho)^2);
avgThryNoSystem = avgThryNoQueue+rho;
lambda0 = lambda*(1-pk);
avgThryWaitTime = avgThryNoQueue/lambda0;
avgThrySysTime = avgThryNoSystem/lambda0;

% M/M/1/inf system uncomment the vaules below if used for mm1inf
% avgThryNoSystem = rho/(1-rho);
% avgThryNoQueue = avgThryNoSystem-rho;
% avgThryWaitTime = avgThryNoQueue/lambda;
% avgThrySysTime = avgThryNoSystem/lambda;

Result = [avgExpNoSystem, avgExpNoQueue, avgExpSysTime, avgExpWaitTime;...
          avgThryNoSystem,avgThryNoQueue, avgThrySysTime, avgThryWaitTime];  
      
      
