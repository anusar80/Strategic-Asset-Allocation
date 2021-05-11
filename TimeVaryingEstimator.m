%% Get data from python
data = readtable('data2.xlsx');
Week = data.Var1;
data = removevars(data, {'Var1'}); % , 'TEDRATE'
data = data(Week.Year>1990,:);
Week = Week(Week.Year>1990,:);
display(data(end-5:end,:))
%% Stocastically detrend the series
dataMatrix = table2array(data);
MA = nan(size(dataMatrix));
win = 12;
for num = 1:size(dataMatrix,1)-win
    for count = 1:size(dataMatrix,2)
        MA(num+win,count) = dataMatrix(num+win,count) - mean(dataMatrix(num:num+win-1,count),1,'omitnan');
    end
end
% Plot detrended series
for num=1:size(dataMatrix,2)
    figure()
    plot(Week,MA(:,num))
    title(data.Properties.VariableNames{num})
end
%% Time-invariant PLS
features = data(1,1:end-1).Properties.VariableNames;
idx = ~any(isnan(MA),2);
week = Week(idx);
window = 100;
cleanData = MA(idx,:);
X = cleanData(:,1:end-1);
y = cleanData(:,end);
[Xloadings,Yloadings,Xscore,Yscore,Beta,Explained,MSE] = plsregress(y, X, 1);
display(features)
display(Yloadings)
display(Explained)
display(MSE)
%
GlobalPLS = 0.6+0.1*normalize(table2array(data(:,1:end-1))*Yloadings);
figure()
area(Week(window:end),GlobalPLS(window:end))
ytickformat('%0.1f')
title('PLS-FCI global')
%% Expanding window loadings
window = 100;
n_obs = size(cleanData,1);
Loadings = nan(n_obs, size(cleanData,2)-1);
Rsqrd = nan(n_obs);

for num = window:n_obs
    X = cleanData(1:num,1:end-1);
    y = cleanData(1:num,end);
    [Xloadings,Yloadings,Xscore,Yscore,Beta,Explained,MSE] = plsregress(y, X, 1);
    Loadings(num,:) = Yloadings;
    Rsqrd(num) = Explained(2);
end
%% Fitted 
TimeVaryingPLS = nan(n_obs,1);
Data = data(idx,:);
for num=window:n_obs
    x = normalize(sum(table2array(Data(1:num,1:end-1)).*Loadings(1:num,:), 2), 'range');
    TimeVaryingPLS(num) = x(end);%cdf('Normal',x(end),0,1);
end
SAA = 1 - TimeVaryingPLS;
figure()
hold on
area(week(window:end),TimeVaryingPLS(window:end))
xline(datetime(2000,3,20));
xline(datetime(2001,9,11));
xline(datetime(2008,3,24));
xline(datetime(2008,9,12));
xline(datetime(2012,7,26));
xline(datetime(2020,2,23));
ytickformat('%0.1f')
ylabel('PLS-FCI expanding window')
title('The Global Financial Cycle')
saveas(gcf,'GFC.png')
%% Loadings over time
for num=1:size(Loadings,2)
    figure
    plot(Week,Loadings(:,num))
    title(features{num})
end
%% Securities held outright by the Fed
SecHeld = readtable('SecHeld');
GCF = TimeVaryingPLS;
QE = diff(SecHeld.SecHeld)/1000>50;
figure()
hold on
area(week(window:end),TimeVaryingPLS(window:end));
plot(SecHeld.Var1(2:end),QE,':','Color',[0.15 0.15 0.15])
title('QE and the global financial cycle')
axis tight
%% Strategic asset allocation
allocation = TimeVaryingPLS;
spxret = [nan; diff(data.SPX(idx))];
Bill = readtable('Bill.xlsx');
riskfree = (Bill.Bill(and(Bill.Var1>=week(window),Bill.Var1<=week(end)))/100+1).^(1/52)-1;
Portfolio = allocation(window:end).*spxret(window:end) + (1-allocation(window:end)).*riskfree;
Holding = 0.5*spxret(window:end) + 0.5*riskfree;
spx = spxret(window:end);
cumPort = ones(length(Portfolio),1);
cumHold = ones(length(Portfolio),1);
cumMarket = ones(length(Portfolio),1);

for num = 1:length(Portfolio)-1
    cumPort(num+1,1) = cumPort(num,1)*(1+Portfolio(num+1));
    cumHold(num+1,1) = cumHold(num,1)*(1+Holding(num+1));
    cumMarket(num+1,1) = cumMarket(num,1)*(1+spx(num+1));
end
%% Cumulative returns
figure()
hold on
plot(week(window:end),cumPort)
% plot(week(window:end),cumHold)
plot(week(window:end),cumMarket)
legend('strategic','hold SPX')
legend('location','best')
ytickformat('%0.1f')
title('Cumulative returns')
saveas(gcf,'cumRet.png')
%% Densities
[fh,xh] = ksdensity(spx); 
[fs,xs] = ksdensity(Portfolio); 
figure
hold on
plot(xh,fh)
plot(xs,fs)
legend('SPX','Strategic')
title('Kernel density of weekly returns')
saveas(gcf,'ksdensity.png')
%% Summary
figure()
scatter(spx, Portfolio,'+')
h=refline([1 0]);
h.Color = [0.15 0.15 0.15];
xline(0)
yline(0)
xlim([-0.162 0.083])
ylim([-0.162 0.083])
title('SPX and strategic portfolio returns')
ylabel('strategic')
xlabel('spx')
saveas(gcf,'compareStrategies.png')
%% Stack tables
Strategic = Portfolio;
SPX = spx;
Time = week(window:end);
Tbl = table(Time,SPX,Strategic);
Tbl = stack(Tbl,2:3);
[grp,mean,std,skew] = grpstats(Tbl.SPX_Strategic,Tbl.SPX_Strategic_Indicator,{'gname','mean','std',@skewness});
display(table(categorical(grp),(mean+1).^52-1,mean,std,skew,sqrt(52)*mean./std, 'VariableNames',{'Strategy', 'Mean (annualized)','Mean (weekly)','Std (weekly)', 'Skew','Sharpe ratio'}))
%% Export Strategic and SPX to excel
CompareRet = table(week(window:end),spx,Strategic,'VariableNames',{'Week','SPX','Strategic'});
writetable(CompareRet,'CompareRet.xlsx')