Neutral = readtable('Neutral.csv');

indexNeutralNo = table2array(Neutral(:, 'Dispersal')) == 0;
indexNeutralMed = table2array(Neutral(:, 'Dispersal')) == 2e-05;
indexNeutralFull = table2array(Neutral(:, 'Dispersal')) == 0.8647;
NeutralNoExp = mle(table2array(Neutral(indexNeutralNo, 'WaitTime')), 'Distribution', 'exp');
NeutralMedExp = mle(table2array(Neutral(indexNeutralMed, 'WaitTime')), 'Distribution', 'exp');
NeutralFullExp = mle(table2array(Neutral(indexNeutralFull, 'WaitTime')), 'Distribution', 'exp');
NeutralNoGam = mle(table2array(Neutral(indexNeutralNo, 'WaitTime')), 'Distribution', 'gamma');
NeutralMedGam = mle(table2array(Neutral(indexNeutralMed, 'WaitTime')), 'Distribution', 'gamma');
NeutralFullGam = mle(table2array(Neutral(indexNeutralFull, 'WaitTime')), 'Distribution', 'gamma');

Regional = readtable('Regional.csv');

indexRegionalNo = table2array(Regional(:, 'Dispersal')) == 0;
indexRegionalMed = table2array(Regional(:, 'Dispersal')) == 2e-05;
indexRegionalFull = table2array(Regional(:, 'Dispersal')) == 0.8647;
RegionalNoExp = mle(table2array(Regional(indexRegionalNo, 'WaitTime')), 'Distribution', 'exp');
RegionalMedExp = mle(table2array(Regional(indexRegionalMed, 'WaitTime')), 'Distribution', 'exp');
RegionalFullExp = mle(table2array(Regional(indexRegionalFull, 'WaitTime')), 'Distribution', 'exp');
RegionalNoGam = mle(table2array(Regional(indexRegionalNo, 'WaitTime')), 'Distribution', 'gamma');
RegionalMedGam = mle(table2array(Regional(indexRegionalMed, 'WaitTime')), 'Distribution', 'gamma');
RegionalFullGam = mle(table2array(Regional(indexRegionalFull, 'WaitTime')), 'Distribution', 'gamma');

Local = readtable('Local.csv');

indexLocalNo = table2array(Local(:, 'Dispersal')) == 0;
indexLocalMed = table2array(Local(:, 'Dispersal')) == 2e-05;
indexLocalFull = table2array(Local(:, 'Dispersal')) == 0.8647;
LocalNoExp = mle(table2array(Local(indexLocalNo, 'WaitTime')), 'Distribution', 'exp');
LocalMedExp = mle(table2array(Local(indexLocalMed, 'WaitTime')), 'Distribution', 'exp');
LocalFullExp = mle(table2array(Local(indexLocalFull, 'WaitTime')), 'Distribution', 'exp');
LocalNoGam = mle(table2array(Local(indexLocalNo, 'WaitTime')), 'Distribution', 'gamma');
LocalMedGam = mle(table2array(Local(indexLocalMed, 'WaitTime')), 'Distribution', 'gamma');
LocalFullGam = mle(table2array(Local(indexLocalFull, 'WaitTime')), 'Distribution', 'gamma');

tiledlayout(3, 3)
nexttile
myplot(table2array(Neutral(indexNeutralNo, 'WaitTime')),...
    NeutralNoExp, NeutralNoGam, 'Neutral, No Dispersal', 1)
nexttile
myplot(table2array(Neutral(indexNeutralMed, 'WaitTime')),...
    NeutralMedExp, NeutralMedGam, 'Neutral, Medium Dispersal', 2)
nexttile
myplot(table2array(Neutral(indexNeutralFull, 'WaitTime')),...
    NeutralFullExp, NeutralFullGam, 'Neutral, Full Dispersal', 3)
nexttile
myplot(table2array(Regional(indexRegionalNo, 'WaitTime')),...
    RegionalNoExp, RegionalNoGam, 'Regional, No Dispersal', 4)
nexttile
myplot(table2array(Regional(indexRegionalMed, 'WaitTime')),...
    RegionalMedExp, RegionalMedGam, 'Regional, Medium Dispersal',5)
nexttile
myplot(table2array(Regional(indexRegionalFull, 'WaitTime')),...
    RegionalFullExp, RegionalFullGam, 'Regional, Full Dispersal', 6)
nexttile
myplot(table2array(Local(indexLocalNo, 'WaitTime')),...
    LocalNoExp, LocalNoGam, 'Local, No Dispersal', 7)
nexttile
myplot(table2array(Local(indexLocalMed, 'WaitTime')),...
    LocalMedExp, LocalMedGam, 'Local, Medium Dispersal', 8)
nexttile
myplot(table2array(Local(indexLocalFull, 'WaitTime')),...
    LocalFullExp, LocalFullGam, 'Local, Full Dispersal', 9)

function retval = myplot(target, targetExp, targetGam, targetTitle, index)
%target = table2array(Regional(indexRegionalMed, 'WaitTime'));
%targetExp = RegionalMedExp;
%targetGam = RegionalMedGam;
%targetTitle = 'Regional, Medium Dispersal';

histogram(target, 'Normalization', 'cdf', 'DisplayStyle', 'stairs')
hold on
x = linspace(0, max(target));
plot(x, expcdf(x, targetExp(1)))
plot(x, gamcdf(x, targetGam(1), targetGam(2)))
if(index == 4)
    ylabel('CDF')
end
if(index == 6) 
    legend({'Observed', 'Exponential', 'Gamma'}, 'Location', 'eastoutside') 
    legend('boxoff')
end
if(index == 8)
    xlabel('Inter-event Time Interval Length')
end
title(targetTitle)
hold off

retval = 0; return;
end
