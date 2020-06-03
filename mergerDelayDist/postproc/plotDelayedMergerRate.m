root = "D:\Dropbox\Projects\20181213_BatseSgrbRedshift\git\mergerDelayDist\build\winx64\intel\19.0.4.245\release\static\serial\B10\romberg\bin";
d = importdata(fullfile(root,"mergerDelayRateB10.txt"));
figure
plot(d.data(:,1),d.data(:,3))
hold on
plot(d.data(:,1),d.data(:,2))