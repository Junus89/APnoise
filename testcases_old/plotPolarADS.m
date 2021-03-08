figure(1)
grid on

YF=load('ADScase1_NB5RPM1400_SPLBMr.txt');
LF=load('Directivity-xzs_NB5RPM1400.txt');
HF=load('Directivity-xzs_HfNB5RPM1400.txt');

DF=1.2*max(HF(:,2))*ones(length(HF),1);

p=polarplot(HF(:,1),DF,'w');
set(p,'visible','off');
hold on
polarplot(HF(:,1),HF(:,2),'ko');
