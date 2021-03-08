figure(1)
hold on
box on

xlim([0,1200]);
ylim([0,90]);
%YF=load('ADScase1_NB5RPM1400_SPLBMm.txt');
YF=load('HansonADScase1_SPLHansonm.txt');
LF=load('fwh.Mic_131_spl_13_NB5RPM1400.txt');
HF=load('fwh.Mic_131_spl_13_HFNB5RPM1400.txt');

stem(YF(:,1),YF(:,2),'b','filled','LineWidth',2)
stem(YF(:,1),YF(:,3),'o','filled','LineWidth',2)
stem(YF(:,1),YF(:,4),'r','filled','LineWidth',2)
plot(LF(:,1),LF(:,2),'g','LineWidth',2)
plot(HF(:,1),HF(:,2),'m','LineWidth',2)


xlabel('Frequency, Hz','interpreter','latex')
ylabel('SPL, dB','interpreter','latex')


set(gcf, 'PaperPositionMode','Auto')   
print -deps 'SPL_ADScase1.eps'


legend('Thickns','Loading','Total','BEMT+Pnoise','High-fidelity','latex')
legend boxoff

%%

figure(2)
hold on
box on

xlim([0,1200]);
ylim([0,100]);
%YF=load('ADScase3_NB5RPM1900_SPLBMm.txt');
YF=load('HansonADScase3_SPLHansonm.txt');
LF=load('fwh.Mic_131_spl_13_NB5RPM1900.txt');
HF=load('fwh.Mic_131_spl_13_HFNB5RPM1900.txt');

stem(YF(:,1),YF(:,2),'b','filled','LineWidth',2)
stem(YF(:,1),YF(:,3),'o','filled','LineWidth',2)
stem(YF(:,1),YF(:,4),'r','filled','LineWidth',2)
plot(LF(:,1),LF(:,2),'g','LineWidth',2)
plot(HF(:,1),HF(:,2),'m','LineWidth',2)


xlabel('Frequency, Hz','interpreter','latex')
ylabel('SPL, dB','interpreter','latex')

set(gcf, 'PaperPositionMode','Auto')   
print -deps 'SPL_ADScase3.eps'


legend('Thickns','Loading','Total','BEMT+Pnoise','High-fidelity','latex')
legend boxoff

%%

figure(3)
hold on
box on

xlim([0,1200]);
ylim([0,100]);
%YF=load('ADScase8_NB7RPM1900_SPLBMm.txt');
YF=load('SPL/ADScase8_SPLH_Mic1.txt');
LF=load('fwh.Mic_131_spl_13_NB7RPM1900.txt');
HF=load('fwh.Mic_131_spl_13_HFNB7RPM1900.txt');

stem(YF(:,1),YF(:,2),'b','filled','LineWidth',2)
stem(YF(:,1),YF(:,3),'o','filled','LineWidth',2)
stem(YF(:,1),YF(:,4),'r','filled','LineWidth',2)
plot(LF(:,1),LF(:,2),'g','LineWidth',2)
plot(HF(:,1),HF(:,2),'m','LineWidth',2)


xlabel('Frequency, Hz','interpreter','latex')
ylabel('SPL, dB','interpreter','latex')

set(gcf, 'PaperPositionMode','Auto')   
print -deps 'SPL_ADScase8.eps'


legend('Thickns','Loading','Total','BEMT+Pnoise','High-fidelity','latex')
legend boxoff
