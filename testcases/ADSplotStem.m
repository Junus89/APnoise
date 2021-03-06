figure(1)
hold on
box on

xlim([0,1200]);
ylim([0,90]);
%YF=load('ADScase1_NB5RPM1400_SPLBMm.txt');
YF=load('SPL/ADScase1_SPLH_Mic0.txt');
LF=load('fwh.Mic_131_spl_13_NB5RPM1400.txt');
LFT=load('fwh.Mic_T_131_spl_13_NB5RPM1400.txt');
LFL=load('fwh.Mic_L_131_spl_13_NB5RPM1400.txt');
HF=load('fwh.Mic_131_spl_13_HFNB5RPM1400.txt');

stem(YF(:,1),YF(:,2),'b','filled','LineWidth',2)
stem(YF(:,1),YF(:,3),'o','filled','LineWidth',2)
stem(YF(:,1),YF(:,4),'r','filled','LineWidth',2)
plot(LFT(:,1),LFT(:,2),'k-','LineWidth',2)
plot(LFL(:,1),LFL(:,2),'p-.','LineWidth',2)
plot(LF(:,1),LF(:,2),'g','LineWidth',2)
plot(HF(:,1),HF(:,2),'m','LineWidth',2)


xlabel('Frequency, Hz','interpreter','latex')
ylabel('SPL, dB','interpreter','latex')


set(gcf, 'PaperPositionMode','Auto')   
print -deps 'SPL_ADScase1.eps'


legend('Thickns','Loading','Total','BPT','BPL','BEMT+Pnoise','High-fidelity','latex')
legend boxoff

%%

figure(2)
hold on
box on

xlim([0,1200]);
ylim([0,100]);
%YF=load('ADScase3_NB5RPM1900_SPLBMm.txt');
YF=load('SPL/ADScase3_SPLH_Mic0.txt');
LF=load('fwh.Mic_131_spl_13_NB5RPM1900.txt');
LFT=load('fwh.Mic_T_131_spl_13_NB5RPM1900.txt');
LFL=load('fwh.Mic_L_131_spl_13_NB5RPM1900.txt');
HF=load('fwh.Mic_131_spl_13_HFNB5RPM1900.txt');

stem(YF(:,1),YF(:,2),'b','filled','LineWidth',2)
stem(YF(:,1),YF(:,3),'o','filled','LineWidth',2)
stem(YF(:,1),YF(:,4),'r','filled','LineWidth',2)
plot(LFT(:,1),LFT(:,2),'k-','LineWidth',2)
plot(LFL(:,1),LFL(:,2),'p-.','LineWidth',2)
plot(LF(:,1),LF(:,2),'g','LineWidth',2)
plot(HF(:,1),HF(:,2),'m','LineWidth',2)


xlabel('Frequency, Hz','interpreter','latex')
ylabel('SPL, dB','interpreter','latex')

set(gcf, 'PaperPositionMode','Auto')   
print -deps 'SPL_ADScase3.eps'


legend('Thickns','Loading','Total','BPT','BPL','BEMT+Pnoise','High-fidelity','latex')
legend boxoff


%%

figure(3)
hold on
box on

xlim([0,1200]);
ylim([0,100]);
%YF=load('ADScase8_NB7RPM1900_SPLBMm.txt');
YF=load('SPL/ADScase8_SPLH_Mic0.txt');
LFT=load('fwh.Mic_T_131_spl_13_NB7RPM1900.txt');
LFL=load('fwh.Mic_L_131_spl_13_NB7RPM1900.txt');
LF=load('fwh.Mic_131_spl_13_NB7RPM1900.txt');
HF=load('fwh.Mic_131_spl_13_HFNB7RPM1900.txt');

stem(YF(:,1),YF(:,2),'b','filled','LineWidth',2)
stem(YF(:,1),YF(:,3),'o','filled','LineWidth',2)
stem(YF(:,1),YF(:,4),'r','filled','LineWidth',2)
plot(LFT(:,1),LFT(:,2),'k-','LineWidth',2)
plot(LFL(:,1),LFL(:,2),'p-.','LineWidth',2)
plot(LF(:,1),LF(:,2),'g','LineWidth',2)
plot(HF(:,1),HF(:,2),'m','LineWidth',2)


xlabel('Frequency, Hz','interpreter','latex')
ylabel('SPL, dB','interpreter','latex')

set(gcf, 'PaperPositionMode','Auto')   
print -deps 'SPL_ADScase8.eps'


legend('Thickns','Loading','Total','BPT','BPL','BEMT+Pnoise','High-fidelity','latex')
legend boxoff
%% ==== NB7 vs NB5 RPM 1900 at the first BPF

figure(4)
hold on
box on

xlim([0,1200]);
ylim([0,100]);
YFc3=load('SPL/ADScase3_SPLH_Mic0.txt');
YFc8=load('SPL/ADScase8_SPLH_Mic0.txt');

stem(YFc3(:,1),YFc3(:,2),'b','filled','LineWidth',2)
stem(YFc3(:,1),YFc3(:,3),'g','filled','LineWidth',2)
stem(YFc3(:,1),YFc3(:,4),'r','filled','LineWidth',2)

stem(YFc8(:,1),YFc8(:,2),'c-.','filled','LineWidth',2)
stem(YFc8(:,1),YFc8(:,3),'p-.','filled','LineWidth',2)
stem(YFc8(:,1),YFc8(:,4),'m-.','filled','LineWidth',2)



xlabel('Frequency, Hz','interpreter','latex')
ylabel('SPL, dB','interpreter','latex')

set(gcf, 'PaperPositionMode','Auto')   
print -deps 'SPL_ADScase5vs8.eps'


legend('NB 5 T','NB5 L','NB5 Total','NB 7 T','NB 7 L','NB 7 Total','latex')
legend boxoff

%% ==== RPM 1400 vs RPM 1900 at the first BPF

figure(5)
hold on
box on

xlim([0,1200]);
ylim([0,100]);
YFc3=load('SPL/ADScase3_SPLH_Mic0.txt');
YFc1=load('SPL/ADScase1_SPLH_Mic0.txt');

stem(YFc1(:,1),YFc1(:,2),'c-.','filled','LineWidth',2)
stem(YFc1(:,1),YFc1(:,3),'p-.','filled','LineWidth',2)
stem(YFc1(:,1),YFc1(:,4),'m-.','filled','LineWidth',2)

stem(YFc3(:,1),YFc3(:,2),'b','filled','LineWidth',2)
stem(YFc3(:,1),YFc3(:,3),'g','filled','LineWidth',2)
stem(YFc3(:,1),YFc3(:,4),'r','filled','LineWidth',2)


xlabel('Frequency, Hz','interpreter','latex')
ylabel('SPL, dB','interpreter','latex')

set(gcf, 'PaperPositionMode','Auto')   
print -deps 'SPL_ADScase1vs3.eps'


legend('RPM 1400 T','RPM 1400 L','RPM 1400 Total','RPM 1900 T','RPM 1900 L','RPM 1900 Total','latex')
legend boxoff
