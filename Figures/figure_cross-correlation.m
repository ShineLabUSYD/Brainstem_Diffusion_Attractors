figure
subplot(1,3,1)
plot(-10:1:10,mean(xcf_lc_corr1(11:31,:,:),3)-mean(xcf_nbm_corr1(11:31,:,:),3),'Color',[.2 .3 1])
hold on
plot(-10:1:10,mean(mean(xcf_lc_corr1(11:31,:,:),3),2)-mean(mean(xcf_nbm_corr1(11:31,:,:),3),2),'LineWidth',10,'Color',[0 0 1])
xlabel('TR')
ylabel('corr #SL with r(LC-nBM,Part)')
set(gcf,'color','w');
plot([-10 10],[-.097 -.097],'-k','LineWidth',1)
plot([-10 10],[.12 .12],'-k','LineWidth',1)
plot([0 0],[-.7 .7],'--k')
ylim([-.7 .7])
title('LC - nBM')


subplot(1,3,2)
plot(-10:1:10,mean(xcf_lc_corr1(11:31,:,:),3),'Color',[1 .4 0])
hold on
plot(-10:1:10,mean(mean(xcf_lc_corr1(11:31,:,:),3),2),'LineWidth',10,'Color',[1 0 0])
xlabel('TR')
ylabel('corr #SL with r(LC,Part)')
set(gcf,'color','w');
plot([-10 10],[-.097 -.097],'-k','LineWidth',1)
plot([-10 10],[.12 .12],'-k','LineWidth',1)
ylim([-.7 .7])
plot([0 0],[-.7 .7],'--k')
title('LC')

subplot(1,3,3)
plot(-10:1:10,mean(xcf_nbm_corr1(11:31,:,:),3),'Color',[.1 .9 .5])
hold on
plot(-10:1:10,mean(mean(xcf_nbm_corr1(11:31,:,:),3),2),'LineWidth',10,'Color',[0 1 .5])
xlabel('TR')
ylabel('corr #SL with r(nBM,Part)')
set(gcf,'color','w');
plot([-10 10],[-.097 -.097],'-k','LineWidth',1)
plot([-10 10],[.12 .12],'-k','LineWidth',1)
ylim([-.7 .7])
plot([0 0],[-.7 .7],'--k')
title('nBM')

