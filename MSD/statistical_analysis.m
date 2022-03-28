%% Statistical Analysis of MSD - paper: 'The structural and functional connectivity between distinct neuromodulatory systems underpins dynamic shifts in brain network topology and attractor landscape topography'


%% 1. Load Peak forward MSD calculation analysis
%Predefine variables
all_nrgLC = zeros(49,19,11);
all_nrgBNM = zeros(49,19,11);
all_nrgBase = zeros(49,19,11);
all_nrgLCBNM = zeros(49,19,11);
%load into all subject variables
for ii= 1:length(subjects)
subject_id = mat2str(subjects(ii,1));
filename3 = sprintf('%s%s',subject_id,'_peakforward.mat');
load(filename3);
%MSD corresponding burst in activity
all_nrgLC(ii,:,:) = nrgLC;
all_nrgBNM(ii,:,:) = nrgBNM;
all_nrgLCBNM(ii,:,:) = nrgLCBNM;
all_nrgBase(ii,:,:) = nrgBase;
end
save('peakforward_all_subjects.mat','all_nrgBase','all_nrgBNM','all_nrgLC','all_nrgLCBNM')
clear all
load('peakforward_all_subjects.mat');

%Number of locations identified in peaks in activity
for ii= 1:length(subjects)
subject_id = mat2str(subjects(ii,1));
filename3 = sprintf('%s%s',subject_id,'_peakforward.mat');
load(filename3);

no_locsLC(ii,:) = length(LClocs);
no_locsBNM(ii,:) = length(BNMlocs);
no_locsLCBNM(ii,:) = length(LCBNMlocs);
no_locsBase(ii,:) = length(Baselocs);
end


%% avg calculated across subjects
avg_nrgLC = nanmean(all_nrgLC);
avg_nrgLC = reshape(avg_nrgLC,19,11);
avg_nrgBNM = nanmean(all_nrgBNM);
avg_nrgBNM = reshape(avg_nrgBNM,19,11);
avg_nrgLCBNM = nanmean(all_nrgLCBNM);
avg_nrgLCBNM = reshape(avg_nrgLCBNM,19,11);
avg_nrgBase = nanmean(all_nrgBase);
avg_nrgBase = reshape(avg_nrgBase,19,11);

%% standard deviation calculated for each subjects across TR
std_nrgLC = std(all_nrgLC,0,3); %variability of MSD energy across time
std_nrgBNM = std(all_nrgBNM,0,3);
std_nrgLCBNM = std(all_nrgLCBNM,0,3);
std_nrgBase = std(all_nrgBase,0,3);

%% Difference - taking away Base MSD energy
% reshape 
reshape_all_nrgBase = reshape(all_nrgBase,49,11,19);
reshape_all_nrgBNM = reshape(all_nrgBNM,49,11,19);
reshape_all_nrgLC = reshape(all_nrgLC,49,11,19);
reshape_all_nrgLCBNM = reshape(all_nrgLCBNM,49,11,19);
diff_all_nrgLC = reshape_all_nrgLC - reshape_all_nrgBase;
diff_all_nrgBNM = reshape_all_nrgBNM - reshape_all_nrgBase;
diff_all_nrgLCBNM = reshape_all_nrgLCBNM - reshape_all_nrgBase;
%% Comparison of MSD with streamline weights between LC and nbM
%Pearson correlation of relative MSD post peaks in LC and nbM, compared to
%streamline weights between LC and nbM
    for jj=1:no_TR %no of TRs
        for kk=1:no_EnergyBins
        rval_nrgLCBNM(:,jj,kk) = corr((all_nrgLCBNM(:,jj,kk)-all_nrgBase(:,jj,kk)),sum_weights_tracks);
        rval_nrgLC(:,jj,kk) = corr((all_nrgLC(:,jj,kk)-all_nrgBase(:,jj,kk)),sum_weights_tracks);
        rval_nrgBNM(:,jj,kk) = corr((all_nrgBNM(:,jj,kk)-all_nrgBase(:,jj,kk)),sum_weights_tracks);
        rval_nrgBase(:,jj,kk) = corr((all_nrgBase(:,jj,kk)-all_nrgBase(:,jj,kk)),sum_weights_tracks);
        end
    end
    
 % run non-parametric permutation testing on the above correlation (5000
 % iterations)
 for jj=1:no_TR %no of TRs
     for kk=1:no_EnergyBins
        pval_all_nrgLCBNMminusBase(:,jj,kk) = perm_1d_corr(sum_weights_tracks,(all_nrgLCBNM(:,jj,kk)-all_nrgBase(:,jj,kk)),5000);
        pval_all_nrgLCminusBase(:,jj,kk) = perm_1d_corr(sum_weights_tracks,(all_nrgLC(:,jj,kk)-all_nrgBase(:,jj,kk)),5000);
        pval_all_nrgBNMminusBase(:,jj,kk) = perm_1d_corr(sum_weights_tracks,(all_nrgBNM(:,jj,kk)-all_nrgBase(:,jj,kk)),5000);
        pval_all_nrgBaseminusBase(:,jj,kk) = perm_1d_corr(sum_weights_tracks,(all_nrgBase(:,jj,kk)-all_nrgBase(:,jj,kk)),5000);
      end
 end
 
rval_nrgLCBNM = reshape(rval_nrgLCBNM,19,11); %reshape for plotting
rval_nrgLC = reshape(rval_nrgLC,19,11); 
rval_nrgBNM = reshape(rval_nrgBNM,19,11); 
rval_nrgBase = reshape(rval_nrgBase,19,11); 

%% Gradient, to determine the change across variables
%xx = (avg_nrgLC'-avg_nrgBase');
%[FX,FY] = gradient(xx);

LC_nrg_diff = all_nrgLC - all_nrgBase; 
BNM_nrg_diff = all_nrgBNM - all_nrgBase;
LCBNM_nrg_diff = all_nrgLCBNM - all_nrgBase;
% LC Gradient
for ii=1:49
    xx = LC_nrg_diff(ii,:,:);
    xx = squeeze(xx);
    [FX,FY] = gradient(xx');
    gradient_LC_FX(ii,:,:) = FX;
    gradient_LC_FY(ii,:,:) = FY;
    clear xx
end
%Correlation of Gradient vs streamline weights
    for kk=1:no_EnergyBins
        for jj=1:no_TR %no of TRs
        rval_gradient_LC_FX(:,kk,jj) = corr(gradient_LC_FX(:,kk,jj),sum_weights_tracks);
        rval_gradient_LC_FY(:,kk,jj) = corr(gradient_LC_FY(:,kk,jj),sum_weights_tracks);
        end
    end
% Permutation Testing
 for kk=1:no_EnergyBins
    for jj=1:no_TR %no of TRs
    pval_gradient_LC_FX(:,kk,jj) = perm_1d_corr(sum_weights_tracks, gradient_LC_FX(:,kk,jj),5000);
    pval_gradient_LC_FY(:,kk,jj) = perm_1d_corr(sum_weights_tracks, gradient_LC_FY(:,kk,jj),5000);
    end
end
pval_gradient_LC_FX = squeeze(pval_gradient_LC_FX);
pval_gradient_LC_FY = squeeze(pval_gradient_LC_FY);
rval_gradient_LC_FX = squeeze(rval_gradient_LC_FX);
rval_gradient_LC_FY = squeeze(rval_gradient_LC_FY);

%nbM gradient
for ii=1:49
    xx = BNM_nrg_diff(ii,:,:);
    xx = squeeze(xx);
    [FX,FY] = gradient(xx');
    gradient_BNM_FX(ii,:,:) = FX;
    gradient_BNM_FY(ii,:,:) = FY;
    clear xx
end
%Correlation of Gradient vs streamline weights
    for kk=1:no_EnergyBins
        for jj=1:no_TR %no of TRs
        rval_gradient_BNM_FX(:,kk,jj) = corr(gradient_BNM_FX(:,kk,jj),sum_weights_tracks);
        rval_gradient_BNM_FY(:,kk,jj) = corr(gradient_BNM_FY(:,kk,jj),sum_weights_tracks);
        end
    end
% Permutation Testing
 for kk=1:no_EnergyBins
    for jj=1:no_TR %no of TRs
    pval_gradient_BNM_FX(:,kk,jj) = perm_1d_corr(sum_weights_tracks, gradient_BNM_FX(:,kk,jj),5000);
    pval_gradient_BNM_FY(:,kk,jj) = perm_1d_corr(sum_weights_tracks, gradient_BNM_FY(:,kk,jj),5000);
    end
end
pval_gradient_BNM_FX = squeeze(pval_gradient_BNM_FX);
pval_gradient_BNM_FY = squeeze(pval_gradient_BNM_FY);
rval_gradient_BNM_FX = squeeze(rval_gradient_BNM_FX);
rval_gradient_BNM_FY = squeeze(rval_gradient_BNM_FY);

%LC + nbM gradient
for ii=1:49
    xx = LCBNM_nrg_diff(ii,:,:);
    xx = squeeze(xx);
    [FX,FY] = gradient(xx');
    gradient_LCBNM_FX(ii,:,:) = FX;
    gradient_LCBNM_FY(ii,:,:) = FY;
    clear xx
end
%Correlation of Gradient vs streamline weights
    for kk=1:no_EnergyBins
        for jj=1:no_TR %no of TRs
        rval_gradient_LCBNM_FX(:,kk,jj) = corr(gradient_LCBNM_FX(:,kk,jj),sum_weights_tracks);
        rval_gradient_LCBNM_FY(:,kk,jj) = corr(gradient_LCBNM_FY(:,kk,jj),sum_weights_tracks);
        end
    end
% Permutation Testing
 for kk=1:no_EnergyBins
    for jj=1:no_TR %no of TRs
    pval_gradient_LCBNM_FX(:,kk,jj) = perm_1d_corr(sum_weights_tracks, gradient_LCBNM_FX(:,kk,jj),5000);
    pval_gradient_LCBNM_FY(:,kk,jj) = perm_1d_corr(sum_weights_tracks, gradient_LCBNM_FY(:,kk,jj),5000);
    end
end
pval_gradient_LCBNM_FX = squeeze(pval_gradient_LCBNM_FX);
pval_gradient_LCBNM_FY = squeeze(pval_gradient_LCBNM_FY);
rval_gradient_LCBNM_FX = squeeze(rval_gradient_LCBNM_FX);
rval_gradient_LCBNM_FY = squeeze(rval_gradient_LCBNM_FY); %no significance found


%% Final Figures Plotted 
% Include meshgrid attractor landscape plot, and correlation between MSD
% results and streamline between LC and nbM
set(gcf,'color','w');
figure
subplot(2,3,1)
baseline = 0.03;
mesh(X,Y,(avg_nrgLC'-avg_nrgBase')+baseline,'EdgeColor', [236 102 102]./255)
xlim([1 xmax])
%zlim([-0.05 0.04])
ylim([1 maxMSd])
view(-15,30)
%view(0,0)   % XZ
zticks([-0.02+baseline,-0.01+baseline,0+baseline,0.005+baseline,0.01+baseline])
zticklabels({'-0.02','-0.01','0','0.005','0.01'})
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
title('LC')
hold on
meshgrid(imagesc(rval_nrgLC.'),imagesc(rval_nrgLC.'),'r') %change colour scheme to have red/blue
colorbar

subplot(2,3,2)
baseline = 0.09;
mesh(X,Y,(avg_nrgBNM'-avg_nrgBase')+baseline,'EdgeColor', [60 184 79]./255)
xlim([1 xmax])
%zlim([-0.08 0.08])
ylim([1 maxMSd])
view(-15,30)
%view(0,0)   % XZ
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
title('BNM')
zticks([-0.02+baseline,0+baseline,0.02+baseline,0.04+baseline])
zticklabels({'-0.02','0','0.02','0.04'})
hold on
imagesc(rval_nrgBNM.') %unsig correlations, plotted at -0.05 of x axis
%meshgrid(imagesc(rval_nrgBNM.'),imagesc(rval_nrgBNM.'),-0.08)
colorbar

subplot(2,3,3)
baseline = 0.06;
mesh(X,Y,(avg_nrgLCBNM'-avg_nrgBase')+baseline,'EdgeColor', 'b')
xlim([1 xmax])
%zlim([2 78])
ylim([1 maxMSd])
view(-15,30) 
%view(0,0)% XZ
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
title('LC+BNM')
zticks([-0.01+baseline,0+baseline,0.01+baseline,0.02+baseline,0.03+baseline])
zticklabels({'-0.05','0','0.01','0.02','0.03'})
hold on
imagesc(rval_nrgLCBNM.') %unsig correlations
colorbar
set(gcf,'color','w')
