%% Final Brainstem Analysis - incorporate correlations with new streamline weights

% copying BFLC weights from server
for ii=1:50 %[1, 10,11,14,15, 16,19,34,43,45,47,49]
    subject = subjects(ii,1);
    path = sprintf('%s%d','X:\PRJ-Brainstem_Tract\BF_fork_v2\',subjects(ii,1));
[filename,destFolder] = copy_files('_Weights_Extracted_BFLC.csv',subject,path,'C:\Users\natas\Documents\PhD\Brainstem_DTI\weight_tracks') 
end

%% Track Weights Analysis
% count number of tracks
for ii=1:length(subjects)
filename = sprintf('%d%s',subjects(ii,1),'_Weights_Extracted_BFLC.csv');
a = dlmread(filename);
no_tracks(ii,:) = length(a); %loads number of tracks between BNM and LC
end
save('number_tcks_new_BNMLC.mat','no_tracks')
% sum weights of tracks
sum_weights_tracks = zeros(50,1);
for ii=1:length(subjects)
filename = sprintf('%d%s',subjects(ii,1),'_Weights_Extracted_BFLC.csv');
a = dlmread(filename);
sum_weights_tracks(ii,:) = sum(a); %loads number of tracks between BNM and LC
end
save('sum_weights_tcks_new_BNMLC.mat','sum_weights_tracks')

%% remove subject - did not have time-series for all ROIs (missing brainstem)
load('subjects_all.mat')
subjects(34,:) = []; %remove

sum_weights_tracks(34,:) = [];

%% Time-series of Cortex + Brainstem Calculated together
%MTD extractions (time-series of cort and brainstem ROI concatenated and
%zscored)
for ii=1:length(subjects)
    filename = sprintf('%d%s',subjects(ii,1),'_ts_all.mat');
    load(filename); %load time-series, been bandpass filtered and zscored per node
    mtd = coupling(ts,20); %window =20
    nROI = 343; %333 + 10 ROIs
    for nn = 1:nROI
    template = find(triu(ones(nROI))-eye(nROI));%finding unique combination of pairs
    end
    time = 900;
    for tt = 1:time
    temp = mtd(:,:,tt);
    mtd_flat(:,tt) = temp(template); %flattens mtd across timepoints
    sprintf('%d%s',tt,'completed mtd flat');
    end
    savefilename_mtd = sprintf('%d%s',subjects(ii,1),'_mtd_cort_brainstem.mat');
    savefilename_mtdflat = sprintf('%d%s',subjects(ii,1),'_mtd_flat_cort_brainstem.mat');
    save([savefilename_mtd],'mtd'); %save raw mtd file based on subnum
    save([savefilename_mtdflat],'mtd_flat'); %save raw mtd file based on subnum
end
% all mtd into large matrix
% take average over time for each participant
mtd_all = zeros(length(subjects),343,343);
for ii=1:length(subjects)
    filename = sprintf('%d%s',subjects(ii,1),'_mtd_cort_brainstem.mat');
    load(filename);
    a = mean(mtd,[3]);
    mtd_all(ii,:,:) = a;
end

%loading all mtd_flat files into one
mtd_flat_all = zeros(length(subjects),58653,900);
for ii=1:length(subjects)
    filename = sprintf('%d%s',subjects(ii,1),'_mtd_flat_cort_brainstem.mat');
    load(filename);
    mtd_flat_all(ii,:,:) =mtd_flat;
end
% Integration calculation
for ii=1:length(subjects)
subnum = sprintf('%d',subjects(ii,1));
dataname = sprintf('%d%s',subjects(ii,1),'_mtd_cort_brainstem.mat'); 
load(dataname);
[ci,q,p,z,hc,f] = integration_plus5(mtd,1,0.75);

save([subnum '_community.mat'],'ci');
save([subnum '_partcoeff.mat'],'p');
save([subnum '_zmodule.mat'],'z');
save([subnum '_modularity.mat'],'q');
save([subnum '_cartographic.mat'],'hc');
save([subnum '_flexibility.mat'],'f');
end

%% Participation Coefficient
% take average over time for each participant
pc_mean_all = zeros(length(subjects),338);
pc_std_all = zeros(length(subjects),338);
pc_all = zeros(length(subjects),338,1200);
for ii=1:length(subjects)
    filename = sprintf('%d%s',subjects(ii,1),'_partcoeff.mat');
    load(filename);
    a = mean(p,2); % 343x900, average across time
    std = std(p,[],2);
    pc_all(ii,:,:) = p; %no cutting for this one (900 time points)
    pc_mean_all(ii,:) = a;
    pc_std_all(ii,:)=std;
    clear a
    clear std
end
save('pc_all_rest.mat','pc_all')
save('pc_mean_all_rest.mat','pc_mean_all')
save('pc_std_all_rest.mat','pc_std_all')

[orig_corr,pval,sig] = permutation_correlate_nonmatch(pc_mean_all,sum_weights_tracks,5000,'Kendall');
sum(sig)
sig_cor_pc_mean_roi_tcks = orig_corr.*sig;
save('sig_corr_pc_mean_roi_tcks.mat','sig_cor_pc_mean_roi_tcks')
[orig_corr,pval,sig] = permutation_correlate_nonmatch(pc_std_all,sum_weights_tracks,5000,'Kendall');
sum(sig)
sig_cor_pc_std_roi_tcks = orig_corr.*sig;
save('sig_corr_pc_std_roi_tcks.mat','sig_cor_pc_mean_roi_tcks')
% See Munn et al. (2021) code for cross-correlation analysis of
% participation coefficient post phasic burst in LC and nbM
%% Figures

figure
scatter(sum_weights_tracks,mean(pc_all,[2]))
xlabel('sum weight tracks between LCBNM')
ylabel('left and right LC pc for each subject')
Fit = polyfit(sum_weights_tracks,mean(pc_all,[2]),1);
f = polyval(Fit,sum_weights_tracks);
hold on
plot(sum_weights_tracks,f,'r')









%% Analysis of no significant difference
% z module
for ii=1:length(subjects)
    filename = sprintf('%d%s',subjects(ii,1),'_zmodule.mat');
    load(filename);
    z_cut = z;
    z_cut(:,1200) = [];
    z_avg = mean(z_cut,[2]);% 333x1200, average across time
    z_std_indiv = std(z_cut,[],2);
    z_all(ii,:) = z_avg;
    z_std(ii,:) = z_std_indiv;% standard deviation of z module
end

% modularity
for ii=1:length(subjects)
    filename = sprintf('%d%s',subjects(ii,1),'_modularity.mat');
    load(filename);
    modularity_all(ii,:) = q;
end
[r_modularity,p_modularity] = corr(modularity_all,no_tracks,'type','Spearman');

% Cartographic Profile
hc_all = zeros(49,101,101,900);
hc_mean_all = zeros(49,10201);
hc_std_all = zeros(49,10201);
for ii=1:length(subjects)
    filename = sprintf('%d%s',subjects(ii,1),'_cartographic.mat');
    load(filename);
    %hc_flat = reshape(hc,10201,900);
    hc_avg = mean(hc,3); %mean across time
    hc_flat_avg =reshape(hc_avg,10201,1);% 101x101x900, average across time
    hc_std = std(hc,[],3); %calc standard deviation
    hc_flat_std =reshape(hc_std,10201,1);%  flatten
    hc_mean_all(ii,:) = hc_flat_avg;
    hc_std_all(ii,:) = hc_flat_std;
    hc_all(ii,:,:,:) = hc; %all subs raw cartograph
end
save('cartographic_all.mat','hc_all')
save('cartographic_mean_all.mat','hc_mean_all')
save('cartographic_std_all.mat','hc_std_all')
[orig_corr_hc,pval_hc,sig_hc] = permutation_correlate_nonmatch(hc_mean_all,sum_weights_tracks,5000,'Pearson');
sum(sig_mean_pearson)
sig_cor_hc_tcks = orig_corr_mean_pearson.*sig_mean_pearson;
save('cartograph_vs_new_tcks.mat','sig_cor_hc_tcks')