function [msdLC,msdBNM,msdLCBNM,msdBase,nrgLC,nrgBNM,nrgLCBNM,nrgBase,LClocs,BNMlocs,LCBNMlocs,Baselocs] = MSD_calculations_peakforward(subject_id,figs)

%% MSD Calculation for Paper: 'The structural and functional connectivity between distinct neuromodulatory systems underpins dynamic shifts in brain network topology and attractor landscape topography'
%figs = 1 = yes figures displayed, 0 = no figures displayed
% calculates based on incorporating acceleration and starting points
% saved as 'peakforward.mat'
% Different standard deviations used, but stdThresh = 2 was most compelling
% results

%Example Automation
%     for ii = 1:length(subjects)
%     subject_id = mat2str(subjects(ii,1));   
%     [msdLC,msdBNM,msdLCBNM,msdBase,nrgLC,nrgBNM,nrgLCBNM,nrgBase,LClocs,BNMlocs,LCBNMlocs,Baselocs] = MSD_calculations_peakforward(subject_id,0)
%       end

%loading preprocessed and denoised timeseries 
filename = sprintf('%s%s',subject_id,'_ts_all_clean.mat');
load(filename); %matrix of time-series
ts_cort = ts_all_clean(:,1:333);
cortSig = ts_cort;
lc_ts = ts_all_clean(:,334); %LC time-series
bnm_ts = ts_all_clean(:,335); %nbM time-series

% How many TR post phasic burst are we interested in
stdThresh = 2; %standard deviation threshold was tested with different levels
TRpost = 20;
% Original code by calculating based on incorporating acceleration and starting points
 [msdLC,msdBNM,msdLCBNM,msdBase,nrgLC,nrgBNM,nrgLCBNM,nrgBase,LClocs,BNMlocs,LCBNMlocs,Baselocs] = loopPeakforward(lc_ts,bnm_ts,cortSig,figs,TRpost,stdThresh);
% Save function outputs
savefilename = sprintf('%s%s',subject_id,'_peakforward.mat');
save([savefilename])

end

