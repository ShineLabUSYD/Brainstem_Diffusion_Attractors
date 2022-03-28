function [LClocs,BNMlocs,LCBNMlocs,Baselocs] = LcBnmPkTimeTash_individualSubjectpeakforward(lc_ts,bnm_ts,TRpost,stdThresh)
% Original code by calculating based on incorporating acceleration and
% starting points

%% Free paramters try playing aroudn with
% Original calculation with acceleration

%% For example process we have concatenated all recordings

%signal length
lengthSignal = 900;

% How many TR post phasic burst are we interested in
%TRpost = 20;

% That means we do not want any locations 20tr before the end of a
% recording
edgvals = lengthSignal-30-1:lengthSignal;


% When we compare the phasic burst to (BSI paper we set delay = 10)
peakDelay = 5; %doesn't matter


%std threshold
%stdThresh = 1.5;


%% First LC-BNM peaks 


sig = lc_ts-bnm_ts; %grab signal

[~,locs] =  findpeaks(sig,'MinPeakDistance',1,'MinPeakheight',stdThresh*std(sig)); %find peaks of acc

%Remove locs near boundaries
locs = locs(~ismember(locs,edgvals));

LClocs = locs;

%% BNM-LC peaks

sig = bnm_ts-lc_ts; %grab signal

[~,locs] =  findpeaks(sig,'MinPeakDistance',1,'MinPeakheight',stdThresh*std(sig)); %find peaks of acc

%Remove locs near boundaries
locs = locs(~ismember(locs,edgvals));
BNMlocs = locs;

%% LC+BNM peaks

sig = bnm_ts+lc_ts; %grab signal

[~,locs] =  findpeaks(sig,'MinPeakDistance',1,'MinPeakheight',stdThresh*std(sig)); %find peaks of acc

%Remove locs near boundaries
locs = locs(~ismember(locs,edgvals));

LCBNMlocs = locs;

%% Baseline times

%compile locations of phasic bursts
allLocs=[LClocs; BNMlocs; LCBNMlocs];

%trick to avoid -TR before a acc peak
nearLocs = repmat(-TRpost:0,numel(allLocs),1)+allLocs;

%then get the indices that arent edge values or near a loc
allindx = 1:numel(sig);
Baselocs = allindx(~ismember(allindx,[nearLocs(:); edgvals(:)]));




end