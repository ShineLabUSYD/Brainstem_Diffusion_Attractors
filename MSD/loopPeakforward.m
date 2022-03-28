function [msdLC,msdBNM,msdLCBNM,msdBase,nrgLC,nrgBNM,nrgLCBNM,nrgBase,LClocs,BNMlocs,LCBNMlocs,Baselocs] = loopPeakforward(lc_ts,bnm_ts,cortSigZscore,figs,TRpost,stdThresh)



%% MSD calculation and Energy plots
%%
% 
%cortSigZscore = zscore(cortSig');
%cortSigZscore = cortSigZscore';
%%
%number of tr into future calc msd
ndt = 20;


maxMSd = 10;
ds = 0:1:maxMSd; % the msd range calculated across


nrgLC = nan(ndt-1,numel(ds));
nrgBNM = nan(ndt-1,numel(ds));
nrgLCBNM = nan(ndt-1,numel(ds));
nrgBase = nan(ndt-1,numel(ds));

%TRpost = 20;
%Get the locations of just peaks in diff signal
[LClocs,BNMlocs,LCBNMlocs,Baselocs] = LcBnmPkTimeTash_individualSubjectpeakforward(lc_ts,bnm_ts,TRpost,stdThresh);

for dt = 2:ndt

%%
% MSD calculation
MSD = mean( (cortSigZscore(1+dt:end,:) - cortSigZscore(1:end-dt,:)).^2,2) ;

msdLC = MSD(LClocs);
msdBNM = MSD(BNMlocs);
msdLCBNM = MSD(LCBNMlocs);
msdBase = MSD(Baselocs);

%% Calculate probability distribution  and energy for each dt and each neuromod

[nrgLCdt,nrgBNMdt,nrgBasedt,nrgLCBNMdt] = lcBnmPdistn(msdLC,msdBNM,msdLCBNM,msdBase,ds);

% Pool results across time
nrgLC(dt-1,:) = nrgLCdt;
nrgBNM(dt-1,:) = nrgBNMdt;
nrgLCBNM(dt-1,:) = nrgLCBNMdt;
nrgBase(dt-1,:) = nrgBasedt;

end


%%
%  Plot estimated MSD energy landscape 
% Fig. 2C
if figs ==1 
    x = 1:ndt-1;
    y = ds;
    [X,Y] = meshgrid(x,y);

    xmax = ndt;

    figure
    subplot(2,2,1)
    mesh(X,Y,nrgBase','EdgeColor', [105,105,105]./255)
    xlabel('TR')
    ylabel('MSD')
    zlabel('MSD  energy')
    xlim([1 xmax])
    %zlim([2 78])
    ylim([1 maxMSd])
    view(-15,30)   % XZ
    title('Baseline')

    subplot(2,2,2)
    mesh(X,Y,nrgLC'-nrgBase','EdgeColor', [236 102 102]./255)
    xlim([1 xmax])
    %zlim([2 78])
    ylim([1 maxMSd])
    view(-15,30)   % XZ
    xlabel('TR')
    ylabel('MSD')
    zlabel('MSD  energy')
    title('LC')


    subplot(2,2,3)
    mesh(X,Y,nrgBNM'-nrgBase','EdgeColor', [60 184 79]./255)
    xlim([1 xmax])
    %zlim([2 78])
    ylim([1 maxMSd])
    view(-15,30)   % XZ
    xlabel('TR')
    ylabel('MSD')
    zlabel('MSD  energy')
    title('BNM')

    subplot(2,2,4)
    mesh(X,Y,nrgLCBNM'-nrgBase','EdgeColor', 'b')
    xlim([1 xmax])
    %zlim([2 78])
    ylim([1 maxMSd])
    view(-15,30)   % XZ
    xlabel('TR')
    ylabel('MSD')
    zlabel('MSD  energy')
    title('LC+BNM')
else 
    sprintf('no figures chosen')
end

end
