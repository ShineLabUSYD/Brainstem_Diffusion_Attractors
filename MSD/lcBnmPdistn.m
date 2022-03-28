function [nrgLCdt,nrgBNMdt,nrgBasedt,nrgLCBNMdt] = lcBnmPdistn(msdLC,msdBNM,msdLCBNM,msdBase,ds)
%ds spatial time-step, derivative, MSD steps
dat = msdLC;
pd = fitdist(dat,'Kernel','BandWidth',4);
yLC = pdf(pd,ds); %calculate proability distribution, at given MSD
nrgLCdt = -1.*log(yLC); %calculate energy, inverse of probability


dat = msdBNM;
pd = fitdist(dat,'Kernel','BandWidth',4);
yBNM = pdf(pd,ds);
nrgBNMdt = -1.*log(yBNM);


dat = msdBase;
pd = fitdist(dat,'Kernel','BandWidth',4);
yBase = pdf(pd,ds);
nrgBasedt = -1.*log(yBase);

dat = msdLCBNM;
pd = fitdist(dat,'Kernel','BandWidth',4);
yLCBNM = pdf(pd,ds);
nrgLCBNMdt = -1.*log(yLCBNM);

end