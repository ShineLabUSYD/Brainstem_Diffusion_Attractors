function [nrg] = PdistGaussKern(dat,ds)

%dat is msd values
pd = fitdist(dat,'Kernel','BandWidth',4);
yLC = pdf(pd,ds);
nrg = -1.*log(yLC);

end