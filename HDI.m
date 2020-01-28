function HDIofMCMC = HDI(sampleVec,credMass)
% Computes highest density interval from a sample of representative values,
% estimated as shortest credible interval.
% Arguments:
% sampleVec is a vector of representative values from a probability distribution
% credMass is a scalar between 0 and 1, indicating the mass within the credible
% interval that is to be estimated.
% Value:
% HDIlim is a vector containing the limits of the HDI 
sortedPts = sort(sampleVec);
ciIdxInc = ceil(credMass*length(sortedPts));
nCIs = length(sortedPts) - ciIdxInc; %number of intervals
ciWidth = zeros(1, nCIs); %width of interval
for i = 1:nCIs
    ciWidth(i) = sortedPts(i + ciIdxInc) - sortedPts(i);
end
[~,id] = min(ciWidth);
HDImin = sortedPts(id);
HDImax = sortedPts(id + ciIdxInc);
HDIlim = [HDImin, HDImax];
HDIofMCMC = HDIlim;
end