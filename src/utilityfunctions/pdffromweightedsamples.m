function [pdfobject, unweightedsamples] = pdffromweightedsamples(samples,weights,nSamples)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%normalize weights
weights = weights./sum(weights);


%sort samples and weights
[sortedsamples,I] = sort(samples);
sortedweights = weights(I);

%get empirical cdf
empcdf = cumsum(sortedweights);

Usamples = rand(nSamples,1);
unweightedsamples = zeros(length(Usamples),1);
for i = 1:length(Usamples)
    u = Usamples(i);
    idx = find((empcdf > u),1,'first');
    unweightedsamples(i) = sortedsamples(idx);
end

pdfobject = fitdist(unweightedsamples, 'kernel');

