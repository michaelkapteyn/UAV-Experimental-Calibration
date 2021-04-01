function [mu_posterior,stddev_posterior] = propagateNormalStatistics(mu_prior,stddev_prior,mu_likelihood,stddev_likelihood)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu_posterior = (mu_prior*stddev_likelihood.^2 + mu_likelihood*stddev_prior^2)/(stddev_prior^2+stddev_likelihood^2);
stddev_posterior = sqrt((stddev_likelihood^2*stddev_prior^2) / (stddev_likelihood^2+stddev_prior^2));
end

