%% Calibration Step 2: Calibrate young's modulus scaling, e
% This script performs the particle-filter update of p(e)
% 
% Input: Data read-in from datafiles/step2
% 
% Output: Recreates output shown in Figure 5, as well as Table 1 column 4, 
% in the paper
%% Read in data:
%%
clear all;
addpath('utilityfunctions')
model_data = readtable('datafiles/step2/load_displacement_model.csv')
measured_data = readtable('datafiles/step2/load_displacement_measurements.csv')
nMeasurements = length(measured_data.displacement);
%% Set up priors on e and k, as well as likelihood distribution
%%
%% Prior on e
mu_e = 1;
Q_e = (0.05/1.96)^2; % 95%CI = +/- 5 percent

%% Prior on k
mu_k = (model_data.load(4)-model_data.load(1))/(model_data.displacement(4)-model_data.displacement(1));
Q_k = ((Etok(1.05)-Etok(1.00))/1.96)^2; % convert Q_e into Q_k via linear relationship defined in Etok

% set up Distributions on x and f
mu_x = measured_data.displacement;
R_x = (0.5/1.96)^2; % 95%CI = +/- 0.5mm

mu_f = measured_data.load;
R_f = ((10e-3)*9.81/1.96)^2; % 95%CI = +/- 10g

% Plot measurements and 95%CI error bars in f=kx space
figure;
hold on;
for i = 1:nMeasurements
    if i <= 4
        errorbar(mu_x(i),mu_f(i),1.96*sqrt(R_f),1.96*sqrt(R_f),1.96*sqrt(R_x),1.96*sqrt(R_x),'r', 'LineWidth',2,'DisplayName',sprintf('Measurement %i',i));
    else
        errorbar(mu_x(i),mu_f(i),1.96*sqrt(R_f),1.96*sqrt(R_f),1.96*sqrt(R_x),1.96*sqrt(R_x),'b','LineWidth',2,'DisplayName',sprintf('Measurement %i',i));
    end

end
xlabel('Tip displacement, x (mm)')
ylabel('Applied tip force, f (N)')
ax = gca;
set(ax, 'FontSize', 16);
xlim([0,16])
ylim([0,10])
title('Measured load-displacement data')
%% Now perform Bayesian update via particle filter algorithm:
%%

%% Draw samples from f and x
nSamples = 1e6; % number of particles
for i = 1:nMeasurements
    sampled_x(:,i) = normrnd(mu_x(i),sqrt(R_x), nSamples,1);
    sampled_f(:,i) = normrnd(mu_f(i),sqrt(R_f), nSamples,1);
%     scatter(sampled_x(:,i),sampled_f(:,i));  %if you want to add samples to previous plot
end

%% Convert each sampled pair {f,x} into an equivalent sample in terms of {e}
for i= 1:nMeasurements
    sampled_k(:,i) = sampled_f(:,i) ./ sampled_x(:,i);
    sampled_e(:,i) = ktoE(sampled_k(:,i));
end

%% set up likelihood functions by fitting a kernel density to p(e) (it is non-Gaussian)
for i= 1:8
    p(i) = fitdist(sampled_e(:,i),'kernel');
end
%% Incorporate each measurement via Bayesian particle filter update
% *Note:* This code block can take some time to run as it uses 1 million samples. 
% 
% Suggest either reducing the number of particles (nParticles) for approximate 
% results, or load the provided .mat file which contains pre-evaluated output 
% variables for 1 million samples
%%
load('datafiles/step2/particlefilteroutput.mat')
% nParticles = 1e6;
% 
% % draw particles from the prior p(e)
% sampled_e = normrnd(mu_e, sqrt(Q_e),nParticles,1);
% % set initial particle weights to be uniform
% w = ones(nParticles,1) ./nParticles;
% w = log(w); % work in logspace for better numerics
% for i = 1:8
%     w = w + log(p(i).pdf(sampled_e));
% %     w = w + log(normpdf(sampled_e, ktoE(mu_f(i)/mu_x(i)), sqrt(R_f)/(mu_x(i)*0.5752) ));
% end
% w = exp(w);
% w = w./sum(w);
% 
% % this function draws uniformly weighted samples from a set of weighted-samples, then fits a kernel density (only used for plotting)
% [posteriorpdf,posteriorsamples] = pdffromweightedsamples(sampled_e,w,1e5);
% save('datafiles/step2/posteriorpdf_new.mat','posteriorpdf')

% % draw 100 samples from posterior distribution of e for use in Step 3:
% Esamples = posteriorpdf.random(100,1)
% csvwrite('datafiles/step3/e_samples.csv',Esamples)
%%
% compute mean and standard deviation of particle filter output (weighted-sample representation)
mean_pf_e = sum(sampled_e.*w) / sum(w);
std_pf_e = sqrt(sum((sampled_e-mean_pf_e).^2.*w) / sum(w));

disp(mean_pf_e)
disp(std_pf_e)
%% Plot prior, likelihoods, and posterior
%%
figure;
hold on;
eplot = linspace(0.7,1.3,1000);
% plot Gaussian prior density
plot(eplot,normpdf(eplot,mu_e,sqrt(Q_e)),'b','linewidth',2)

% plot each likelihood function
for i = 1:8
    plot(eplot,p(i).pdf(eplot), 'k');
end
% plot posterior density
plot(eplot,posteriorpdf.pdf(eplot),'k') %fitted kernel density (used for simplicty in Fig. 5)
histogram(posteriorsamples,15,'Normalization','pdf') % histogram (used in Table 1)
stem(mean_pf_e,40,'b','Linewidth',2,'marker','none')

xlabel("Young's modulus scale factor, e")
ylabel('p(e)')
ax = gca;
set(ax, 'FontSize', 16);
xlim([0.8, 1.2])
ylim([0,40])
title('Prior/posterior densities and likelihood functions for e')
%% Plot prior and posterior distributions in f=kx space
%%
priorkCI = [mu_k-1.96*sqrt(Q_k),mu_k+1.96*sqrt(Q_k)];
posteriorkCI = [Etok(mean_pf_e-1.96*std_pf_e),Etok(mean_pf_e+1.96*std_pf_e)];

% Plot in f=kx space
figure;
hold on;
plot(model_data.displacement,model_data.load,'Linewidth',2,'DisplayName','Prior Mean');
plot(model_data.displacement,Etok(mean_pf_e)*model_data.displacement,'k','Linewidth',2,'DisplayName','Posterior Mean');

% plot prior and posterior confidence bands
fill([0,model_data.displacement(6), model_data.displacement(6)],[0,priorkCI(1)*model_data.displacement(6),priorkCI(2)*model_data.displacement(6)],1,'facecolor', [0, 0.4470, 0.7410], 'edgecolor', 'none', 'facealpha', 0.1,'DisplayName','PriorCI');
fill([0,model_data.displacement(6), model_data.displacement(6)],[0,posteriorkCI(1)*model_data.displacement(6),posteriorkCI(2)*model_data.displacement(6)],1,'facecolor', 'black', 'edgecolor', 'none', 'facealpha', 0.2,'DisplayName','PosteriorCI');

title('prior and posterior distributions in f=kx space')
xlabel('Tip displacement, x (mm)')
ylabel('Applied tip force, f (N)')
ax = gca;
set(ax, 'FontSize', 16);
xlim([0,16])
ylim([0,10])