function [exp_ic, exp_rate, omega, zeta, components] = fit_damped_cosine(time,measured_signal,frequencies,optionPlot)

% simple exponential decay function
dampedcosine = @(a,x,fn) a(1).*exp(-a(2).*x).*cos(2*pi*fn*(x-x(1)-a(3)));

parameter_guess = [0.0061,1.0,0.0,...
                   0.0137,1.0,0.0];
               


% output
composite_signal = @(a,x) dampedcosine(a(1:3),x, frequencies(1)) + dampedcosine(a(4:6),x, frequencies(2));

% fit coefficients to data nlinfit function
coeff = nlinfit(time,measured_signal,composite_signal,parameter_guess);

% return fitted coefficients:
exp_ic = [coeff(1), coeff(4)];
exp_rate = [coeff(2), coeff(5)];
omega = 2*pi*[frequencies(1),frequencies(2)];
zeta = [abs(coeff(2))./(2*pi*frequencies(1)), abs(coeff(5))./(2*pi*frequencies(2))];

% optionally: plot the fitted function
if optionPlot== 1, plot(time,composite_signal(coeff,time),'r'); end

components(:,1) = dampedcosine(coeff(1:3),time, frequencies(1));
components(:,2) = dampedcosine(coeff(4:6),time, frequencies(2));
end