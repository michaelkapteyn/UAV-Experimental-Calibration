function [E] = ktoE(k)
k0 = 0.677; %prior slope from the model
translation_factor = 0.8496; % coefficient of proportionality between k/k0 to e

E = ((k/k0 -(1-translation_factor)) /translation_factor);
end

