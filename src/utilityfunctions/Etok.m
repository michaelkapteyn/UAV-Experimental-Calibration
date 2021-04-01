function [k] = Etok(E)
k0 = 0.677; %prior slope from the model
translation_factor = 0.8496; % coefficient of proportionality between k/k0 to E

k = k0*(translation_factor*E + (1-translation_factor));
end

