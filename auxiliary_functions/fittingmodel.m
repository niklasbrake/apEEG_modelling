function [full_model,AP_model] = fittingmodel(fitType)
% [full_model,AP_model] = FITTINGMODEL(fitType) returns matlab function that
%	implement the equation for the EEG spectrum described in the following
% 	paper:
%
%	Brake, N. et al. A neurophysiological basis for aperiodic eeg and the
%		background spectral trend. Nat. Commun. 15, (2024).
%
%	fitType can be either 'eq1', 'eq5', or 'eq6' and corresponds to the
%	equations in the above paper.
%	full_model returns a funciton that includes gaussin bumps to capture spectral peaks.
%	AP_model is the model for the background spectra trend.

	if(nargin == 0 || strcmp(fitType,'eq6'))
		fitAP = @(x,f) eq6(x(1),x(2),x(3),x(4),2*pi*f);
		full_model = @(f,x) fitAP(x(1:4),f) + fitPeaks(x(5:end),f);
		AP_model = @(f,x) fitAP(x(1:4),f);
	elseif(strcmp(fitType,'eq1'))
		fitAP = @(x,f) eq1(x(1),x(2),x(3),x(4),2*pi*f);
		full_model = @(f,x) fitAP(x(1:4),f) + fitPeaks(x(5:end),f);
		AP_model = @(f,x) fitAP(x(1:4),f);
	elseif(strcmp(fitType,'eq5'))
		fitAP = @(x,f) eq5(x(1),x(2),x(3),2*pi*f);
		full_model = @(f,x) fitAP(x(1:3),f) + fitPeaks(x(4:end),f);
		AP_model = @(f,x) fitAP(x(1:3),f);
	end
end
function y = eq1(tau1,tau2,A1,A2,f)
	y = log10(10^A1./(1+tau1^2*f.^2) + 10^A2./(1+tau2^2*f.^2));
end
function y = eq5(tau1,offset,mag,f)
	y0 = tau1 ./ (1+tau1^2*f.^2);
	y = mag + log10(exp(offset)+y0);
end
function y = eq6(tau1,tau2,offset,mag,f)
	y0 = (tau1-tau2)^2 ./ ((1+tau1^2*f.^2).*(1+tau2^2*f.^2));
	y = mag + log10(exp(offset)+y0);
end
function y = fitPeaks(x,f)
	y = zeros(size(f));
	for i = 1:3:length(x)
		y = y+fitGauss(x(i),x(i+1),x(i+2),f);
	end
end
function y = fitGauss(mu,hgt,sig,f)
	y = hgt * exp(-(f-mu).^2 / (2*sig^2));
end
