%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for the modulo denoising via kNN and unwrapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

addpath(genpath('lib'))  
% choose the method below
method = 'localPoly';
disp(method)
disp('Started')

%% Parameters
sigma = 0.12; % noise level 

n = 600;

l = 2; beta = 2.4; C_lp = 0.1;
h = C_lp*(log(n)./n).^(beta/(2*beta+1)); %% length of rectangular window

%% Ground truth
ff = @(x) 4 + 4.*x .* cos(2*pi*x) .^2 - 2.*sin(2*pi*x).^2 + 0.7; 

a = 0;b = 1; 

x = (a:((b-a)/(n-1)):b)'; % nx1 vector  
f_clean = ff(x); 
f_mod1_clean = mod(f_clean,1); 

f_noise = f_clean + sigma * randn(n,1); 
y = mod(f_noise,1);  
z = exp(1i*2*pi*y);  % nx1 vector
        
gest_localPoly = localPoly_denoise(z,x,h,l);
gest_localPoly_proj = project_manifold(gest_localPoly);
f_mod1_denoised = extract_modulo(gest_localPoly_proj);  
    

%% Unwrapping
f_unwrapped = unwrap_1D(f_mod1_denoised);
y_noisy_unwrapped = unwrap_1D(y);


%% Plot the results
subplot(6,1,1);
plot(x,f_clean); % Unwrapped Clean mod 1 samples
title('Clean ground truth','Interpreter','latex') 

subplot(6,1,2);
plot(x,f_mod1_clean); % Clean mod 1 samples
title('Clean $\bmod 1$','Interpreter','latex') 

subplot(6,1,3);
plot(x,y); % Noisy mod 1 samples
title('Noisy $\bmod 1$','Interpreter','latex') 

subplot(6,1,4);
plot(x,f_mod1_denoised); % denoised mod 1 samples
title('Denoised $\bmod 1$','Interpreter','latex') 

subplot(6,1,5);
plot(x,f_unwrapped); % Unwrapped Clean mod 1 samples
title('Unwrapped denoised','Interpreter','latex') 

subplot(6,1,6);
plot(x,y_noisy_unwrapped); % Unwrapped Clean mod 1 samples
title('Unwrapped noisy','Interpreter','latex') 

saveas(gcf,'figures/ex2','epsc')

