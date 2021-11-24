%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for the modulo denoising via kNN and unwrapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
clc

addpath(genpath('lib'))  

% choose the method below
methods = ["localPoly" "kNN" "TRS" "UCQP"];
sigma = 0.12; % noise level 
n_MC = 5; % number of Monte-Carlo runs

for method = methods
    disp(method)
    disp('Started')
    
    %% Parameters
    
    range_n = 100:100:1000; % range of n values
    c = 0.04;%0.07;
    
    if strcmp(method,'kNN')
        C_kNN = 0.09; %0.07;%
        k = ceil(C_kNN*(range_n.^(2/3)).*(log(range_n).^(1/3))); %% number of neigbours
    elseif strcmp(method,'UCQP')
        C_UCQP =  c;
        lambda = C_UCQP*(range_n.^(10/3)).^(1/4);
    elseif strcmp(method,'TRS')
        C_TRS =  c;
        lambda = C_TRS*(range_n.^(10/3)).^(1/4);
    elseif strcmp(method,'localPoly')
        %l = 2; beta = 2.4; C_lp = 0.05;
        l = 2; beta = 2.4; C_lp = 0.1;

        h = C_lp*(log(range_n)./range_n).^(beta/(2*beta+1)); %% length of rectangular window
    end
    
    
  
    %% Initialization mean errors and std's
    err_wrap_around = zeros(size(range_n)); std_err_wrap_around = zeros(size(range_n));
    err_wrap_around_noisy = zeros(size(range_n));std_err_wrap_around_noisy = zeros(size(range_n));
    err_unwrapped = zeros(size(range_n)); std_err_unwrapped = zeros(size(range_n));
    err_unwrapped_noisy = zeros(size(range_n)); std_err_unwrapped_noisy = zeros(size(range_n));
    
    
    err_wrap_around_temp = ones(n_MC,1);
    err_wrap_around_noisy_temp = ones(n_MC,1);
    err_unwrapped_temp = ones(n_MC,1);
    err_unwrapped_noisy_temp = ones(n_MC,1);
    
    
    %% Ground truth
    ff = @(x) 4 + 4.*x .* cos(2*pi*x) .^2 - 2.*sin(2*pi*x).^2 + 0.7; 
    
    a = 0;b = 1; 
    
    for index = 1:length(range_n)
    
        n = range_n(index);
    
        x = (a:((b-a)/(n-1)):b)'; % nx1 vector  
        f_clean = ff(x); 
        f_mod1_clean = mod(f_clean,1); 
        
            
        if strcmp(method,'kNN')
            nb_nb = k(index);
        else
            %% For graph-based methods
            r = 1;
            % Form graph G with "connectivity radius" = r and its Laplacian
            A = zeros(n);
            for i=1:n
                for j=1:n
                    if i~=j & abs(i-j) <= r
                        A(i,j)=1;
                    end
                end
            end
            d = A*ones(n,1);
            L = diag(d) - A;
        end
    
        for iter= 1:n_MC
            disp(iter)
            f_noise = f_clean + sigma * randn(n,1); 
            y = mod(f_noise,1);  
            z = exp(1i*2*pi*y);  % nx1 vector
            
            %% Denoising
            if strcmp(method,'kNN')
                gest_kNN = kNN_denoise(z,x,nb_nb);
                gest_kNN_proj = project_manifold(gest_kNN);
                f_mod1_denoised = extract_modulo(gest_kNN_proj);
                
            elseif strcmp(method,'UCQP')
                reg_param = lambda(index);
                gest_ucqp = UCQP_denoise(z,L,reg_param,n);
                gest_ucqp_proj = project_manifold(gest_ucqp);
                f_mod1_denoised = extract_modulo(gest_ucqp_proj);      
            elseif strcmp(method,'TRS')
                reg_param = lambda(index);
                gest_trs = TRS_denoise(z,L,reg_param,n);
                gest_trs_proj = project_manifold(gest_trs);
                f_mod1_denoised = extract_modulo(gest_trs_proj);
            elseif strcmp(method,'localPoly')
                gest_localPoly = localPoly_denoise(z,x,h(index),l);
                gest_localPoly_proj = project_manifold(gest_localPoly);
                f_mod1_denoised = extract_modulo(gest_localPoly_proj);  
            end
         
            err_wrap_around_temp(iter) = MS_wrap_around_error(f_mod1_denoised, f_mod1_clean);
            err_wrap_around_noisy_temp(iter) = MS_wrap_around_error(y, f_mod1_clean);
    
            %% Unwrapping
            f_unwrapped = unwrap_1D(f_mod1_denoised);
    
            err_unwrapped_temp(iter) = MS_error(f_unwrapped, f_clean);
    
            y_noisy_unwrapped = unwrap_1D(y);
            err_unwrapped_noisy_temp(iter) = MS_error(y_noisy_unwrapped, f_clean);
    
        end
    
    
        err_wrap_around(index) = mean(err_wrap_around_temp);
        std_err_wrap_around(index) = std(err_wrap_around_temp);
    
        err_wrap_around_noisy(index) = mean(err_wrap_around_noisy_temp);
        std_err_wrap_around_noisy(index) = std(err_wrap_around_noisy_temp);
    
        err_unwrapped(index) = mean(err_unwrapped_temp);
        std_err_unwrapped(index) = std(err_unwrapped_temp);
    
        err_unwrapped_noisy(index) = mean(err_unwrapped_noisy_temp);
        std_err_unwrapped_noisy(index) = std(err_unwrapped_noisy_temp);
    
    end

    %%%%%%%%%%%%%%% RMSE Wrap around error %%%%%%%%%%%%%%% 
    %close all;
    figure(1);
    if strcmp(method,"localPoly")
        plot(range_n, err_wrap_around, '-bo', 'MarkerSize', 4, 'markerfacecolor','b','DisplayName','localPoly');hold on;
        errorbar(range_n, err_wrap_around, std_err_wrap_around,'b', 'HandleVisibility','off' )
        legend('FontSize',20,'Interpreter','latex')
        legend off; legend show;

    elseif strcmp(method,"kNN")
        plot(range_n, err_wrap_around, '-rs', 'MarkerSize', 4, 'markerfacecolor','r','DisplayName','kNN');hold on;
        errorbar(range_n, err_wrap_around, std_err_wrap_around,'r', 'HandleVisibility','off' )
        legend('FontSize',20,'Interpreter','latex')
        legend off; legend show;

    elseif strcmp(method,"TRS")
        plot(range_n, err_wrap_around, '-cd', 'MarkerSize', 4, 'markerfacecolor','c','DisplayName','TRS');hold on;
        errorbar(range_n, err_wrap_around, std_err_wrap_around,'c', 'HandleVisibility','off' )
        legend('FontSize',20,'Interpreter','latex')
        legend off; legend show;

    elseif strcmp(method,"UCQP")
        plot(range_n, err_wrap_around, '-m^', 'MarkerSize', 4, 'markerfacecolor','m','DisplayName','UCQP');hold on;
        errorbar(range_n, err_wrap_around, std_err_wrap_around,'m', 'HandleVisibility','off' )
        legend('FontSize',20,'Interpreter','latex')
        legend off; legend show;

    end

    xlim([95,1005])
    xlabel('$n$','Interpreter','latex', 'FontSize', 25)
    ylabel('Wrap around RMSE','Interpreter','latex', 'FontSize', 25)
    
    %%%%%%%%%%%%%%% RMSE unwrapped %%%%%%%%%%%%%%% 
    
    %close all;
    figure(2);
    if strcmp(method,"localPoly")
        plot(range_n, err_unwrapped, '-bo', 'MarkerSize', 4, 'markerfacecolor','b','DisplayName','localPoly');hold on;
        errorbar(range_n, err_unwrapped, std_err_unwrapped ,'b', 'HandleVisibility','off')
        legend('FontSize',20,'Interpreter','latex')
        legend off; legend show;
    elseif strcmp(method,"kNN")
        plot(range_n, err_unwrapped, '-rs', 'MarkerSize', 4, 'markerfacecolor','r','DisplayName','kNN');hold on;
        errorbar(range_n, err_unwrapped, std_err_unwrapped ,'r', 'HandleVisibility','off')
        legend('FontSize',20,'Interpreter','latex')
        legend off; legend show;
    elseif strcmp(method,"TRS")
        plot(range_n, err_unwrapped, '-cd', 'MarkerSize', 4, 'markerfacecolor','c','DisplayName','TRS');hold on;
        errorbar(range_n, err_unwrapped, std_err_unwrapped ,'c', 'HandleVisibility','off')
        legend('FontSize',20,'Interpreter','latex')
        legend off; legend show;
        
    elseif strcmp(method,"UCQP")
        plot(range_n, err_unwrapped, '-m^', 'MarkerSize', 4, 'markerfacecolor','m','DisplayName','UCQP');hold on;
        errorbar(range_n, err_unwrapped, std_err_unwrapped ,'m', 'HandleVisibility','off')
        legend('FontSize',20,'Interpreter','latex')
        legend off; legend show;
    end
    xlim([95,1005])
    xlabel('$n$','Interpreter','latex', 'FontSize', 25)
    ylabel('RMSE','Interpreter','latex', 'FontSize', 25)
    
end



%% Plot the results for the paper

%%%%%%%%%%%%%%% RMSE Wrap around error %%%%%%%%%%%%%%% 
%close all;
figure(1);
plot(range_n, err_wrap_around_noisy, '-k<', 'MarkerSize', 4, 'markerfacecolor','k','DisplayName','noisy')
errorbar(range_n, err_wrap_around_noisy, std_err_wrap_around_noisy,'k' , 'HandleVisibility','off')

legend off; legend show;
xlabel('$n$','Interpreter','latex', 'FontSize', 25)
ylabel('Wrap around RMSE','Interpreter','latex', 'FontSize', 25)
legend('FontSize',20,'Interpreter','latex')
hold off
saveas(gcf,'figures/wraparound_RMSE_ex2','epsc')
%%%%%%%%%%%%%%% RMSE unwrapped %%%%%%%%%%%%%%% 

%close all;
figure(2);
%plot(range_n, err_unwrapped_noisy, '-k<', 'MarkerSize', 4, 'markerfacecolor','k','DisplayName','noisy');hold on;
%errorbar(range_n, err_unwrapped_noisy, std_err_unwrapped_noisy,'k' )

xlabel('$n$','Interpreter','latex', 'FontSize', 25)
ylabel('RMSE','Interpreter','latex', 'FontSize', 25)
legend('FontSize',20,'Interpreter','latex')

hold off
saveas(gcf,'figures/RMSE_ex2','epsc')
disp('ended')

