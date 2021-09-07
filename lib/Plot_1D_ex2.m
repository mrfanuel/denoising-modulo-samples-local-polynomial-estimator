
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results ex2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Plot wrap-around RMSE %%%%%%%%%%
figure;

% ----- localPoly -----
method = 'localPoly'
nameFile = strcat(strcat('/ex2_paper_err_mod1_',method),'.mat');
nameFile = strcat(folder,nameFile);
load(nameFile)

plot(range_n, err_wrap_around, '-bo', 'MarkerSize', 4, 'markerfacecolor','b','DisplayName',method);hold on;
errorbar(range_n, err_wrap_around, std_err_wrap_around,'b' )
hold on;
% unwrapped noisy
plot(range_n, err_wrap_around_noisy, '-k<', 'MarkerSize', 4, 'markerfacecolor','k','DisplayName','noisy')
errorbar(range_n, err_wrap_around_noisy, std_err_wrap_around_noisy,'k' )

xlabel('$n$','Interpreter','latex', 'FontSize', 25)
ylabel('Wrap around MSE','Interpreter','latex', 'FontSize', 25)
hold on;

% ----- TRS -----
method = 'TRS';
nameFile = strcat(strcat('/ex2_paper_err_mod1_',method),'.mat');
nameFile = strcat(folder,nameFile);
load(nameFile)

plot(range_n, err_wrap_around, '-ro', 'MarkerSize', 4, 'markerfacecolor','r','DisplayName',method);hold on;
errorbar(range_n, err_wrap_around, std_err_wrap_around,'r' )
hold on;


% ----- UCQP -----
method = 'UCQP';
nameFile = strcat(strcat('/ex2_paper_err_mod1_',method),'.mat');
nameFile = strcat(folder,nameFile);
load(nameFile)

plot(range_n, err_wrap_around, '-go', 'MarkerSize', 4, 'markerfacecolor','g','DisplayName',method);hold on;
errorbar(range_n, err_wrap_around, std_err_wrap_around,'g' )

% legend

legend


%%%%%%%%%%%% Plot  RMSE %%%%%%%%%%%%
figure;

% ----- localPoly -----
method = 'localPoly'
nameFile = strcat(strcat('/ex2_paper_err_unwrapped_',method),'.mat');
nameFile = strcat(folder,nameFile);
load(nameFile);

plot(range_n, err_unwrapped, '-bo', 'MarkerSize', 4, 'markerfacecolor','b','DisplayName',method);hold on;
errorbar(range_n, err_unwrapped, std_err_unwrapped ,'b')

hold on; 
% unwrapped noisy
plot(range_n, err_unwrapped_noisy, '-k<', 'MarkerSize', 4, 'markerfacecolor','k','DisplayName','noisy');hold on;
errorbar(range_n, err_unwrapped_noisy, std_err_unwrapped_noisy,'k')

xlabel('$n$','Interpreter','latex', 'FontSize', 25)
ylabel('MSE','Interpreter','latex', 'FontSize', 25)
hold on;

% ----- TRS -----
method = 'TRS'
nameFile = strcat(strcat('/ex2_paper_err_unwrapped_',method),'.mat');
nameFile = strcat(folder,nameFile);
load(nameFile);

plot(range_n, err_unwrapped, '-ro', 'MarkerSize', 4, 'markerfacecolor','r','DisplayName',method);hold on;
errorbar(range_n, err_unwrapped, std_err_unwrapped ,'r')
hold on;

% ----- UCQP -----
method = 'UCQP'
nameFile = strcat(strcat('/ex2_paper_err_unwrapped_',method),'.mat');
nameFile = strcat(folder,nameFile);
load(nameFile);

plot(range_n, err_unwrapped, '-go', 'MarkerSize', 4, 'markerfacecolor','g','DisplayName',method);hold on;
errorbar(range_n, err_unwrapped, std_err_unwrapped ,'g')

% legend

legend