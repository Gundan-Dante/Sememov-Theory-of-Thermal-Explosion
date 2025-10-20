clc; clearvars; format long;

% function new_figure = Temperature_Vs_Concentration(Solution_Matrix)
%     disp(class(Solution_Matrix));
%     new_figure = figure("Name",'ligma');
%     plot(flip(Solution_Matrix(:, 2)), Solution_Matrix(:, 1));
%     xlabel('\nu');
%     ylabel('\theta');
%     grid on;
%     return
% end

% function new_figure = Concentration_Vs_Temperature(Solution_Matrix)
%     new_figure = figure("Name","sugma");
%     plot(Solution_Matrix(:, 1), Solution_Matrix(:, 2));
%     xlabel('\theta');
%     ylabel('\nu');
%     grid on;
%     return
% end

% function Alpha_Crit(Parameter_Matrix, Ini_Val_Matrix)
%     %Unpacking:
%     [alphaLB, alphaUB, bet, gam] = deal(Parameter_Matrix{:});
%     for n = 1:20
%         alpha_av = (alphaUB + alphaLB)/2;
%         odeballs = @(tau, Y) [
%             1/gam * ( Y(2) * exp(Y(1)/(1+bet * Y(1))) - alpha_av * Y(1));
%             -Y(2) * exp(Y(1)/(1+bet * Y(1)));
%         ]
%         [tau, Y] = ode45(odeballs, [0, 1], Ini_Val_Matrix);
%         %check the sol matrix values through slicing
%         %implement checking for zero values here
%         zero_count = 0;
%         for nu_instance = Y(:,2)'
%             if nu_instance < 0.0001
%                 zero_count = zero_count + 1;
%             end
%         end
%         zero_count
%         if zero_count > 0.05 * length(Y(:, 2))
%             alphaLB = alpha_av
%         else
%             alphaUB = alpha_av
%         end
%     end
%     Concentration_Vs_Temperature(Y)
%     Y(:,2);
% end

function Alpha_Crit_Optimized_Gamma(Parameter_Matrix, Ini_Val_Matrix, ORDER)
    %Unpacking:
    [alphaLB_gamma, alphaUB_gamma, bet, gam] = deal(Parameter_Matrix{:});
    tol = 40;
    for n = 1:(ORDER-1)
        % disp('calculation for the following order >>' + string(n))
        span = linspace(alphaLB_gamma, alphaUB_gamma, 11);
        for alpha = span
            odeballs = @(tau, Y) [
                1/gam * ( Y(2) * exp(Y(1)/(1+bet * Y(1))) - alpha * Y(1));
                -Y(2) * exp(Y(1)/(1+bet * Y(1)));
                ];
            [tau, Y] = ode45(odeballs, [0, 1], Ini_Val_Matrix);
            val = max(Y(:, 1));
            % abs(val - 1/gam)
            if abs(val - 50) > tol
                alphaUB_gamma = alpha;
                break
            end
            alphaLB_gamma = alpha;
        end
    end
    disp('Here is the gamma method alpha boundary')
    alphaUB_gamma
end

function Alpha_Crit_Optimized_NARUTO(Parameter_Matrix, Ini_Val_Matrix, ORDER)
    %Unpacking:
    [alphaLB, alphaUB, bet, gam] = deal(Parameter_Matrix{:});
    tol_zero = 0.01;
    for n = 1:(ORDER-1)
        % disp('calculation for the following order >>' + string(n))
        span = linspace(alphaLB, alphaUB, 11);
        for alpha = span
            odeballs = @(tau, Y) [
                1/gam * ( Y(2) * exp(Y(1)/(1+bet * Y(1))) - alpha * Y(1));
                -Y(2) * exp(Y(1)/(1+bet * Y(1)));
                ];
            [tau, Y] = ode45(odeballs, [0, 1], Ini_Val_Matrix);
            zero_count = 0;
            for conc_instance = Y(:, 2)'
                if conc_instance < tol_zero
                    zero_count = zero_count + 1;
                end
            end
            if zero_count > 0.05 * length(Y(:, 2)')
                alphaLB = alpha;
            else
                alphaUB = alpha;
                break
            end
        end
    end
    disp('Zero counting method alpha upper boundary ')
    alphaUB
end

% function new_figure = integration_ode(alpha, beta, gamma, Theta_null, C_null)
%     Y0 = [Theta_null, C_null];
%     odefun = @(tau, Y) [
%         1/gamma * (Y(2) * exp(Y(1)/(1+beta * Y(1))) - alpha * Y(1));
%         -Y(2) * exp(Y(1)/(1+beta * Y(1)))
%         ];
%     [tau, Y] = ode45(odefun, [0, 1], Y0);
%     [val, idt] = max(Y(:, 1));
%     % peak_instance = tau(idt)
%     % peak_temp = val
%     % subplot(2, 1, 1);
%     % plot(tau, Y(:, 1));
%     % xlabel('\tau');
%     % ylabel('\theta');
%     % title(['\alpha=', num2str(alpha), ', \beta=', num2str(beta), ', \gamma=', num2str(gamma)]);
%     % text(0.5, 0.7, sprintf('$\\theta_{\\max} = %.2f,\\ \\tau = %.2f$', peak_temp, peak_instance), ...
%     % 'Interpreter', 'latex');
%     % grid on;
%     % subplot(2, 1, 2);
%     % plot(tau, Y(:, 2));
%     % xlabel('\tau');
%     % ylabel('\nu');
%     % grid on;
%     %Temperature_Vs_Concentration(Y);
%     new_figure = Concentration_Vs_Temperature(Y);
%     clc;
%     Y
% end

% Specify the parameters:
% Theta_0 = 1;
% C_0 = 1;
% key = true;
% while key == true
%     par1 = input('choose alpha from set : [2, 3]>>');
%     par2 = input('choose beta from set : [0.01, 0.1] or [-0.01, -0.1]>>');
%     par3 = input('choose gamma from set : [0.01, 0.1]>>');
%     figurus = integration_ode(par1, par2, par3, Theta_0, C_0);
%     % output_folder = "C:\Users\aseit\OneDrive\Рабочий стол\For  Matlab Figures"; %unnecessary part needed to save the graph
%     % filename = fullfile(output_folder, 'ode45 - crit.png');
%     % saveas(figurus, filename);
%     space = input('If you will to continue the program, press ENTER, to terminate it, type in 0>> ');
%     if space == 0;
%         key = false
%     end
% end

% Alpha crit determination
Parameter_Matrix_To_Handle = {2.3, 2.4, 0.025, 0.02};
Ini_Val_Matrix_To_Handle = [1, 1];
order_handle = input('Input the order of precision for critical alpha value >>');
Alpha_Crit_Optimized_Gamma(Parameter_Matrix_To_Handle, Ini_Val_Matrix_To_Handle, order_handle)
Alpha_Crit_Optimized_NARUTO(Parameter_Matrix_To_Handle, Ini_Val_Matrix_To_Handle, order_handle)

% NOTES FOR PROF -------------------------------------------------------

% This simple code finds the precise value of alpha parameter (see the
% semenov theory document) needed to set the premixed fuel in
% explosive (dangerous, hence unwanted) or steady burning regime (wanted).

% The model utilizes 0D burning numerical solution.

% Two decision methods were implemented - 1/gamma and zero-counting
% strategy

% Both algorithms fail after magnitude 7: they simply don't distinguish
% between explosive and non-explosive regimes, since it is me who defines
% the tolerances for desicion making process.

% So, main question, what is the approximate proximity of theta_max to
% 1/gamma that actually defines the explosive behavior?

% Or is it the way that the graph looks that defines it?
