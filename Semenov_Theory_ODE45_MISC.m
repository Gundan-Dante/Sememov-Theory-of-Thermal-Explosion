clc; clearvars; close all; format long;
function new_figure = Concentration_Vs_Temperature(Solution_Matrix)
    new_figure = figure("Name","sugma");
    plot(Solution_Matrix(:, 1), Solution_Matrix(:, 2));
    xlabel('\theta');
    ylabel('\nu');
    grid on;
    return
end
function new_figure = integration_ode(alpha, beta, gamma, Theta_null, C_null)
    Y0 = [Theta_null, C_null];
    odefun = @(tau, Y) [
        1/gamma * (Y(2) * exp(Y(1)/(1+beta * Y(1))) - alpha * Y(1));
        -Y(2) * exp(Y(1)/(1+beta * Y(1)))
        ];
    [tau, Y] = ode45(odefun, [0, 1], Y0);
    new_figure = Concentration_Vs_Temperature(Y)
    figure_2 = figure()
    subplot(2, 1, 1);
    plot(tau, Y(:, 1));
    xlabel('\tau');
    ylabel('\theta');
    title(['\alpha=', num2str(alpha), ', \beta=', num2str(beta), ', \gamma=', num2str(gamma)]);
    grid on;
    subplot(2, 1, 2);
    plot(tau, Y(:, 2));
    xlabel('\tau');
    ylabel('\nu');
    grid on;
end
% Specify the parameters:
Theta_0 = 1;
C_0 = 1;
par1 = input('choose alpha from set : [2, 3]>>');
% par2 = input('choose beta from set : [0.01, 0.1] or [-0.01, -0.1]>>');
% par3 = input('choose gamma from set : [0.01, 0.1]>>');
bear = integration_ode(par1, 0.025, 0.02, Theta_0, C_0);
% output_folder = "C:\Users\aseit\OneDrive\Рабочий стол\For  Matlab Figures"; %unnecessary part needed to save the graph
% filename = fullfile(output_folder, 'for bear.png');
% saveas(bear, filename);
% 2.329400539398193
% 2.329400444030761
