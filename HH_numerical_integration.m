close all;

% Numerically integrates Hodgkin-Huxley model and plots solution for V(t).
% Options available for several different stimuli. Also plots gating
% variables in a lower subplot to see how they evolve in time with
% voltage/the action potential.


% gating function parameters taken from Foundations of Computational Neuroscience pg
% 32: same values Hodgkin and Huxley got experimentally from squid
% constant values on pg 32 of this textbook


%% Resources 


% assumed Ca concentration decays exponentially
% assume A current channel has 3 activation subunits / one inactivation
% assume A current inactivation gate behaves similarly to V-Na inactivation


% 1) Common Channel conductances
% https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104094&ver=0&trm=conductance+of+Ca+channels&org=


%%
close all; clear all; clc;
% define constants

v_rest = -65;

% max conductances
g_na = 120;
g_k = 36;
g_l = 0.3;
% g_ca_gated_k = 1.25;
g_ca_gated_k = 0;
g_ca = 1/500;
g_A_current = 0;



% Equilibrium Potentials
E_k = -77;
E_na = 50;
E_l = -54.4;
E_ca = 134; % https://www.cvphysiology.com/Arrhythmias/A007

% Calcium-Gated time constant: describes Ca2+ leaving from cell;
tau_ca = 20;
tau_A = 0.5;


% pg 17 of Foundations of Comp Neuroscience + pg 520 of Hodgkin-Huxley
% paper
c_m = 1; 
phi = 1; % scales for temperature

% input current (stimulus) - Give some function I(t)

% step current:
% I_ex =  @(t) 25 * (t >= 75); % Ca bursting
I_ex =  @(t) 10 * (t >= 25); % A current

x0 = [v_rest 0.3 0.1 0.6 0 0]; % initial conditions; approximated off graph pg 33 FCN
tspan = [0 100];

% c is the gating variable for Ca gate potassium channels
% x = [V n m h Ca AcurrentActivation]
diff_eqns = @(t,x) [
1/c_m .*(I_ex(t) - g_na * x(3)^3 * x(4) * (x(1) - E_na) - g_k*x(2)^4 * (x(1) - E_k) - g_l * (x(1) - E_l) - g_ca_gated_k * x(5)*(x(1) - E_k) - g_A_current * x(6)^3 * x(4) * (x(1) - E_k) );
phi * (alpha_n(x(1)) * (1-x(2)) - beta_n(x(1)) * x(2));
phi * (alpha_m(x(1)) * (1-x(3)) - beta_m(x(1)) * x(3));
phi * (alpha_h(x(1)) * (1-x(4)) - beta_h(x(1)) * x(4));

% exponential decay in Ca concentration plus bursts of Ca
% 0 V-Ca open at v_rest and all open at 30

 phi * (x(2)^4 .* g_ca * (E_ca - x(1)) - 1/tau_ca * x(5));
  
 1/tau_A *((x(1) > (v_rest + 4)) -  x(6))];

%    0.05 * ( (x(1) > (v_rest + 3)) - 1/tau_A * x(6) )];




[t_vals, y_vals] = ode45(diff_eqns, tspan, x0);

figure(1)
subplot(211)
plot(t_vals, y_vals(:, 1), 'r-'), hold on;
plot(t_vals, ones(1, length(t_vals)) * -55, 'k--'), hold on;
plot(t_vals, ones(1, length(t_vals)) * (v_rest), 'b--');
ylabel("Membrane Potential (mV)");
xlabel("time (ms)");
legend({"V(t)", "threshold", 'V_{rest}'});
title("Bursting Action Potential under Influence of Ca^{2+}-gated K^+ Channels");
xlim(tspan)
set(gca, 'fontsize', 15);
ylim([-80, 75])

subplot(212)
I_vals = I_ex(t_vals);
plot(t_vals, I_vals, 'm');
ylabel("I_{ext}(t)")
title("Applied Current I(t)");
ylim([min(I_vals)-1 1.25*max(I_vals)]);
set(gca, 'fontsize', 15);

figure(2)
subplot(211);
plot(t_vals, y_vals(:, 2)), hold on;
plot(t_vals, y_vals(:, 3)), hold on;
plot(t_vals, y_vals(:, 4)), hold on;
plot(t_vals, y_vals(:, 6));
title("Gating Variables")
legend({"n = K^+ Activation", "m = Na^+ Activation", "h = Na^+ Inactivation", "k = A-current activation"});
xlim(tspan);
ylim([0,1]);

set(gca, 'fontsize', 15);

subplot(212)
plot(t_vals, y_vals(:, 5));
title('Calcium Concentration')
xlabel('time (ms)')
ylabel('[Ca^{2+}]')
set(gca, 'fontsize', 15);
ylim([0,1])


gcf