%% Comparativa J–V y P–V: modelo simple vs. realista
% curvas con valores iniciales y con valores óptimos encontrados
clc; clear; close all;

% --- Parámetros y constantes ---
q      = 1.602176634e-19;   % C
k      = 1.3880649e-23;     % J/K
T      = 300;               % K
perm   = 1.035918e-10;      % F/cm
Dn     = 2;                 % cm^2/s
tau_n  = 371e-6;            % s
Dp     = 11.6;              % cm^2/s
tau_p  = 3710e-6;           % s
JL     = 50e-3;             % A/cm²
nint   = 9.696e9;           % cm^-3
Sf_eff = 3e4;               % cm/s (frontal, lado n)    
Sbsf   = 100;               % cm/s (trasera, lado p)    

% Conjunto de parámetros 
% Iniciales de la tabla
Wn0 = 0.2e-4;    Nd0 = 1e20;     Wp0 = 179.8e-4;  Na0 = 4.7e15;

% Óptimos reportados por la optimización:
%  para modelo simple (máximo en V ≈ 0.511 V)
WnS = 7.500e-05; NdS = 5.000e20; WpS = 2.000e-02; NaS = 1.000e16; VoptS = 0.511;
%  para modelo realista (máximo en V ≈ 0.520 V)
WnR = 7.500e-05; NdR = 5.000e20; WpR = 2.000e-02; NaR = 1.000e16; VoptR = 0.520;

% Rango de voltaje 
Vmax = 0.70;
V    = linspace(0, Vmax, 2001);

%% Curvas con parámetros INICIALES 
% Modelo simple
[J0_s, P0_s] = JV_simple(V, Wn0, Nd0, Wp0, Na0, q,k,T,perm,Dn,tau_n,Dp,tau_p,JL,nint);
% Modelo realista
[J0_r, P0_r] = JV_realista(V, Wn0, Nd0, Wp0, Na0, q,k,T,perm,Dn,tau_n,Dp,tau_p,JL,nint,Sf_eff,Sbsf);

%% Curvas con parámetros ÓPTIMOS 
% Modelo simple (óptimos para simple)
[JS_opt, PS_opt] = JV_simple(V, WnS, NdS, WpS, NaS, q,k,T,perm,Dn,tau_n,Dp,tau_p,JL,nint);
% Modelo realista (óptimos para realista)
[JR_opt, PR_opt] = JV_realista(V, WnR, NdR, WpR, NaR, q,k,T,perm,Dn,tau_n,Dp,tau_p,JL,nint,Sf_eff,Sbsf);

%% Cálculo de máximos y valores en V óptimo reportado
% Simple
[~, idxP0s] = max(P0_s);   [~, idxPS] = max(PS_opt);
J0s_atVopt  = interp1(V, J0_s, VoptS, 'linear', 'extrap');
Js_atVopt   = interp1(V, JS_opt, VoptS, 'linear', 'extrap');
P0s_atVopt  = interp1(V, P0_s, VoptS, 'linear', 'extrap');
Ps_atVopt   = interp1(V, PS_opt, VoptS, 'linear', 'extrap');

% Realista
[~, idxP0r] = max(P0_r);   [~, idxPR] = max(PR_opt);
J0r_atVopt  = interp1(V, J0_r, VoptR, 'linear', 'extrap');
Jr_atVopt   = interp1(V, JR_opt, VoptR, 'linear', 'extrap');
P0r_atVopt  = interp1(V, P0_r, VoptR, 'linear', 'extrap');
Pr_atVopt   = interp1(V, PR_opt, VoptR, 'linear', 'extrap');

%% Máscaras solo cuadrante positivo
maskJ0s  = J0_s   > 0;  maskJS   = JS_opt > 0;
maskP0s  = P0_s   > 0;  maskPS   = PS_opt > 0;
maskJ0r  = J0_r   > 0;  maskJR   = JR_opt > 0;
maskP0r  = P0_r   > 0;  maskPR   = PR_opt > 0;

%% FIGURAS: MODELO SIMPLE 
% J–V (simple)
figure('Name','J-V (modelo simple)','Color','w');
plot(V(maskJ0s),  1e3*J0_s(maskJ0s),  'b-', 'LineWidth',2); hold on;
plot(V(maskJS),   1e3*JS_opt(maskJS), 'r--','LineWidth',2);
plot(VoptS, 1e3*Js_atVopt, 'ro', 'MarkerFaceColor','r');    % marca V óptimo
xlabel('Voltaje V (V)'); ylabel('Densidad de corriente (mA/cm²)');
title('J–V (modelo simple): inicial vs. óptimo'); grid on; box on;
legend('Inicial','Óptimo','V_{opt} (simple)','Location','southwest');
xlim([0 Vmax]); ylim([0 inf]);

% P–V (simple)
figure('Name','P-V (modelo simple)','Color','w');
plot(V(maskP0s),  1e3*P0_s(maskP0s),  'b-', 'LineWidth',2); hold on;
plot(V(maskPS),   1e3*PS_opt(maskPS), 'r--','LineWidth',2);
plot(VoptS, 1e3*Ps_atVopt, 'ro', 'MarkerFaceColor','r');    % marca V óptimo
xlabel('Voltaje V (V)'); ylabel('Potencia (mW/cm²)');
title('P–V (modelo simple): inicial vs. óptimo'); grid on; box on;
legend('Inicial','Óptimo','V_{opt} (simple)','Location','southwest');
xlim([0 Vmax]); ylim([0 inf]);

fprintf('--- MODELO SIMPLE ---\n');
fprintf('P_{max} inicial  = %.4e W/cm² en V = %.4f V\n', P0_s(idxP0s), V(idxP0s));
fprintf('P_{max} óptimo   = %.4e W/cm² en V = %.4f V\n', PS_opt(idxPS), V(idxPS));
fprintf('En V_{opt}=%.3f V -> P_inicial=%.4e, P_opt=%.4e | J_inicial=%.4e, J_opt=%.4e (A/cm²)\n', ...
        VoptS, P0s_atVopt, Ps_atVopt, J0s_atVopt, Js_atVopt);

%% FIGURAS: MODELO REALISTA 
% J–V (realista)
figure('Name','J-V (modelo realista)','Color','w');
plot(V(maskJ0r),  1e3*J0_r(maskJ0r),  'b-', 'LineWidth',2); hold on;
plot(V(maskJR),   1e3*JR_opt(maskJR), 'r--','LineWidth',2);
plot(VoptR, 1e3*Jr_atVopt, 'ro', 'MarkerFaceColor','r');    % marca V óptimo
xlabel('Voltaje V (V)'); ylabel('Densidad de corriente (mA/cm²)');
title('J–V (modelo realista): inicial vs. óptimo'); grid on; box on;
legend('Inicial','Óptimo','V_{opt} (realista)','Location','southwest');
xlim([0 Vmax]); ylim([0 inf]);

% P–V (realista)
figure('Name','P-V (modelo realista)','Color','w');
plot(V(maskP0r),  1e3*P0_r(maskP0r),  'b-', 'LineWidth',2); hold on;
plot(V(maskPR),   1e3*PR_opt(maskPR), 'r--','LineWidth',2);
plot(VoptR, 1e3*Pr_atVopt, 'ro', 'MarkerFaceColor','r');    % marca V óptimo
xlabel('Voltaje V (V)'); ylabel('Potencia (mW/cm²)');
title('P–V (modelo realista): inicial vs. óptimo'); grid on; box on;
legend('Inicial','Óptimo','V_{opt} (realista)','Location','southwest');
xlim([0 Vmax]); ylim([0 inf]);

fprintf('--- MODELO REALISTA ---\n');
fprintf('P_{max} inicial  = %.4e W/cm² en V = %.4f V\n', P0_r(idxP0r), V(idxP0r));
fprintf('P_{max} óptimo   = %.4e W/cm² en V = %.4f V\n', PR_opt(idxPR), V(idxPR));
fprintf('En V_{opt}=%.3f V -> P_inicial=%.4e, P_opt=%.4e | J_inicial=%.4e, J_opt=%.4e (A/cm²)\n', ...
        VoptR, P0r_atVopt, Pr_atVopt, J0r_atVopt, Jr_atVopt);

%% Funciones locales 
function [J, P] = JV_simple(V, Wn, Nd, Wp, Na, q,k,T,perm,Dn,tau_n,Dp,tau_p,JL,nint)
    Ln   = sqrt(Dn * tau_n);
    Lp   = sqrt(Dp * tau_p);
    Vint = (k*T/q) * log( (Na*Nd) / (nint^2) );
    Xn   = sqrt( (2*perm*Vint/q) * ( Na / (Nd*(Na+Nd)) ) );
    Xp   = Xn * Nd / Na;

    exp_term = expm1(q * V / (k * T)); % = exp(.) - 1

    JD = q * (nint^2) .* exp_term .* ...
        ( (Dn./(Ln*Na)) .* coth( (Wn + Xp)/Ln ) ...
        - (Dp./(Lp*Nd)) .* coth( (-Xn - Wp)/Lp ) );

    J = JL - JD;
    P = V .* J;
end

function [J, P] = JV_realista(V, Wn, Nd, Wp, Na, q,k,T,perm,Dn,tau_n,Dp,tau_p,JL,nint,Sf_eff,Sbsf)
    Ln   = sqrt(Dn * tau_n);
    Lp   = sqrt(Dp * tau_p);
    Vint = (k*T/q) * log( (Na*Nd) / (nint^2) );
    Xn   = sqrt( (2*perm*Vint/q) * ( Na / (Nd*(Na+Nd)) ) );
    Xp   = Xn * Nd / Na;

    we = max(Wn - Xn, 0) / Lp;   % huecos en región n (usa Lp)
    wb = max(Wp - Xp, 0) / Ln;   % electrones en región p (usa Ln)

    exp_term = expm1(q * V / (k * T)); % = exp(.) - 1

    alpha_e = (Sf_eff * Lp) / Dp;   % lado n
    alpha_b = (Sbsf   * Ln) / Dn;   % lado p

    den_e = cosh(we) + alpha_e .* sinh(we);
    den_b = cosh(wb) + alpha_b .* sinh(wb);
    tiny  = 1e-30;
    den_e = den_e + sign(den_e).*tiny;
    den_b = den_b + sign(den_b).*tiny;

    term1 = (Dp/(Lp*Nd)) * (sinh(we) + alpha_e .* cosh(we)) ./ den_e; 
    term2 = (Dn/(Ln*Na)) * (sinh(wb) + alpha_b .* cosh(wb)) ./ den_b;

    JD = q * (nint^2) .* exp_term .* (term1 + term2);
    J  = JL - JD;
    P  = V .* J;
end

function y = coth(x)
    y = cosh(x)./sinh(x);
end
