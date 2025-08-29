%% Cálculo de J(V) y P(V)
clc; clear; close all;

%Parámetros y constantes
q    = 1.602176634e-19;   % C
k    = 1.3880649e-23;     % J/K
T    = 300;               % K
perm = 1.035918e-10;      % F/cm
Wn   = 0.2e-4;            % cm
Nd   = 1e20;              % cm^-3
Dn   = 2;                 % cm^2/s
tau_n= 371e-6;            % s
Wp   = 179.8e-4;          % cm
Na   = 4.7e15;            % cm^-3
Dp   = 11.6;              % cm^2/s
tau_p= 3710e-6;           % s
JL   = 50e-3;             % A/cm²
nint = 9.696e9;           % cm^-3
Sf_eff = 3e4;             % cm/s (frontal, lado n)
Sbsf   = 100;             % cm/s (trasera, lado p)

%Magnitudes derivadas 
Vint = (k*T/q) * log(Na*Nd/nint^2);
Xn   = sqrt((2*perm*Vint/q) * (Na/(Nd*(Na+Nd))));  
Xp   = Xn * Nd/Na;                                  
Ln   = sqrt(Dn * tau_n);
Lp   = sqrt(Dp * tau_p);

we = max(Wn - Xn, 0) / Lp;   % huecos en n (usa Lp)
wb = max(Wp - Xp, 0) / Ln;   % electrones en p (usa Ln)

% Rango de voltaje 
Vmax = 0.7;                          
V    = linspace(0, Vmax, 2025);

% Exponencial robusta
exp_term = expm1(q * V / (k * T)); % = exp(.) - 1

Nj = (nint^2 / Na) .* exp_term;   % en el lado p (electrones)
Pj = (nint^2 / Nd) .* exp_term;   % en el lado n (huecos)

% Corriente simple
JD_simple = q * nint^2 .* exp_term .* ...
           ( (Dn/(Ln*Na)) * coth((Wn + Xp)/Ln) ...
           - (Dp/(Lp*Nd)) * coth((-Xn - Wp)/Lp) );
J  = JL - JD_simple;

% Corriente realista 
alfae = Sf_eff*Lp/Dp;
alfab = Sbsf  *Ln/Dn;
term1 = (Dp/(Lp*Nd)) * (sinh(we) + alfae*cosh(we)) / (cosh(we) + alfae*sinh(we));
term2 = (Dn/(Ln*Na)) * (sinh(wb) + alfab*cosh(wb)) / (cosh(wb) + alfab*sinh(wb));

JD_realista = q * nint^2 .* exp_term .* (term1 + term2);
Jr = JL - JD_realista;

% Potencia areal 
P  = V .* J;
Pr = V .* Jr;

% Valores de prueba V = 0.45 V 
V_test = 0.45;
exp_t  = expm1(q * V_test / (k * T));

JD_simple_t   = q * nint^2 * exp_t * ...
               ( (Dn/(Ln*Na)) * coth((Wn + Xp)/Ln) ...
               - (Dp/(Lp*Nd)) * coth((-Xn - Wp)/Lp) );
JD_realista_t = q * nint^2 * exp_t * (term1 + term2);

J_simple_t = JL - JD_simple_t;
J_realista_t = JL - JD_realista_t;
fprintf('V = 0.45 V -> J_simple = %.4e A/cm²,   J_realista = %.4e A/cm²\n', ...
        J_simple_t, J_realista_t);

% Gráficas (solo cuadrante positivo)
% máscaras para mantener solo Y>0 (y X ya es >=0)
maskJ  = J  > 0;  maskJr = Jr > 0;
maskP  = P  > 0;  maskPr = Pr > 0;

% J-V
figure('Name','J-V','Color','w');
plot(V(maskJ),  1e3*J(maskJ),  'b',  'LineWidth',2); hold on;
plot(V(maskJr), 1e3*Jr(maskJr), 'r--','LineWidth',2);
xlabel('Voltaje V (V)'); ylabel('Densidad de corriente (mA/cm²)');
title('Curvas J–V: simple vs. realista');
legend('Simple','Realista','Location','southwest'); grid on; box on;
xlim([0 Vmax]); ylim([0 inf]);  % solo positivos

% P-V
figure('Name','P-V','Color','w');
plot(V(maskP),  1e3*P(maskP),  'b',  'LineWidth',2); hold on;
plot(V(maskPr), 1e3*Pr(maskPr), 'r--','LineWidth',2);
xlabel('Voltaje V (V)'); ylabel('Potencia (mW/cm²)');
title('Curvas P–V: simple vs. realista');
legend('Simple','Realista','Location','southwest'); grid on; box on;
xlim([0 Vmax]); ylim([0 inf]);  % solo positivos

% Máximos de P y Pr
[~, idxP ] = max(P);   fprintf('P_simple max = %.4e W/cm² en V = %.4f V\n', P(idxP),  V(idxP));
[~, idxPr] = max(Pr);  fprintf('P_real   max = %.4e W/cm² en V = %.4f V\n', Pr(idxPr), V(idxPr));

% coth(x)
function y = coth(x)
    y = cosh(x)./sinh(x);
end

