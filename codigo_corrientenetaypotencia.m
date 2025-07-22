clc;clear;close all;

%Parámetros y constantes
q = 1.602176634e-19;        % Carga elemental (C) 
k = 1.3880649e-23;          % Constante de Boltzmann (J/K)
T = 300;                    % Temperatura (K)
perm = 1.035918e-10;        % Permitividad del semiconductor (F/cm)
Wn = 0.2e-4;                % Ancho región n (cm)
Nd = 1e20;                  % Dopaje n (cm^-3)
Dp = 11.6;                  % Difusión huecos (cm^2/s)
tau_p = 3710e-6;            % Vida huecos (s)
Wp = 179.8e-4;              % Ancho región p (cm)
Na = 4.7e15;                % Dopaje p (cm^-3)
Dn = 2;                     % Difusión electrones (cm^2/s)
tau_n = 371e-6;             % Vida electrones (s)
JL = 50e-3;                 % Corriente fotogenerada (A/cm²)
nint = 9.696e9;             % Concentración intrínseca (cm^-3)
Sf_eff = 3e4;               % velocidad efectiva de recombinación superficial en la región BSF (cm/s)
Sbsf = 100;                 % velocidad efectiva de recombinación superficial en la superficie frontal (cm/s)

% Parámetros auxiliares
Vint = (k*T/q)*log(Na*Nd/nint^2); 
Xn = sqrt((2*perm*Vint/q)*(Na/(Nd*(Na+Nd))));
Xp = Xn * Nd/Na;
Ln = sqrt(Dn * tau_n);
Lp = sqrt(Dp * tau_p);

% Rango de voltaje
V = linspace(0, 0.55, 500);
J = zeros(size(V));

%CONDICIONES DE FROTERA SIMPLES

% Cálculo de la corriente neta, J 
for i = 1:length(V)
    exp_term = exp(q * V(i) / (k * T)) - 1;
    JD = q * nint^2 * exp_term * ...
         ((Dn / (Ln * Na)) * coth((Wn + Xp) / Ln) - ...
          (Dp / (Lp * Nd)) * coth((-Xn - Wp) / Lp));
    J(i) = JL - JD;  
end

% Potencia
P = J .* V;  % W/cm²

% Calcular y mostrar J para V = 0.45
V_test = 0.45;
exp_test = exp(q * V_test / (k * T)) - 1;
JD_test = q * nint^2 * exp_test * ...
          ((Dn / (Ln * Na)) * coth((Wn + Xp) / Ln) - ...
           (Dp / (Lp * Nd)) * coth((-Xn - Wp) / Lp));
J_test = JL - JD_test;
fprintf('Para V = 0.45 V, J = %.4e A/cm²\n', J_test);

%CONDICIONES DE FRONTERA MÁS REALISTAS

% Cálculo de la corriente neta, Jr

for i = 1:length(V)
    exp_term = exp(q * V(i) / (k * T)) - 1;
    % Jnr
    num_n = (nint^2 / Na) * exp_term;
    Bnr_factor = ((exp(Wp/Ln)/Ln) + (Sbsf/Dn)*exp(Wp/Ln)) / ...
                 ((exp(-Wp/Ln)/Ln) - (Sbsf/Dn)*exp(-Wp/Ln));
    denom_n = exp(Xp/Ln) + Bnr_factor * exp(-Xp/Ln);
    Anr = num_n / denom_n;
    Bnr = Bnr_factor * Anr;

    Jnr = q * Dn * ((Anr/Ln)*exp(Xp/Ln) - (Bnr/Ln)*exp(-Xp/Ln));
    % Jpr
    num_p = (nint^2 / Nd) * exp_term;
    Bpr_factor = ((exp(-Wn/Lp)/Lp) - (Sf_eff/Dp)*exp(-Wn/Lp)) / ...
                 ((exp(Wn/Lp)/Lp) + (Sf_eff/Dp)*exp(Wn/Lp));
    denom_p = exp(-Xn/Lp) + Bpr_factor * exp(Xn/Lp);
    Apr = num_p / denom_p;
    Bpr = Bpr_factor * Apr;

    Jpr = -q * Dp * ((Apr/Lp)*exp(-Xn/Lp) - (Bpr/Lp)*exp(Xn/Lp));

    % Corriente de recombinación total
    Jdr = Jpr + Jnr;
    Jr(i) = JL - Jdr;
end

% Potencia
Pr = Jr .* V;

% Cálculo de Jr para V = 0.45 V
V_test = 0.45;
exp_test = exp(q * V_test / (k * T)) - 1;

num_n = (nint^2 / Na) * exp_test;
Bnr_factor = ((exp(Wp/Ln)/Ln) + (Sbsf/Dn)*exp(Wp/Ln)) / ...
             ((exp(-Wp/Ln)/Ln) - (Sbsf/Dn)*exp(-Wp/Ln));
denom_n = exp(Xp/Ln) + Bnr_factor * exp(-Xp/Ln);
Anr = num_n / denom_n;
Bnr = Bnr_factor * Anr;
Jnr_test = q * Dn * ((Anr/Ln)*exp(Xp/Ln) - (Bnr/Ln)*exp(-Xp/Ln));

num_p = (nint^2 / Nd) * exp_test;
Bpr_factor = ((exp(-Wn/Lp)/Lp) - (Sf_eff/Dp)*exp(-Wn/Lp)) / ...
             ((exp(Wn/Lp)/Lp) + (Sf_eff/Dp)*exp(Wn/Lp));
denom_p = exp(-Xn/Lp) + Bpr_factor * exp(Xn/Lp);
Apr = num_p / denom_p;
Bpr = Bpr_factor * Apr;
Jpr_test = -q * Dp * ((Apr/Lp)*exp(-Xn/Lp) - (Bpr/Lp)*exp(Xn/Lp));

JD_test = Jpr_test + Jnr_test;
Jr_test = JL - JD_test;
fprintf('Para V = 0.45 V, Jr = %.4e A/cm²\n', Jr_test);

%Gráfica de la curva característica de la célula solar
figure;
plot(V, J*1e3, 'b', 'LineWidth', 2); hold on;
plot(V, Jr*1e3, 'r--', 'LineWidth', 2);
xlabel('Voltaje V (V)');
ylabel('Densidad de corriente (mA/cm²)');
title('Curvas J-V vs Jr-V');
legend('Condición simple (J)', 'Condición realista (Jr)', 'Location', 'southwest');
grid on;
axis([0 0.55 0 inf]);

%Gráfica de P contra el V
figure;
plot(V, P*1e3, 'b', 'LineWidth', 2); hold on;
plot(V, Pr*1e3, 'r--', 'LineWidth', 2);
xlabel('Voltaje V (V)');
ylabel('Potencia (mW/cm²)');
title('Curvas P-V vs Pr-V');
legend('Condición simple (P)', 'Condición realista (Pr)', 'Location', 'southwest');
grid on;
axis([0 0.55 0 inf]);

