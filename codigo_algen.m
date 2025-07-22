clear
clc
close all
format long

% -------------------------
% Constantes del modelo
% -------------------------
q = 1.602176634e-19;       
k = 1.3880649e-23;        
T = 300;                    
perm = 1.035918e-10;        
Dp = 11.6;                 
tau_p = 3710e-6;            
Dn = 2;                    
tau_n = 371e-6;             
JL = 50e-3;                
nint = 9.696e9;             
Sf_eff = 3e4;              
Sbsf = 100;                

% -------------------------
% Límites
% -------------------------
lb = [0.05e-4, 1e18, 100e-4, 1e14, 0.1];
ub = [0.75e-4, 5e20, 200e-4, 1e16, 0.52];
Nvars = length(lb);

% -------------------------
% Función objetivo
% -------------------------
objfun = @(x) -potencia_simple(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint);
%objfun = @(x) -potencia_realista(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint, Sf_eff, Sbsf);

% -------------------------
% Parámetros del AG
% -------------------------
Npop = 50;
Ngen = 300;
pmut = 0.1;
Nsic = 10;
txt = 1;
resfin = 1;

% -------------------------
% Ejecutar algoritmo genético
% -------------------------
[Beste, Bestf, conv, ~, ~] = algen(objfun, lb, ub, Npop, Ngen, pmut, Nsic, txt);

fprintf('\n🔍 Resultados del AG\n');
fprintf('Potencia máxima = %.4f mW/cm²\n', -Bestf*1e3);
fprintf('Wn = %.2e cm, Nd = %.2e cm⁻³\n', Beste(1), Beste(2));
fprintf('Wp = %.2e cm, Na = %.2e cm⁻³\n', Beste(3), Beste(4));
fprintf('V  = %.3f V\n', Beste(5));

% -------------------------
% Curva de convergencia
% -------------------------
figure;
plot(conv, 'LineWidth', 2);
grid on;
xlabel('Iteración');
ylabel('Función de costo');
title('Convergencia del AG');
