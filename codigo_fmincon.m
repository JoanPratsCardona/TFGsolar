clear; clc; close all

%% Constantes
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

%% Límites de optimización
xmin = [0.05e-4, 1e18, 100e-4, 1e14, 0.1];   % [Wn, Nd, Wp, Na, V]
xmax = [0.75e-4, 5e20, 200e-4, 1e16, 0.52];
x0 = (xmin + xmax)/2;   % Punto inicial

%% Función objetivo (negativa para maximizar)
objfun = @(x) -potencia_realista(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint, Sf_eff, Sbsf);
% objfun = @(x) -potencia_simple(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint);

%% Ejecución de fmincon
[xsol, Jsol] = fmincon(objfun, x0, [], [], [], [], xmin, xmax);
%% Mostrar resultados
fprintf('Potencia máxima: %.6e W/cm²\n', -Jsol);
fprintf('Wn  = %.3e cm\n', xsol(1));
fprintf('Nd  = %.3e cm⁻³\n', xsol(2));
fprintf('Wp  = %.3e cm\n', xsol(3));
fprintf('Na  = %.3e cm⁻³\n', xsol(4));
fprintf('V   = %.3f V\n', xsol(5));
