%% Optimización local
clear; clc; close all

% Constantes
q    = 1.602176634e-19;   % C
k    = 1.3880649e-23;     % J/K
T    = 300;               % K
perm = 1.035918e-10;      % F/cm
Dn   = 2;                 % cm^2/s
tau_n= 371e-6;            % s
Dp   = 11.6;              % cm^2/s
tau_p= 3710e-6;           % s
JL   = 50e-3;             % A/cm²
nint = 9.696e9;           % cm^-3
Sf_eff = 3e4;             % cm/s (frontal, lado n)
Sbsf   = 100;             % cm/s (trasera, lado p)

% Límites de optimización
xmin = [0.05e-4, 1e18, 100e-4, 1e14, 0.1];   % [Wn, Nd, Wp, Na, V]
xmax = [0.75e-4, 5e20, 200e-4, 1e16, 0.52];
x0   = (xmin + xmax)/2;                      % Punto inicial


for i = 1:2
    rng(1,'twister'); % fijamos semilla y algoritmo (reproducible)
    % Selección de función objetivo 
    if i == 1
        nombre_fun = 'potencia_realista';
        objfun = @(x) -potencia_realista(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint, Sf_eff, Sbsf);
    elseif i == 2
        nombre_fun = 'potencia_simple';
        objfun = @(x) -potencia_simple(  x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint);
    end

    for alg = 1:4
        % Selección de algoritmo fmincon y Descenso del gradiente
        if     alg == 1
            algoritmo = 'interior-point';
            options   = optimoptions('fmincon','Algorithm','interior-point','Display','off');
        elseif alg == 2
            algoritmo = 'sqp';
            options   = optimoptions('fmincon','Algorithm','sqp','Display','off');
        elseif alg == 3
            algoritmo = 'active-set';
            options   = optimoptions('fmincon','Algorithm','active-set','Display','off');
        else
            algoritmo = 'Descenso del gradiente';
        end

        % Optimización
        if alg == 4
            txt = 0; %txt = 1, mostramos la información de cada iteración del método, cambiamos a cualquier otro valor para no mostrarla
            [Jsol,xsol] = gradmet(objfun, x0 , 50, 1, -100, 1e-8, 1e-8, xmin, xmax, txt, 1); 
        else
        [xsol, Jsol] = fmincon(objfun, x0, [], [], [], [], xmin, xmax, [], options);
        end
        % Resultados 
        fprintf('\nFunción: %s | Algoritmo: %s\n', nombre_fun, algoritmo);
        fprintf('Potencia máxima: %.6e W/cm²\n', -Jsol);
        fprintf('Wn = %.3e cm | Nd = %.3e cm^-3 | Wp = %.3e cm | Na = %.3e cm^-3 | V = %.3f V\n', ...
                 xsol(1), xsol(2), xsol(3), xsol(4), xsol(5));
    end
end

