%% Optimización global (Algoritmo genético) + refinamiento (SQP y descenso del gradiente)
clear; clc; close all

% Constantes
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

% Límites de optimización
xmin = [0.05e-4, 1e18, 100e-4, 1e14, 0.10];   % [Wn, Nd, Wp, Na, V]
xmax = [0.75e-4, 5e20, 200e-4, 1e16, 0.52];

% Comprobamos disponibilidad de herramientas
tieneGA = ~isempty(which('ga'));   
tieneFmincon = ~isempty(which('fmincon'));
tieneAlgen   = exist('algen','file')==2;
tieneGradmet = exist('gradmet','file')==2;   

% Bucle sobre las dos funciones objetivo
for i = 1:2
    % Selección de función objetivo
    if i == 1
        nombre_fun = 'potencia_realista';
        objfun = @(x) -potencia_realista(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint, Sf_eff, Sbsf);
    else
        nombre_fun = 'potencia_simple';
        objfun = @(x) -potencia_simple(  x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint);
    end

    fprintf('\n==============================\n');
    fprintf('>> Optimizando (global + local): %s\n', nombre_fun);

    % Algoritmo genético (global) 
    if tieneGA
        rng(1,'twister'); % fijamos semilla y algoritmo (reproducible)
        optsGA = optimoptions('ga', ...
            'Display','off', ...
            'PopulationSize', 60, ...
            'MaxGenerations', 200, ...
            'FunctionTolerance', 1e-8, ...
            'UseParallel', false);

        [x_ga, J_ga] = ga(objfun, 5, [], [], [], [], xmin, xmax, [], [], optsGA);
        algoritmo = 'Algoritmo genético (ga)';
    elseif tieneAlgen
        rng(1,'twister');
        Npop = 60; Ngen = 200; pmut = 0.10; Nsic = 10; txt = 0;
        [x_ga, J_ga] = algen(objfun, xmin, xmax, Npop, Ngen, pmut, Nsic, txt);
        algoritmo = 'Algoritmo genético (propio)';
    else
        error(['No está disponible el Global Optimization Toolbox (ga) ni se encontró algen.m.\n' ...
               'Instala el toolbox o añade algen.m al path.']);
    end

    fprintf('\nFunción: %s | Algoritmo: %s\n', nombre_fun, algoritmo);
    fprintf('Potencia máxima (global): %.6e W/cm²\n', -J_ga);
    fprintf('Wn = %.3e cm | Nd = %.3e cm^-3 | Wp = %.3e cm | Na = %.3e cm^-3 | V = %.3f V\n', ...
             x_ga(1), x_ga(2), x_ga(3), x_ga(4), x_ga(5));

    % Refinamiento local: fmincon (SQP) 
    if tieneFmincon
        optsFC = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
            'MaxIterations', 500, 'MaxFunctionEvaluations', 2e4, ...
            'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);

        [x_sqp, J_sqp] = fmincon(objfun, x_ga, [], [], [], [], xmin, xmax, [], optsFC);

        fprintf('\nFunción: %s | Algoritmo: Refinamiento local (SQP)\n', nombre_fun);
        fprintf('Potencia máxima (ref. SQP): %.6e W/cm²\n', -J_sqp);
        fprintf('Wn = %.3e cm | Nd = %.3e cm^-3 | Wp = %.3e cm | Na = %.3e cm^-3 | V = %.3f V\n', ...
                 x_sqp(1), x_sqp(2), x_sqp(3), x_sqp(4), x_sqp(5));
    else
        fprintf('\n(Refinamiento local con fmincon omitido: Optimization Toolbox no disponible)\n');
        x_sqp = x_ga; J_sqp = J_ga; % para usar como arranque del gradiente si se desea
    end

    % Refinamiento local: descenso por gradiente (gradmet) 
    if tieneGradmet
        % Arrancamos desde la mejor solución conocida (SQP si existe, si no GA)
        x_start = x_sqp;
        txt = 0;  % txt = 1 para ver el progreso por iteración
        [J_gd, x_gd] = gradmet(objfun, x_start, 50, 1, -100, 1e-8, 1e-8, xmin, xmax, txt, 1);

        fprintf('\nFunción: %s | Algoritmo: Refinamiento local (Descenso del gradiente)\n', nombre_fun);
        fprintf('Potencia máxima (ref. GD): %.6e W/cm²\n', -J_gd);
        fprintf('Wn = %.3e cm | Nd = %.3e cm^-3 | Wp = %.3e cm | Na = %.3e cm^-3 | V = %.3f V\n', ...
                 x_gd(1), x_gd(2), x_gd(3), x_gd(4), x_gd(5));
    else
        fprintf('\n(Refinamiento por gradiente omitido: no se encontró gradmet.m)\n');
    end
end
