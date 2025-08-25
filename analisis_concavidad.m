%% Concavidad (Hessiana) y análisis + vistas 2×2 de P(Wn,Nd,Wp,Na,V) 
clear; clc; close all;

% Símbolos y parámetros 
syms Wn Nd Wp Na V 
x = [Wn, Nd, Wp, Na, V];

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

% Magnitudes derivadas
Ln   = sqrt(Dn * tau_n);
Lp   = sqrt(Dp * tau_p);
Vint = (k*T/q) * log(Na*Nd/nint^2);
Xn   = sqrt((2*perm*Vint/q) * (Na/(Nd*(Na+Nd))));
Xp   = Xn * Nd/Na;

% Exponencial 
exp_term = exp(q * V / (k * T)) - 1;

% Corriente simple
JD_simple = q * nint^2 .* exp_term .* ...
           ( (Dn/(Ln*Na)) * coth((Wn + Xp)/Ln) ...
           - (Dp/(Lp*Nd)) * coth((-Xn - Wp)/Lp) );
J  = JL - JD_simple;

% Potencial 
P = V * J;                     

% Hessiana simbólica de P respecto a [Wn, Nd, Wp, Na, V]
H = hessian(P, x);

% Límites y punto de evaluación
xmin = [0.05e-4, 1e18, 100e-4, 1e14, 0.10];   % [Wn, Nd, Wp, Na, V]
xmax = [0.75e-4, 5e20, 200e-4, 1e16, 0.52];
xval = [0.40e-4, 2.24e19, 150e-4, 1.00e15, 0.40];

% Evaluación numérica de la Hessiana
H_num = double(subs(H, {Wn, Nd, Wp, Na, V}, num2cell(xval)));

% Autovalores globales
format long e
eigenvals = eig(H_num);
disp('Autovalores de la Hessiana completa (5x5)');
disp(eigenvals.');
format short g

% Clasificación global en xval
if all(eigenvals >= 0)
    disp('=> La función es CONVEXA en ese punto (en R^5).');
elseif all(eigenvals <= 0)
    disp('=> La función es CÓNCAVA en ese punto (en R^5).');
else
    disp('=> La función es INDEFINIDA en ese punto (en R^5).');
end

% Corriente realista 
we    = (Wn - Xn)/Lp;
wb    = (Wp - Xp)/Ln;
alfae = Sf_eff*Lp/Dp;
alfab = Sbsf  *Ln/Dn;
term1 = (Dp/(Lp*Nd)) * (sinh(we) + alfae*cosh(we)) / (cosh(we) + alfae*sinh(we));
term2 = (Dn/(Ln*Na)) * (sinh(wb) + alfab*cosh(wb)) / (cosh(wb) + alfab*sinh(wb));
JD_realista = q * nint^2 .* exp_term .* (term1 + term2);
Jr = JL - JD_realista;

% Potencial realista
Pr = V * Jr;

% Hessiana simbólica realista
Hr = hessian(Pr, x);

% Evaluación numérica realista
Hr_num = double(subs(Hr, {Wn, Nd, Wp, Na, V}, num2cell(xval)));

% Autovalores realista
format long e
eigenvals_r = eig(Hr_num);
disp('Autovalores de la Hessiana (modelo REALISTA, 5x5)');
disp(eigenvals_r.');
format short g

% Clasificación global en xval (realista)
if all(eigenvals_r >= 0)
    disp('=> [REALISTA] La función es CONVEXA en ese punto (en R^5).');
elseif all(eigenvals_r <= 0)
    disp('=> [REALISTA] La función es CÓNCAVA en ese punto (en R^5).');
else
    disp('=> [REALISTA] La función es INDEFINIDA en ese punto (en R^5).');
end
%% Análisis 2×2 (todas las parejas)  -  Modelo simple
varNames = {'Wn','Nd','Wp','Na','V'}; 
pairs = nchoosek(1:5, 2); % 15 parejas 
fprintf('Análisis 2x2 por parejas (bloques de la Hessiana en xval)\n');
for r = 1:size(pairs,1) 
    idx = pairs(r,:); 
    H2 = H_num(idx, idx); 
    ev2 = eig(H2); a = ev2(1);
    b = ev2(2); 
    if all(ev2 >= 0) 
        cls = 'convexa en ese plano'; 
    elseif all(ev2 <= 0) 
        cls = 'cóncava en ese plano';
    else 
        cls = 'indefinida en ese plano'; 
    end 
    fprintf('(%s, %s): ev = [%.3e, %.3e] => %s\n', ... 
        varNames{idx(1)}, varNames{idx(2)}, a, b, cls); 
end

%% Análisis 2×2 (todas las parejas)  -  Modelo realista
fprintf('Análisis 2x2 por parejas (bloques de la Hessiana REALISTA en xval)\n');
for r = 1:size(pairs,1)
    idx = pairs(r,:);
    H2r = Hr_num(idx, idx);
    ev2r = eig(H2r); ar = ev2r(1);
    br = ev2r(2);
    if all(ev2r >= 0)
        clsr = 'convexa en ese plano';
    elseif all(ev2r <= 0)
        clsr = 'cóncava en ese plano';
    else
        clsr = 'indefinida en ese plano';
    end
    fprintf('REALISTA (%s, %s): ev = [%.3e, %.3e] => %s\n', ...
        varNames{idx(1)}, varNames{idx(2)}, ar, br, clsr);
end

%% Gráficas 2x2 (todas las parejas) sobre [xmin, xmax]  -  Modelo simple
graficas = true; 
if graficas
Pfun = matlabFunction(P, 'Vars', {Wn, Nd, Wp, Na, V});
P0 = Pfun(xval(1), xval(2), xval(3), xval(4), xval(5));
varNames = {'Wn','Nd','Wp','Na','V'};
pairs    = nchoosek(1:5, 2);
npts = 50;   
for pp = 1:size(pairs,1)
    i = pairs(pp,1);
    j = pairs(pp,2);
    Xi = linspace(xmin(i), xmax(i), npts);
    Xj = linspace(xmin(j), xmax(j), npts);
    [XI, XJ] = meshgrid(Xi, Xj);
    Pgrid = zeros(size(XI));
    for a = 1:npts
        for b = 1:npts
            xv = xval;
            xv(i) = XI(a,b);
            xv(j) = XJ(a,b);
            Pgrid(a,b) = Pfun(xv(1), xv(2), xv(3), xv(4), xv(5));
        end
    end
    figure('Color','w','Name',sprintf('P vs %s y %s', varNames{i}, varNames{j}));
    surf(XI, XJ, Pgrid, 'EdgeAlpha', 0.2);
    xlabel(varNames{i}); ylabel(varNames{j}); zlabel('P (W/cm^2)');
    title(sprintf('Potencia (modelo simple): %s-%s en [xmin, xmax]', varNames{i}, varNames{j}));
    grid on; box on; view(35,30); shading interp;
    hold on; plot3(xval(i), xval(j), P0, 'ro', 'MarkerFaceColor','r');
end
end

%% Gráficas 2x2 (todas las parejas) sobre [xmin, xmax]  -  Modelo realista
graficas_realista = true;
if graficas_realista
Prfun = matlabFunction(Pr, 'Vars', {Wn, Nd, Wp, Na, V});
P0r = Prfun(xval(1), xval(2), xval(3), xval(4), xval(5));
varNames = {'Wn','Nd','Wp','Na','V'};
pairs    = nchoosek(1:5, 2);
npts = 50;
for pp = 1:size(pairs,1)
    i = pairs(pp,1);
    j = pairs(pp,2);
    Xi = linspace(xmin(i), xmax(i), npts);
    Xj = linspace(xmin(j), xmax(j), npts);
    [XI, XJ] = meshgrid(Xi, Xj);
    Pgrid = zeros(size(XI));
    for a = 1:npts
        for b = 1:npts
            xv = xval;
            xv(i) = XI(a,b);
            xv(j) = XJ(a,b);
            Pgrid(a,b) = Prfun(xv(1), xv(2), xv(3), xv(4), xv(5));
        end
    end
    figure('Color','w','Name',sprintf('Pr vs %s y %s', varNames{i}, varNames{j}));
    surf(XI, XJ, Pgrid, 'EdgeAlpha', 0.2);
    xlabel(varNames{i}); ylabel(varNames{j}); zlabel('P_r (W/cm^2)');
    title(sprintf('Potencia (modelo realista): %s-%s en [xmin, xmax]', varNames{i}, varNames{j}));
    grid on; box on; view(35,30); shading interp;
    hold on; plot3(xval(i), xval(j), P0r, 'ro', 'MarkerFaceColor','r');
end
end
