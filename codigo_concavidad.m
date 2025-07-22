clc;clear;close all;

syms Wn Nd Wp Na V real
x = [Wn, Nd, Wp, Na, V];

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

Ln = sqrt(Dn * tau_n);
Lp = sqrt(Dp * tau_p);
Vint = (k*T/q)*log(Na*Nd/nint^2);
Xn = sqrt((2*perm*Vint/q)*(Na/(Nd*(Na+Nd))));
Xp = Xn * Nd/Na;

exp_term = exp(q * V / (k * T)) - 1;

JD = q * nint^2 * exp_term * ...
     ((Dn / (Ln * Na)) * coth((Wn + Xp) / Ln) - ...
      (Dp / (Lp * Nd)) * coth((-Xn - Wp) / Lp));
J = JL - JD;
P = J * V;
H = hessian(P, x);
xval = [0.2e-4, 1e19, 150e-4, 1e15, 0.4]; % punto en el dominio válido

% Sustituir variables
H_num = double(subs(H, {Wn, Nd, Wp, Na, V}, xval));

% Verificar definitud
eigenvals = eig(H_num);
disp('Autovalores de la Hessiana:')
disp(eigenvals)

if all(eigenvals >= 0)
    disp('La función es convexa en ese punto.')
elseif all(eigenvals <= 0)
    disp('La función es cóncava en ese punto.')
else
    disp('La función no es ni convexa ni cóncava en ese punto.')
end
