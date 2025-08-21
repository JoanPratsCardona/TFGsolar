function J = grad(f, x, o)
% grad: Calcula el gradiente numérico de una función escalar
%
% Entradas:
%   f : función a derivar (handle)
%   x : punto en el cual derivar (vector fila o columna)
%   o : esquema de derivación
%       o = 1  → derivada hacia adelante (f(x+h) - f(x)) / h
%       o ≠ 1 → derivada centrada   (f(x+h) - f(x-h)) / 2h
%
% Salida:
%   J : vector gradiente de f evaluado en x

h = 1e-6;                     % Paso de derivación
ndim = length(x);             % Dimensión del vector x
J = zeros(1, ndim);           % Inicializar gradiente

for i = 1:ndim
    aa = zeros(1, ndim);
    aa(i) = 1;
    if o == 1
        J(i) = (f(x + h * aa) - f(x)) / h;
    else
        J(i) = (f(x + h * aa) - f(x - h * aa)) / (2 * h);
    end
end
end
