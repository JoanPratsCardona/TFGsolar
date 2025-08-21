function H = hessi(f, x)
% hessi: Calcula la matriz Hessiana numérica de una función escalar
%
% Entradas:
%   f : función escalar (handle)
%   x : punto en el cual calcular la Hessiana (vector fila o columna)
%
% Salida:
%   H : matriz Hessiana (segunda derivada) de f evaluada en x

h = 1e-6;                    % Paso de derivación
ndim = length(x);            % Dimensión del vector x
H = zeros(ndim);             % Inicializar Hessiana

for i = 1:ndim
    for j = 1:ndim
        aa = zeros(1, ndim); aa(i) = 1;
        bb = zeros(1, ndim); bb(j) = 1;

        H(i, j) = ( ...
            f(x + h * (aa + bb)) ...
          - f(x + h * (-aa + bb)) ...
          - f(x + h * (aa - bb)) ...
          + f(x + h * (-aa - bb)) ) / (4 * h^2);
    end
end
end
