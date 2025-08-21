function [Jsol, xsol, Jit, xit, i, timf] = gradmet(f, x0, Im, Tm, e1, e2, e3, xmin, xmax, txt, met)
% Algoritmo de descenso: Gradiente / Newton
% Entradas:
%   f     = función objetivo (handle)
%   x0    = punto inicial
%   Im    = número máximo de iteraciones
%   Tm    = tiempo máximo (segundos)
%   e1,e2,e3 = criterios de parada
%   xmin, xmax = límites de las variables
%   txt   = 1 para mostrar estado, 0 para silencioso
%   met   = 1: Gradiente Descendente, 2: Newton

% Salidas:
%   Jsol  = mínimo encontrado de f
%   xsol  = punto que minimiza f
%   Jit   = evolución del valor de f
%   xit   = evolución del punto x
%   i     = número de iteraciones realizadas
%   timf  = tiempo de ejecución

if txt == 1
    if met == 2
        disp('Algoritmo de Newton');
    else
        disp('Algoritmo de Descenso por Gradiente');
    end
end

tic;
xit(1,:) = x0;
h(1) = 0.1;
Jit(1) = f(x0);

for i = 1:Im
    if met == 2
        coefg = inv(hessi(f, xit(i,:)));  
        gdit(i,:) = grad(f, xit(i,:), 1); 
    else
        gdit(i,:) = grad(f, xit(i,:), 1); % gradiente numérico
        coefg = h(i);
    end

    % Paso del algoritmo
    delta = (coefg * gdit(i,:)')';
    xnext = xit(i,:) - delta;

    % Aplicar límites
    xnext = max(min(xnext, xmax), xmin);
    xit(i+1,:) = xnext;
    Jit(i+1) = f(xnext);

    % Ajuste de h si es gradiente
    if met == 1
        if Jit(i+1) >= Jit(i)
            h(i+1) = h(i) / 2;
        else
            h(i+1) = h(i);
        end
    end

    % Mostrar información
    if txt == 1
        disp(['Iteración: ' num2str(i)]);
        disp(['Punto actual: ' num2str(xit(i+1,:))]);
        disp(['Valor de la función: ' num2str(Jit(i+1))]);
        disp(['Norma del gradiente: ' num2str(norm(gdit(i,:),2))]);
        if met == 1
            disp(['Valor de h: ' num2str(h(i))]);
        end
        disp(['Tiempo: ' num2str(toc) ' s']);
        disp(' ');
    end

    % Criterios de parada
    if Jit(i+1) < e1
        disp('Parada: J < e1');
        break;
    end
    if norm(gdit(i,:),2) < e2
        disp('Parada: ||gradJ|| < e2');
        break;
    end
    if abs(Jit(i+1) - Jit(i)) < e3
        disp('Parada: |J(i+1) - J(i)| < e3');
        break;
    end
    if toc > Tm
        disp('Parada: tiempo > Tm');
        break;
    end
end

if i == Im
    disp('Parada: se alcanzó el número máximo de iteraciones');
end

[Jsol, ind] = min(Jit);
xsol = xit(ind,:);
timf = toc;

end
