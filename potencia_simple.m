function P = potencia_simple(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint)
% x = [Wn, Nd, Wp, Na, V]
    Wn = x(1); Nd = x(2); Wp = x(3); Na = x(4); V = x(5);

    Ln = sqrt(Dn * tau_n);
    Lp = sqrt(Dp * tau_p);
    Vint = (k*T/q) * log( (Na*Nd) / (nint^2) );
    Xn   = sqrt( (2*perm*Vint/q) * ( Na / (Nd*(Na+Nd)) ) );
    Xp   = Xn * Nd / Na;

    % Exponencial robusta
    exp_term = expm1(q * V / (k * T)); % = exp(.) - 1
    
    % Corriente del diodo (modelo simple)
    JD = q * (nint^2) .* exp_term .* ...
        ( (Dn./(Ln*Na)) .* coth( (Wn + Xp)/Ln ) ...
        - (Dp./(Lp*Nd)) .* coth( (-Xn - Wp)/Lp ) );

    % Corriente neta y potencia
    J = JL - JD;
    P = V .* J;
end

% coth(x)
function y = coth(x)
    y = cosh(x)./sinh(x);
end

