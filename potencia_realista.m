function P = potencia_realista(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint, Sf_eff, Sbsf)
% x = [Wn, Nd, Wp, Na, V]
    Wn = x(1); Nd = x(2); Wp = x(3); Na = x(4); V = x(5);

    Ln = sqrt(Dn * tau_n);
    Lp = sqrt(Dp * tau_p);
    Vint = (k*T/q) * log( (Na*Nd) / (nint^2) );
    Xn   = sqrt( (2*perm*Vint/q) * ( Na / (Nd*(Na+Nd)) ) );
    Xp   = Xn * Nd / Na;

    % Longitudes adimensionales (no negativas)
    we = max(Wn - Xn, 0) / Lp;   % huecos en región n (usa Lp)
    wb = max(Wp - Xp, 0) / Ln;   % electrones en región p (usa Ln)

    % Exponencial robusta
    exp_term = expm1(q * V / (k * T)); % = exp(.) - 1

    alpha_e = (Sf_eff * Lp) / Dp;   % lado n
    alpha_b = (Sbsf   * Ln) / Dn;   % lado p
    den_e = cosh(we) + alpha_e .* sinh(we);   % lado n
    den_b = cosh(wb) + alpha_b .* sinh(wb);   % lado p
    tiny  = 1e-30; %garantizamos que el denominador no sea 0
    den_e = den_e + sign(den_e).*tiny;
    den_b = den_b + sign(den_b).*tiny;

    term1 = (Dp/(Lp*Nd)) * (sinh(we) + alpha_e .* cosh(we)) ./ den_e; 
    term2 = (Dn/(Ln*Na)) * (sinh(wb) + alpha_b .* cosh(wb)) ./ den_b;

    % Corriente neta y potencia
    JD = q * (nint^2) .* exp_term .* (term1 + term2);
    J = JL - JD;
    P = V .* J;
end

