function P = potencia_realista(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint, Sf_eff, Sbsf)
    Wn = x(1); Nd = x(2); Wp = x(3); Na = x(4); V = x(5);

    Ln = sqrt(Dn * tau_n);
    Lp = sqrt(Dp * tau_p);
    Vint = (k*T/q)*log(Na*Nd/nint^2);
    Xn = sqrt((2*perm*Vint/q)*(Na/(Nd*(Na+Nd))));
    Xp = Xn * Nd/Na;

    exp_term = exp(q * V / (k * T)) - 1;

    % Jnr
    num_n = (nint^2 / Na) * exp_term;
    Bnr_factor = ((exp(Wp/Ln)/Ln) + (Sbsf/Dn)*exp(Wp/Ln)) / ...
                 ((exp(-Wp/Ln)/Ln) - (Sbsf/Dn)*exp(-Wp/Ln));
    denom_n = exp(Xp/Ln) + Bnr_factor * exp(-Xp/Ln);
    Anr = num_n / denom_n;
    Bnr = Bnr_factor * Anr;
    Jnr = q * Dn * ((Anr/Ln)*exp(Xp/Ln) - (Bnr/Ln)*exp(-Xp/Ln));

    % Jpr
    num_p = (nint^2 / Nd) * exp_term;
    Bpr_factor = ((exp(-Wn/Lp)/Lp) - (Sf_eff/Dp)*exp(-Wn/Lp)) / ...
                 ((exp(Wn/Lp)/Lp) + (Sf_eff/Dp)*exp(Wn/Lp));
    denom_p = exp(-Xn/Lp) + Bpr_factor * exp(Xn/Lp);
    Apr = num_p / denom_p;
    Bpr = Bpr_factor * Apr;
    Jpr = -q * Dp * ((Apr/Lp)*exp(-Xn/Lp) - (Bpr/Lp)*exp(Xn/Lp));

    Jr = JL - (Jnr + Jpr);
    P = Jr * V;
end
