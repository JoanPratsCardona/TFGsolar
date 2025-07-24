function P = potencia_simple(x, q, k, T, perm, Dn, tau_n, Dp, tau_p, JL, nint)
    Wn = x(1); Nd = x(2); Wp = x(3); Na = x(4); V = x(5);
    
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
end
