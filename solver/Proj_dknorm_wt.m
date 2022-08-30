%%********************************************************************
%% Compute the projection of one vector a under weighted dual norm 
%% of k-norm with the center 0 and radius r.
%% x = Proj_knorm_wt(a,c,b,eta)
%% where b:=r*ones(length(a),1), c:=k*r, eta is a given weighted vector of the vector a.
%% Based on Helgason, R., Kennington, J., Lall, H.: A polynomially bounded algorithm 
%% for a singly constrained quadratic program. Math. Program. 18, 338-343 (1980).
%%******************************************************************

function x = Proj_dknorm_wt(a,c,b,eta)

a_abs = abs(a)./eta;
x = zeros(size(a));
if min(b)<=0 || c<0
    disp(' ---b and c should be non-negative--- ');
    return;
end

%% If a in B^r_{(k)*}, then x = a
if a_abs <= b
    if sum(a_abs) <= c
        x = a;
        return;
    end
end

%% If a is not in B^r_inf but its projection on B^r_inf is in B^kr_1 
xtemp = min(a_abs,b);
if sum(xtemp) <  c
    x = sign(a).*eta.*xtemp;
    return;
end
eta_sqr = eta.^2;

x = HKLeq(eta_sqr.*a_abs,eta_sqr,c,b);

x = sign(a).*eta.*x;
