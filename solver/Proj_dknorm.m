%%*************************************************************************
%% Compute the projection of a on the ball with center 0 and radius r 
%% in sense of the dual norm of k-norm:
%% x = Proj_knorm(a,c,b)
%% where b:=r*ones(length(a),1) and c:=k*r.
%% Based on Helgason, R., Kennington, J., Lall, H.: A polynomially bounded algorithm 
%% for a singly constrained quadratic program. Math. Program. 18, 338-343 (1980).
%%*************************************************************************
function x = Proj_dknorm(a,c,b)
a_abs=abs(a);
x=zeros(size(a));
if min(b)<=0 || c<0
    disp(' ---b and c should be non-negative--- ');
    return;
end

%% If a in B^r_{(k)*}, then x = a
if a_abs<=b
    if sum(a_abs)<=c
        x=a;
        return;
    end
end

%% If a is not in B^r_inf but its projection on B^r_inf is in B^kr_1 
xtemp=min(a_abs,b);
if sum(xtemp)< c
    x=sign(a).*xtemp;
    return;
end
x=HKLeq(a_abs,ones(size(a)),c,b);
x=sign(a).*x;
