%%********************************************************************
%% The function HKLeq is for solving the following problem: 
%%
%% min{ 1/2 x^TDiag(d)x-a^T*x : sum(x) = c, 0<=x<=b } 
%% 
%% Based on Helgason, R., Kennington, J., Lall, H.: A polynomially bounded algorithm 
%% for a singly constrained quadratic program. Math. Program. 18, 338-343 (1980).
%% The codes are available at http://www.dingchao.info/codes/
%%*******************************************************************
function x = HKLeq(a,d,c,b)
n = length(a); 
x = zeros(size(a)); 

y = [a-d.*b;a];
y = sort(y);

L = sum(b); R = 0;
if c > L || c < 0
    disp(' ---No feasible solution exists--- ')
    return;
end
l = 1; r = 2*n; 

while r-l > 1
    m = floor(0.5*(l+r)); 
    temp = (a-y(m))./d;
    temp = max(min(temp,b),0);
    C = sum(temp);
    if C == c
        lambda = y(m);
        break;
    end
    if C > c
        l = m; L = C;
    end
    if C < c
        r = m; R = C;
    end
    if r-l == 1
       lambda = y(l)+(y(r)-y(l))*(c-L)/(R-L);
    end
end

t = a-b.*d; t1 = (a-lambda)./d;

x(t>=lambda) = b(t>=lambda);
x(t<lambda & a>=lambda) = t1(t<lambda & a>=lambda);
