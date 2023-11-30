function [x obj] = EProjSimplex_diag_yb(D, p)
%   min_u  1/2 u' * D * u  - u' * p
%   s.t. u'*1 =1, uj >= 0. 


% clc;clear;
% U = [1  0 0; 0 1 0 ; 0 0 1]
% p = [1 1 1.1]'; 
% EProjSimplex_new(p');

u = eps + diag(D); u1 = 1./ u ;
eta = max(abs(p)) + 1; 

tmp = eta * u1 + p.*u1;
ind = find(tmp > 0);
obj = sum(tmp(ind)) - 1 ;   % º¯ÊýÖµ
g_ = sum(u1(ind));          % gradient 

iter = 1; NITER = 200;
while abs( obj ) >  10^ - 10 && iter < NITER
    iter = iter + 1;
    eta = eta - obj / g_;
    
    %%% update obj
    tmp = eta * u1 + p.*u1;
    ind = find(tmp > 0);
    obj = sum(tmp(ind)) - 1 ;
    
    %%% update g_
    g_ = sum(u1(ind));   
end

    s = eta * u1 + p.*u1 ;
    x = max(s,0);
end



% [x1 ~] = EProjSimplex_new(p');
% EProjSimplex_new(p')
