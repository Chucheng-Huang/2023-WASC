% max Tr(P'AP)+Tr(P'B)
% s.t. P'P=I_r
% P   d*r
% A   d*d
% B   d*r

function [P,obj]=GPI_WLPP(A, B, Niter)
    [d,r] = size(B);
    tmp = abs(rand(d,r));
    [u,~,v]= svd(tmp);
    if d >= r
        P = u(:,1:r)*v';
    else
        P = u*(v(:,1:d))';
    end

%     I = eye(d);
%     AA = A + eps * I;
    obj =[];
    for i =1: Niter
        M = 2*A*P + B;
        [U,sigma,V] = svd(M);
        if d >= r
            P = U(:,1:r)*V';
        else
            P = U*(V(:,1:d))';
        end
        obj(i) = trace(P'*A*P + P'*B);

        if i >= 6 
            if   abs(obj(i) - obj(i-5))/ abs(obj(i-5)) < 1e-6 
                break;
            end
        end 
    end
end
