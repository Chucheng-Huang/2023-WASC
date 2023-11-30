function [P,W,R,Y,loss,lambda1] = WASC(X, B, r, mu, lambda, InitY, NITER)
%   Objective Funtion:
%   max_{P,W,R,Y}  tr(P'*X'*B'*W*W'*B'*X*P) + tr(R'*P'*X'*Y) - lambda*||W||_F^2
%   s.t.  P^TS_tP=I_r,R'*R = I_c, W = diag(w), w*1^T = 1,w>=0.

    [n,d] = size(X);
    m = size(B,2);

    % Init W
    w = 1/m*ones(1,m);
    W = diag(w);

    % Init Y
    Y = InitY;
    c = size(Y,2);

    % Init R
    R = rand(r,c);
    [u,~,v] = svd(R);
    if r >= c
        R = u(:,1:c)*v';
    else
        R = u*(v(:,1:r))';
    end

    % Init P
    St = X'*X+10^-7*eye(d);
    AA = St^(-1/2);

    loss = [];
    for iter = 1:NITER

        % Update P
        M = X'*B*W*W'*B'*X;
        T1 = AA'*M*AA;
        N = X'*Y*R';
        T2 = mu*AA'*N;
        [Q,~] = GPI_WASC(T1, T2, 60);
        P = AA*Q;

        % Udpate W
        PM = (B'*X*P*P'*X'*B).*eye(m);
        PMmax = max(diag(PM));
        lambda1 = PMmax;
        PP = zeros(1,m);
        [w, ~] = EProjSimplex_diag_yb(-(PM - lambda*eye(m)),PP');
        W = diag(w);

        % Update R
        tmp = P'*X'*Y;
        [u,~,v] = svd(tmp);
        if r >= c
            R = u(:,1:c)*v';
        else
            R = u*(v(:,1:r))';
        end


        % Update Y
        D = R'*P'*X';
        [~,Label] = max(D);
        Y = zeros(n,c);
        for i =1:1:length(Label)
            ind = Label(i);
            Y(i,ind) = 1;
        end
        
        loss(iter) = trace(P'*X'*B*W*W'*B'*X*P) + mu * trace(R'*P'*X'*Y) - lambda * (norm(W, 'fro'))^2;
       

        if iter > 10 && abs(loss(iter) - loss(iter - 5)) / abs(loss(iter)) < 10^-6
            return 
        end   
    end
end





