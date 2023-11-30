clc; clear;

%% Data
load Yeast\data.txt;
load Yeast\gnd.txt; 
X = data;
% load texture.mat;
% gnd = Y;

%% Parameter
[n,d] = size(X);
c = length(unique(gnd));
m = 0.2 * n;    % anchors
k = 5;          % neighbors of anchor graph

if d < 20
    rL  = [2:1:d];   %dimentionality
elseif d < 100
    rL = [5:5:d];
else
  % 当维度大于100时，用PCA预处理加速
    [eigvector, eigvalue, elapse] = PCA_dencai(X, 100); 
    projection=eigvector;%(:,1:dim);
    X = X* projection;
    rL = [10:10:d];
end
muL = [10^-6,10^-5,10^-4,10^-3,10^-2,10^-1];
lambdaL = [0.2,0.4,0.6,0.8];
niter = 50;

%% Anchor Generation
tic
RandIDX = randperm(n,floor(m));
locAnchor = X(RandIDX,:);
B0 = ConstructA_NP(X', locAnchor', k);
[B] = norm_AnchorGraphB(B0);
t1 = toc;      

%% InitY
G  = rand(c,n);
[~,InitLabel] = max(G);
InitY = zeros(n,c);
for i =1:1:length(InitLabel)
    ind = InitLabel(i);
    InitY(i,ind) = 1;
end

%% WASC
optAcc = 0;
step = 0;
tag = 1;
Record  = [];
for ri = 1:1:length(rL)
    r = rL(ri);
    for mui = 1:1:length(muL)
        mu = muL(mui);
        for i = 1:1:length(lambdaL)
            lambda = lambdaL(i);
            step = step + 1;
            tic;
            [P, W, R, Y_WASC, loss,lambda1] = WASC(X, B, r, mu, lambda, InitY, niter); 
            [~,Label_WLPP] = max(Y_WASC');
            t2 = toc;
            
            tempACC = ACC2(gnd, Label_WLPP', c);
            tempNMI = NMI(gnd, Label_WLPP');
    
            tmp = [step,r,mu,lambda,lambda1,t2,tempNMI,tempACC];
            Record = [Record;tmp];
            if tempACC > optAcc
                Paras{tag} = {tmp,P,W,R,loss,mu,lambda,lambda1};
                optAcc = tempACC; 
            end
            fprintf('Step (anchor_WASC): id=%d, r = %d, mu = %f, lambda = %f, lambda1 = %f, times = %f s, Nmi = %f,acc = %f\n', step, r, mu, lambda, lambda1, t2, tempNMI, tempACC);
    
        end
    end
end
tmp = Paras{tag}{1};
fprintf('opt_Step (anchor_WASC): id = %d,r = %d, mu = %f, lambda = %f, lambda1 = %f,times = %f s,Nmi = %f, acc = %f \n \n',tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),tmp(7),tmp(8));

