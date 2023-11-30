function [B] = norm_AnchorGraphB(Z)
    d = sum(Z);
    d_norm = 1./(sqrt(d)+ eps);
    B = Z*diag(d_norm);
end

