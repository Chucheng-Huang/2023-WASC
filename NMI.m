function MIhat = NMI( A, B ) 



if length( A ) ~= length( B)
    error('length( A ) must == length( B)');
end
total = length(A);
A_ids = unique(A);
B_ids = unique(B);


MI = 0;

for i=1:1:length(A_ids)
    for j=1:1:length(B_ids)
         idAOccur = find( A == A_ids(i) );
         idBOccur = find( B == B_ids(j) );
         idABOccur = intersect(idAOccur,idBOccur); 
         
         px = length(idAOccur)/total;
         py = length(idBOccur)/total;
         pxy = length(idABOccur)/total;
         
         MI = MI + pxy*log2(pxy/(px*py)+eps); % eps : the smallest positive number

    end
end

% Normalized Mutual information
Hx = 0; % Entropies
%for idA = A_ids
for i=1:1:length(A_ids)
    idAOccurCount = length( find( A == A_ids(i) ) );
    Hx = Hx - (idAOccurCount/total) * log2(idAOccurCount/total + eps);
    %Hx = Hx + (idAOccurCount/total) * log2(idAOccurCount/total + eps);
end
Hy = 0; % Entropies
%for idB = B_ids
for i=1:1:length(B_ids)
    idBOccurCount = length( find( B == B_ids(i) ) );
    Hy = Hy - (idBOccurCount/total) * log2(idBOccurCount/total + eps);
    %Hy = Hy + (idBOccurCount/total) * log2(idBOccurCount/total + eps);
end

MIhat = 2 * MI / (Hx+Hy);
%MIhat = MI / sqrt(Hx*Hy);
end
