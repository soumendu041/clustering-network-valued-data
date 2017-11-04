% Parameters for the block model graphon (See Levina et al. (neighbourhood smoothing paper) for definition)
K = 2;
c = 0.3;
str=2; %Set to 'fro' for frobenius
% Latent variables
n = 100;
Xi = rand(n,1);

[U,V] = meshgrid(Xi,Xi);

% Edge probabilities
P1 = arrayfun(@graphon1,U,V,K*ones(size(U)),c*ones(size(U)));
P2 = arrayfun(@graphon2,U,V);
%P3 = arrayfun(@graphon3,U,V);
%P4 = arrayfun(@graphon2,U,V);

k=3;
rho=.3;
eta=.5;
p=1; q=p*eta;
[~,~,P3,~]=create_block_model(n,rho,(p-q)*eye(k)+q*ones(k,k),ones(1,k)/k);
%p=1.3;q=.13;
diff = .1;
p1=p*(1+diff);
q1=q*(1+diff);
[~,~,P4,~]=create_block_model(n,rho,(p1-q1)*eye(k)+q1*ones(k,k),ones(1,k)/k);

% T networks sampled from the above four graphons with some specified
% frequency, i.e. a mixture model

T = 50;
Networks = zeros(n,n,T);

freq = [0.0,0.0,0.5,0.5]; % mixing distribution

pd = makedist('Multinomial','Probabilities',freq);
class = zeros(T,1); % true class labels

for t = 1:T
    class(t) = random(pd);
    P = 'Pr';
    P = strrep(P,'r',num2str(class(t)));
    temp = arrayfun(@binornd,ones(n),eval(P));
    Networks(:,:,t) = triu(temp)+triu(temp,1)'; 
end

% Estimation of the matrix

Est = zeros(n,n,T);
nmotifs = 3; 
countstat = zeros(nmotifs,T);

for t = 1:T
    Est(:,:,t) = USVT(Networks(:,:,t)); % substitute method by intended function
    countstat(:, t) = moments(Networks(:,:,t),nmotifs,'approx');
end

% Eucledian distance matrix between estimated graphons (Frobenius distance) 

DM_usvt = zeros(t,t);
DM_counts = zeros(t,t);
for i = 1:(T-1)
    for j = (i+1):T
       
        DM_usvt(i,j) = norm(Est(:,:,i) - Est(:,:,j), str);
        DM_counts(i,j) = norm(countstat(:,i) - countstat(:,j), str);
        DM_usvt(j,i) = DM_usvt(i,j);
        DM_counts(j,i) = DM_counts(i,j);
    end
end

% Now we do spectral clustering on distance matrix

[V_counts, Lambda_counts] = eigs(DM_counts,K);
comm_counts = mykmeans1(V_counts,K);
err_counts = cluster_acc(comm_counts,class)

[V_usvt, Lambda_usvt] = eigs(DM_usvt,K);
comm_usvt = mykmeans1(V_usvt,K);
err_usvt = cluster_acc(comm_usvt,class)