function F=USVT(A)
n=size(A,1);
[V,E]=eigs(A,floor(sqrt(n)));
e=diag(E);
S=find(abs(e)>sqrt(n*mean(A(:)))*4);
length(S);
F=V(:,S)*E(S,S)*V(:,S)';
end