function P = levina_smooth(A,q)
A2 = A^2;
A2 = A2 - diag(diag(A2));
n = size(A,1);
for i = 1:n
    for j = 1:n
        vi = A2(i,:);
        vj = A2(j,:);
        vi(j) = 0;
        vj(i) = 0;
        D(i,j) = max(abs(vi - vj))/n;
        D(j,i) = D(i,j);
    end  
end

for i = 1:n
    Ni = find(D(i,:) < quantile(D(i,:),q));
    Ni = setdiff(Ni,i);
    if(isempty(Ni))
    	for j = (i+1):n
            Nj = find(D(j,:) < quantile(D(j,:),q));
            Nj = setdiff(Nj,j);
            if(isempty(Nj))
        		P(i,j) = 0;
        		P(j,i) = P(i,j);
    	    else
	    	    P(i,j) = mean(A(i,Nj));
           	    P(j,i) = P(i,j);
            end
        end
    else
        for j = (i+1):n
            Nj = find(D(j,:) < quantile(D(j,:),q));
            Nj = setdiff(Nj,j);
            if(isempty(Nj))
        		P(i,j) = mean(A(Ni,j));
        		P(j,i) = P(i,j);
    	    else
	    	    P(i,j) = 0.5*(mean(A(Ni,j))+mean(A(i,Nj)));
           	    P(j,i) = P(i,j);
            end
        end
    end
end        
%keyboard
end

