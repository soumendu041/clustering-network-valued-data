function x = nmi(C,Chat)
N = length(C);
CF = confusionmat(C,Chat);
p1 = sum(CF,1)/N;
p2 = sum(CF,2)/N;
p12 = reshape(CF,[],1)/N;
%keyboard
H1 = entropy(p1);
H2 = entropy(p2);
H12 = entropy(p12);
x = 2*(1 - H12/(H1 + H2));
end