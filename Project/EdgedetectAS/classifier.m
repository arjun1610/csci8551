function value = classifier(P)
%% calculates the Threshold value for the Input P
global epsilon
epsilon=0.1;
P = P(:);
i=1;
Threshold(i)=mean(P);
Mean1= mean(P(P<Threshold(i)));
Mean2= mean(P(P>=Threshold(i)));
i=i+1;
Threshold(i)= (Mean1+Mean2)/2;
value=Threshold(i);
%change the epsilon value to get more/less pixels (noise issue)
while abs(Threshold(i)-Threshold(i-1))>=epsilon
    Mean1= mean(P(P<Threshold(i)));
    Mean2= mean(P(P>=Threshold(i)));
    i=i+1;
    Threshold(i)= (Mean1+Mean2)/2;
    value=Threshold(i);
end
