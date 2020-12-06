function P_num = Pilot_Num(D, Sam) 
[A,B]=size(D);
P_num = zeros(A,B,4);
S = [1,1i,-1,-1i];

%%
for r = 1:4
    D1 = D*S(1,r);
    for b = 1:B
        for a = 1:A
            [~,n] = min(abs(Sam-D1(a,b)));
            P_num(a,b,r) = n;
        end            
    end
end
