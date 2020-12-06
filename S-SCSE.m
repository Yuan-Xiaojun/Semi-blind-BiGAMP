function Pa_P = S-SCSE(Pa, P_num)   %
%% A is row,   B is column,  I is permutation case,   R is phase case
Pa_P = Pa;
[A,B,R] = size(P_num);
V = sort(randperm(B));      %  produce  a   array   with  1  to B
F = perms(V);         %all  permutation  case   in V
[I,~] = size(F);
Permu_Pro = ones(I,1);
Permu_Pro2 = ones(I,B);
Cp= ones(R,B,B);    % save the probablity of every column convert to all other column 
Cp2= zeros(R,B,B);
Pa_P(1:A,:,:) = 0;

%% the probability that b column change to p column
for b = 1: B     % b column
    for p = 1:B     % b column change to p column
        for r = 1:R     % all phrase
            for a = 1:A
                Cp(r,p,b) = Cp(r,p,b)*Pa(a,b,P_num(a,p,r));
            end
        end
    end
end
Cp = max(eps,Cp);
Cp1 = sum(Cp);

%%   normalized
for b = 1:B        
    for p = 1:B
        for r = 1:R
            Cp2(r,p,b) = Cp(r,p,b)/Cp1(1,p,b);
        end
    end
end

%%   the probability of all permutation
for i = 1:I
    for b = 1:B
        Permu_Pro(i,1) = Permu_Pro(i,1)*Cp1(1,F(i,b),b);
    end
end
%%
%%  normalized
Q = sum(Permu_Pro);
Permu_Pro(:,1) = Permu_Pro(:,1)/Q;
%%  calculate the probability of each note
for b = 1:B
    for a = 1:A
        for i = 1:I
            for r = 1:R
                Pa_P(a,b,P_num(a,F(i,b),r)) = Pa_P(a,b,P_num(a,F(i,b),r)) + Permu_Pro(i,1)*Cp2(r,F(i,b),b);
            end
        end
    end
end

%%   eliminated the exterior massage
% for b = 1:B
%     for i = 1 : I
%         Permu_Pro2(i,b) = Permu_Pro(i,1)/Cp1(1,F(i,b),b);
%     end
% end

%%  normalized
% Q = sum(Permu_Pro2);
% for b = 1:B
%     Permu_Pro2(:,b) = Permu_Pro2(:,b)/Q(1,b);
% end
%%  calculate the probability of each note
% for b = 1:B
%     for a = 1:A
%         for i = 1:I
%             for r = 1:R
%                 Pa_P(a,b,P_num(a,F(i,b),r)) = Pa_P(a,b,P_num(a,F(i,b),r)) + Permu_Pro2(i,b)*Cp2(r,F(i,b),b);
%             end
%         end
%     end
% end




