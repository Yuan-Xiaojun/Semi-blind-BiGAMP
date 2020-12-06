function [Pa_P] = SCSE(Pa, P_num)   %
Pa_P = Pa;
[A,B,R]=size(P_num);
Pa_P(1:A,:,:) = 0;
Cp= ones(R,B,B);    % save the probablity of every column convert to all other column

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
%Cp = max(eps,Cp);
Cp1 = sum(sum(Cp));

%% normalized
for b = 1: B
    Cp(:,:,b) = Cp(:,:,b)/Cp1(1,1,b);
end

%%  calculate the probability of each note
for b = 1:B
    for a = 1:A
        for p = 1:B
            for r = 1:R
                Pa_P(a,b,P_num(a,p,r)) = Pa_P(a,b,P_num(a,p,r)) + Cp(r,p,b);
            end
        end
    end
end
%%
% Pa1 = sum(Pa_P(1:A, :, :), 3);
% for b = 1: B
%    for a = 1 : A
%        Pa_P(a, b, :) = Pa_P(a, b, :)/Pa1(a,b); 
%    end
% end







%%

% for l = 1:B
%     CP = ones(R,B);
%     for r=1:R
%         for b = 1:B
%             for a = 1:A  
%                 CP(r,b) = CP(r,b)*Pa(a,l, P_num(a,b,r));
%             end
%         end
%     end
%     
%     for a = 1:A
%         Pall = 0;
%         for r=1:R
%             for b = 1:B
%                 Pa_P(a, l, P_num(a,b,r)) = Pa_P(a ,l ,P_num(a,b,r)) + CP(r,b); 
%             end
%         end
%         for i = 1:R
%             Pall = Pall + Pa_P(a,l,i);
%         end
%         for i=1:R
%             Pa_P(a,l,i) = Pa_P(a,l,i)/Pall;
%         end
%         
%     end
% end

