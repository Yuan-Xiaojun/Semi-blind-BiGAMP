function  [Aph, pMatrix2] = find_permu(A,P, Sam)
    [M,N]=size(A);
    Q = sqrt(length(Sam));
    tempBhat = zeros(Q*P,N);
    
    permu = zeros(N,N);
    A(A == 0) = Sam(1);
    if Q == 2
       n2=Sam(1)*ones(1,N)./A(1,:);
    else
       n2=Sam(9)*ones(1,N)./A(1,:);
    end
    
    Phase=diag(n2);
    Aph = A*Phase;
    
     for m=1:M
         for n=1:N
             [~,I]=min(abs(Sam-Aph(m,n)));
             Aph(m,n)=Sam(I);
          end
     end
     
       if length(Sam) == 4         
           for  n1=1:N
                for  m1=1:P               
                    tempBhat((2*m1-1),n1)=real(Aph(m1,n1));
                    tempBhat((2*m1),n1)=imag(Aph(m1,n1));
                end
            end
            tempBhat(tempBhat>0) = 1;
            tempBhat(tempBhat<0) = 0;
        elseif length(Sam) == 16
            for  n=1:N
                for  m=1:P
                    if  Aph(m,n) == Sam(1)
                        tempBhat((4*m-3),n) = 1;
                        tempBhat((4*m-2),n) = 0;
                        tempBhat((4*m-1),n) = 1;
                        tempBhat((4*m),n) = 0;
                   elseif  Aph(m,n) == Sam(2) 
                        tempBhat((4*m-3),n) = 0;
                        tempBhat((4*m-2),n) = 1;
                        tempBhat((4*m-1),n) = 1;
                        tempBhat((4*m),n) = 0;
                    elseif  Aph(m,n) == Sam(3)
                        tempBhat((4*m-3),n) = 0;
                        tempBhat((4*m-2),n) = 1;
                        tempBhat((4*m-1),n) = 0;
                        tempBhat((4*m),n) = 1;
                    elseif  Aph(m,n) == Sam(4)
                        tempBhat((4*m-3),n) = 1;
                        tempBhat((4*m-2),n) = 0;
                        tempBhat((4*m-1),n) = 0;
                        tempBhat((4*m),n) = 1;
                    elseif  Aph(m,n) == Sam(5)
                        tempBhat((4*m-3),n) = 1;
                        tempBhat((4*m-2),n) = 1;
                        tempBhat((4*m-1),n) = 1;
                        tempBhat((4*m),n) = 0;
                   elseif  Aph(m,n) == Sam(6) 
                        tempBhat((4*m-3),n) = 0;
                        tempBhat((4*m-2),n) =1;
                        tempBhat((4*m-1),n) = 1;
                        tempBhat((4*m),n) = 1;
                   elseif  Aph(m,n) == Sam(7) 
                        tempBhat((4*m-3),n) = 0;
                        tempBhat((4*m-2),n) = 0;
                        tempBhat((4*m-1),n) = 0;
                        tempBhat((4*m),n) = 1;
                   elseif  Aph(m,n) == Sam(8)                         
                        tempBhat((4*m-3),n) = 1;
                        tempBhat((4*m-2),n) = 0;
                        tempBhat((4*m-1),n) = 0;
                        tempBhat((4*m),n) = 0;
                   elseif  Aph(m,n) == Sam(9)                        
                        tempBhat((4*m-3),n) = 1;
                        tempBhat((4*m-2),n) = 1;
                        tempBhat((4*m-1),n) = 1;
                        tempBhat((4*m),n) = 1;
                   elseif  Aph(m,n)== Sam(10)                      
                        tempBhat((4*m-3),n) = 0;
                        tempBhat((4*m-2),n) = 0;
                        tempBhat((4*m-1),n) = 1;
                        tempBhat((4*m),n) = 1;
                   elseif  Aph(m,n) == Sam(11)   
                        tempBhat((4*m-3),n) = 0;
                        tempBhat((4*m-2),n) = 0;
                        tempBhat((4*m-1),n) = 0;
                        tempBhat((4*m),n) = 0;
                   elseif  Aph(m,n) == Sam(12) 
                        tempBhat((4*m-3),n) = 1;
                        tempBhat((4*m-2),n) = 1;
                        tempBhat((4*m-1),n) = 0;
                        tempBhat((4*m),n) = 0;
                   elseif  Aph(m,n) == Sam(13) 
                        tempBhat((4*m-3),n) = 1;
                        tempBhat((4*m-2),n) = 0;
                        tempBhat((4*m-1),n) = 1;
                        tempBhat((4*m),n) = 1;
                   elseif  Aph(m,n) == Sam(14) 
                        tempBhat((4*m-3),n) = 0;
                        tempBhat((4*m-2),n) = 0;
                        tempBhat((4*m-1),n) = 1;
                        tempBhat((4*m),n) = 0;
                   elseif  Aph(m,n) == Sam(15) 
                        tempBhat((4*m-3),n) = 0;
                        tempBhat((4*m-2),n) = 1;
                        tempBhat((4*m-1),n) = 0;
                        tempBhat((4*m),n) = 0;
                     else   
                        tempBhat((4*m-3),n) = 1;
                        tempBhat((4*m-2),n) = 1;
                        tempBhat((4*m-1),n) = 0;
                        tempBhat((4*m),n) = 1;
                    end
                end  
            end    
        end

for n=1:N
    n3 = 0;
  for m=(1+Q) : (Q+ceil(log2(N)))       
       n3 = n3 + 2^(m-1-Q)*tempBhat(m,n); 
  end 
  if (0<n3)&& (n3<=N)
     permu(n,n3)=1;
  end
  if (n3==0)
     permu(n,N)=1;
  end
end
Aph = Aph * permu ;
pMatrix2 = Phase*permu;

