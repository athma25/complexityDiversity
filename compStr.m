%Probability of diversification (from Doebeli and Ispolatov, 2010)

%PARAMETERS
%m: dimension of phenotype space
%[-d,d]: range of off diagonal elements
%a: stength of frquency dependence
%N: number of simulations
%A: competition matrix
%K: selection matrix
%p: probability vector

m=5;
d=.2;
N=10000;
a=.2:.01:.8;

A=zeros(m,m);
K=zeros(m,m);
p=zeros(length(a));

for k=1:length(a)
    %simulation for a(k)
    for l=1:N  
        %initialising A and K
        check=0;
        while(check<10)
            for i=1:m
                for j=1:m
                    if(i==j)
                        A(i,j)=a(k);
                        K(i,j)=1;
                    end
                    if(i<j)
                        A(i,j)=2*d*rand(1)-d;
                        K(i,j)=2*d*rand(1)-d;
                    end
                    if(j<i)
                        A(i,j)=A(j,i);
                        K(i,j)=K(j,i);
                    end
                end
            end
            Ea=eig(A);
            Ek=eig(K);
            check=0;
            for n=1:m
                if(Ea(n)>0)
                    check=check+1;
                end
                if(Ek(n)>0)
                    check=check+1;
                end
            end
        end 
        %we have s.p.d. A and K
        
        %diagonalising A 
        [U,D]=eig(A);
        A=transpose(U)*A*U;
        K=transpose(U)*K*U;
        %making A identity
        for i=1:m
            D(i,i)=1/sqrt(D(i,i));
        end
%        A=transpose(D)*A*D;
        K=transpose(D)*K*D;
        
        %diagonalising K
        [U,D]=eig(K);
        K=transpose(U)*K*U;
        i=1;
        check=0;
        while(i<=m && check==0)
            if(K(i,i)<1)
                p(k)=p(k)+1;
                check=1;
            end
            i=i+1;
        end
    end
    p(k)=p(k)/N;
end
plot(a,p)
        
    
