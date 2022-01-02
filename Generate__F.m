function [F_opt,power_opt,flag] = Generate__F(N,K, H)

 
h=sort(abs(H)); 
cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable Gamma(N,N,K) hermitian
 

    
    minimize trace(sum(Gamma,3))
  
    subject to
    
    
         for k=1:K

              for i=k+1:K
                  Gamma(:,:,i) <=   Gamma(:,:,k);
              end
         end
     
cvx_end

if cvx_status(1)=='S' || cvx_status(3)=='a' 
    flag=1;

   
   

 for k=1:K
    [t1,t2]=eig(Gamma(:,:,k));
    [X,location]=max( abs(diag(t2)) );
    if length(location)==1
        F_opt(:,k)=t1(:,location)*t2(location,location)^(1/2);
    else
        for i=1:500
            b2(:,i)=t1*t2^(1/2)*sqrt(1/2)*(randn(N,1) + sqrt(-1)*  randn(N,1));
            Obj(i)=norm(b2(:,i),2);
        end
        [value locat]=min(Obj);
        F_opt(:,k)=b2(:,locat);
    end
end
power_opt=trace(F_opt*F_opt');


%%%%%  test rate constraints  %%%%%

%%%%%          end            %%%%%

else
    flag=0;
    F_opt=ones(N,K);
    power_opt=0;
    x_opt=0;
end

end

