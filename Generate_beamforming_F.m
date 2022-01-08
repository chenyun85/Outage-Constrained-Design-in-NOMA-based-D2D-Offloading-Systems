function [F_opt,power_opt,x_opt,flag] = Generate_beamforming_F(N,K, H, H_error,...
             prob, noise_maxpower, B, rate_min)

 
h=sort(abs(H)); 
cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable Gamma(N,N,K) hermitian
    variable x(K) 
    variable y(K)
    expressions    PHI(N,N,K)  constant(K) constraint_1(K)...
                   constraint_2(K)  constraint_3(N,N,K);
    
    for k=1:K
        f_k=[Gamma(:,:,k:K)];    %NOMA
        PHI(:,:,k)=1/(2^(rate_min/B)-1)*Gamma(:,:,k)-sum(f_k,3);
        constant(k)=(h(:,k))' * PHI(:,:,k) * h(:,k)...
                    - noise_maxpower;
        
        constraint_1(k)=H_error(k)^2*trace(PHI(:,:,k))-sqrt(2*log(1/prob))*x(k)...
                        +log(prob)*y(k)+constant(k);
        element_1=H_error(k)^2*vec(PHI(:,:,k));
        element_2=sqrt(2)*H_error(k)*PHI(:,:,k)*h(:,k);
        constraint_2(k)=norm([element_1;element_2],2);
        constraint_3(:,:,k)=y(k)*eye(N,N)+H_error(k)^2*PHI(:,:,k);
    end
    
    minimize trace(sum(Gamma,3))
  
    subject to
    
         real(constraint_1)>=0;
         real(constraint_2-x)<=0;
         y>=0;
    
         for k=1:K
             constraint_3(:,:,k)  == hermitian_semidefinite(N);
             Gamma(:,:,k) == hermitian_semidefinite(N);
              for i=k+1:K
                  Gamma(:,:,i) <=   Gamma(:,:,k);
              end
         end
     
cvx_end

if cvx_status(1)=='S' || cvx_status(3)=='a' 
    flag=1;

   
   %  test rate constraints  %
for k=1:K
    Matrix(:,:,k)=h(:,k)*(h(:,k))';
    PHI(:,:,k)=1/(2^(rate_min/B)-1)*Gamma(:,:,k)-sum(Gamma,3);
    Obj_original(k)=trace(PHI(:,:,k)*Matrix(:,:,k))-noise_maxpower;
end


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
for k=1:K
    F_opt_noma=F_opt(:,k:K);
    PHI(:,:,k)=1/(2^(rate_min/B)-1)*F_opt(:,k)*F_opt(:,k)'-F_opt_noma*F_opt_noma';
    Obj_new(k)=trace(PHI(:,:,k)*Matrix(:,:,k))-noise_maxpower;

    element_1=H_error(k)^2*vec(PHI(:,:,k));
    element_2=sqrt(2)*H_error(k)*PHI(:,:,k)*(h(:,k));
    x_opt(k)=norm([element_1;element_2],2);
end

else
    flag=0;
    F_opt=ones(N,K);
    power_opt=0;
    x_opt=0;
end

end

