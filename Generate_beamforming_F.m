function [F_opt,power_opt,s_opt,flag] = Generate_beamforming_F(N,L, H, H_error,...
             prob, noise_maxpower, B, rate_min)

 
h=sort(abs(H)); 
cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable Gamma(N,N,L) hermitian
    variable s(L) 
    variable t(L)
    expressions    PHI(N,N,L)  constant(L) constraint_1(L)...
                   constraint_2(L)  constraint_3(N,N,L);
    
    for l=1:L
        f_l=[Gamma(:,:,l:L)];    %NOMA
        PHI(:,:,l)=1/(2^(rate_min/B)-1)*Gamma(:,:,l)-sum(f_l,3);
        constant(l)=(h(:,l))' * PHI(:,:,l) * h(:,l)...
                    - noise_maxpower;
        
        constraint_1(l)=H_error(l)^2*trace(PHI(:,:,l))-sqrt(2*log(1/prob))*s(l)...
                        +log(prob)*t(l)+constant(l);
        element_1=H_error(l)^2*vec(PHI(:,:,l));
        element_2=sqrt(2)*H_error(l)*PHI(:,:,l)*h(:,l);
        constraint_2(l)=norm([element_1;element_2],2);
        constraint_3(:,:,l)=t(l)*eye(N,N)+H_error(l)^2*PHI(:,:,l);
    end
    
    minimize trace(sum(Gamma,3))
  
    subject to
    
         real(constraint_1)>=0;
         real(constraint_2-s)<=0;
         t>=0;
    
         for l=1:L
             constraint_3(:,:,l)  == hermitian_semidefinite(N);
             Gamma(:,:,l) == hermitian_semidefinite(N);
              for i=l+1:L
                  Gamma(:,:,i) <=   Gamma(:,:,l);
              end
         end
     
cvx_end

if cvx_status(1)=='S' || cvx_status(3)=='a' 
    flag=1;

   
   %  test rate constraints  %
for l=1:L
    Matrix(:,:,l)=h(:,l)*(h(:,l))';
    PHI(:,:,l)=1/(2^(rate_min/B)-1)*Gamma(:,:,l)-sum(Gamma,3);
    Obj_original(l)=trace(PHI(:,:,l)*Matrix(:,:,l))-noise_maxpower;
end


for l=1:L
    [t1,t2]=eig(Gamma(:,:,l));
    [X,location]=max( abs(diag(t2)) );
    if length(location)==1
        F_opt(:,l)=t1(:,location)*t2(location,location)^(1/2);
    else
        for i=1:500
            b2(:,i)=t1*t2^(1/2)*sqrt(1/2)*(randn(N,1) + sqrt(-1)*  randn(N,1));
            Obj(i)=norm(b2(:,i),2);
        end
        [value locat]=min(Obj);
        F_opt(:,l)=b2(:,locat);
    end
end
power_opt=trace(F_opt*F_opt');


%%%%%  test rate constraints  %%%%%
for l=1:L
    F_opt_noma=F_opt(:,l:L);
    PHI(:,:,l)=1/(2^(rate_min/B)-1)*F_opt(:,l)*F_opt(:,l)'-F_opt_noma*F_opt_noma';
    Obj_new(l)=trace(PHI(:,:,l)*Matrix(:,:,l))-noise_maxpower;

    element_1=H_error(l)^2*vec(PHI(:,:,l));
    element_2=sqrt(2)*H_error(l)*PHI(:,:,l)*(h(:,l));
    s_opt(l)=norm([element_1;element_2],2);
end

else
    flag=0;
    F_opt=ones(N,L);
    power_opt=0;
    s_opt=0;
end

end

