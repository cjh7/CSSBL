clear;
rng(100);
% cd('C:\Users\jihoon7\Desktop\GSSBL - Copy - Copy_new\sim\propose\Result_Kaveh_case_modify_80_3');
for total=1:5
   
    G=2;
    M=12;
    N=33;
    leng=1;
    L1=50;
    L2=50;
    L=L1+L2;
    if total==1
        brho1= 0.1;
        brho2= 0.1;
    elseif total==2
        brho1= 0.3;
        brho2= 0.3;
    elseif total==3
        brho1= 0.6;
        brho2= 0.6;
    elseif total==4
        brho1= 0.9;
        brho2= 0.9;
    elseif total==5
        brho1= 0.95;
        brho2= 0.95;
     elseif total==6
        brho1= 0.4;
        brho2= 0.4;       
     elseif total==7
        brho1= 0.5;
        brho2= 0.5; 
     elseif total==8
        brho1= 0.7;
        brho2= 0.7;   
     elseif total==9
        brho1= 0.8;
        brho2= 0.8;
     elseif total==10
        brho1= 0.2;
        brho2= 0.2;           
        
  end 
  Btrue1 = [1 brho1 brho1 brho1 brho1 brho1; brho1 1 brho1 brho1 brho1 brho1; brho1 brho1 1 brho1 brho1 brho1; brho1 brho1 brho1 1 brho1 brho1; brho1 brho1 brho1 brho1 1 brho1; brho1 brho1 brho1 brho1 brho1 1];
  R1 = chol(Btrue1);
  Btrue2 = [1 brho2 brho2 ; brho2 1 brho2 ; brho2 brho2 1 ];
  R2 = chol(Btrue2);
    
    
for iteration=1:20
rng(iteration*total)

% data generation + sparse+ spatial correlated+ group (different sparsity)

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\Code')
%reorder = randi([1 33],1,33);
phi=csvread('phi_assembly_Kaveh_case.csv');
phi = phi./(ones(M,1)*sqrt(sum(phi.^2)));
%phi=phi(:,reorder);


var=1;
var_small=0.01;

%%%%%%%%independent fault
total_index=[1:N];
fault_ind=randsample([10:N],2);
fault_ind=sort(fault_ind);
%%%%%%%%%%1st group

Wgen = zeros(N,L1);
indice1=[1,2,3,4,5,6];
for i=1:L1
    r = normrnd(0,var,[1,6]);
    nonzeroW=r*R1*(sqrt(var));
    Wgen(indice1,i)=nonzeroW;
end

indice2=[7,8,9];
for i=1:L2
    r = normrnd(0,var,[1,3]);
    nonzeroW=r*R2*(sqrt(var_small));
    Wgen(indice2,i)=nonzeroW;
end    

% Wgen(fault_ind(1),:)= normrnd(0,var,[1,L1]);

for i=[setdiff(total_index,[indice1,indice2])]
    Wgen(i,:)=normrnd(0,var_small,[1,L1]);
end      
x1=Wgen;

%%%%%%%%%%%2nd group

Wgen = zeros(N,L2);
indice1=[1,2,3,4,5,6];
for i=1:L1
    r = normrnd(0,var,[1,6]);
    nonzeroW=r*R1*(sqrt(var));
    Wgen(indice1,i)=nonzeroW;
end

indice2=[7,8,9];
for i=1:L2
    r = normrnd(0,var,[1,3]);
    nonzeroW=r*R2*(sqrt(var));
    Wgen(indice2,i)=nonzeroW;
end

Wgen(fault_ind(1),:)= normrnd(0,var,[1,L1]);


for i=[setdiff(total_index,[indice1,indice2,fault_ind(1)])]
    Wgen(i,:)=normrnd(0,var_small,[1,L2]);
end           

x2=Wgen;

%%%%%%%%%%%%%combine both groups

x=[x1,x2];
y=phi*x + normrnd(0,0.001,[M,L1+L2]);

%%%%%%%%%% hyperparameters
a=1e-4;
b=1e-4;

% c &d 도 이용해야 함

c=1e-4;
d=1e-4;
%%%%%%%%%%% initializiation random variables

Z=rand(L,G);
Z= diag(  1./sum(Z,2) ) * Z;

lambda=1;
gamma=ones((N),G);
alpha_all=ones((N),L);
mu=zeros(N,1,L);
B1 = eye(1);
B2 = eye(1);

%%%%%%%% Phi matrix & Y (measurement; data)%%%%%%%%%%%%%%%

for l=1:L
    Phi(:,:,l) = phi;
end

for l=1:L  
    Y_all(:,:,l) = y(:,l);
end
%%%%%%%%%%% initialziation covariance structure
 for l=1:L
     for i=1:1
         cov_spatial(:,:,i) =gamma(i,:) * ( Z(l,:).') *B1;
     end
     for i=2:2
         cov_spatial(:,:,i) =gamma(i,:) * ( Z(l,:).') *B2;
     end
     for i=3:(N)
         cov_independent(i) =gamma(i,:) * ( Z(l,:).') *1;
     end

    for ii=1:(N)
        if ii==1
             covariance=cov_spatial(:,:,ii);
        elseif ii==2 
              covariance=blkdiag(covariance,cov_spatial(:,:,ii));
        else
            covariance=blkdiag(covariance,cov_independent(ii));
        end
    end
 covariance_all(:,:,l)=covariance;
 end

%%%%%%%%%%% initialziation mu 
for l=1:L
    Phi_delta = Phi(:,:,l) *  inv( covariance_all(:,:,l));
    V_temp= 1/lambda*eye(M) + Phi_delta * Phi(:,:,l)';
    Sigma(:,:,l) =  inv( covariance_all(:,:,l)) -Phi_delta' * (V_temp \Phi_delta);
    mu(:,:,l) = lambda * Sigma(:,:,l) * Phi(:,:,l)' * Y_all(:,:,l);
end


%%%%%%%%%%%%%%%%%%%%%%%% Algorithm stats

%%%%%% total number of iteration (maximum)
iter=5000;


for i1 = 1:iter   
  
%%%%%% To check the stopping criterion  
mu_previous=mu;   
if i1==1
    mu_previous=zeros(N,1,L);
else
    mu_previous=mu;
end 
  
%%%%%% Estimate Lambda  
temp0=zeros(M,1,L);
temp2=0;
for l=1:L
temp0(:,:,l)=Phi(:,:,l)*mu(:,:,l);
temp1=sum(diag(Phi(:,:,l)* Sigma(:,:,l)*Phi(:,:,l)'));
temp2=temp2+temp1;                  
end
resid=Y_all-temp0;
lambda=( L*M + a )/( b +  norm(resid(:), 'fro')^2  + (temp2) ); 


%%%%%% Estimate Covariance_all 
 for l=1:L
     for i=1:1
         cov_spatial(:,:,i) =gamma(i,:) * ( Z(l,:).') *B1;
     end
     for i=2:2
         cov_spatial(:,:,i) =gamma(i,:) * ( Z(l,:).') *B2;
     end
     for i=3:(N)
         cov_independent(i) =gamma(i,:) * ( Z(l,:).') *1;
     end

    for ii=1:(N)
        if ii==1
             covariance=cov_spatial(:,:,ii);
        elseif ii==2 
              covariance=blkdiag(covariance,cov_spatial(:,:,ii));
        else
            covariance=blkdiag(covariance,cov_independent(ii));
        end
    end
 covariance_all(:,:,l)=covariance;
 end

%%%%%% Estimate mu  
for l=1:L
    Phi_delta = Phi(:,:,l) *  inv( covariance_all(:,:,l));
    V_temp= 1/lambda*eye(M) + Phi_delta * Phi(:,:,l)';
    Sigma(:,:,l) =  inv( covariance_all(:,:,l)) -Phi_delta' * (V_temp \Phi_delta);
    mu(:,:,l) = lambda * Sigma(:,:,l) * Phi(:,:,l)' * Y_all(:,:,l);
end  

%%%%%% Estimate gamma  
temp0_gamma=zeros((N),G);
for i=1:2
  temp0_gamma(i,:) = 2*a-1 + sum(Z)*1*leng;
end
for i=3:(N)
  temp0_gamma(i,:) = 2*a-1 + sum(Z)*1*1;
end

temp1_gamma=zeros((N),G);
for g=1:G
    for i=1:1
        temp2_gamma=0;
    for l=1:L
       temp2_gamma=temp2_gamma+Z(l,g)*((mu((i-1)*leng+1:(i)*leng,:,l))'*B1*(mu((i-1)*leng+1:(i)*leng,:,l))+trace(B1*Sigma((i-1)*leng+1:(i)*leng,(i-1)*leng+1:(i)*leng,l))); 
    end
    temp1_gamma(i,g)=temp2_gamma;
    end
end

for g=1:G
    for i=2:2
        temp2_gamma=0;
    for l=1:L
   temp2_gamma=temp2_gamma+Z(l,g)*((mu((i-1)*leng+1:(i)*leng,:,l))'*B2*(mu((i-1)*leng+1:(i)*leng,:,l))+trace(B2*Sigma((i-1)*leng+1:(i)*leng,(i-1)*leng+1:(i)*leng,l))); 
    end
     temp1_gamma(i,g)=temp2_gamma;
    end
end

for g=1:G
    for i=3:(N)
        temp2_gamma=0;
    for l=1:L
       temp2_gamma=temp2_gamma+Z(l,g)*((mu((2*leng)+(i-2),l))'*(mu((2*leng)+(i-2),l))+trace(1*Sigma((2*leng)+(i-2),(2*leng)+(i-2),l)));  
    end
     temp1_gamma(i,g)=temp2_gamma;
    end
end

temp3_gamma=(temp1_gamma)+2*b;

for i=1:(N) 
for g=1:G
    gamma(i,g)=  temp0_gamma(i,g)/temp3_gamma(i,g);
    ln_gamma(i,g)= psi( temp0_gamma(i,g) ) -log(temp3_gamma(i,g));
end
end

%%%%%%%%%%%%%%% Estimate Z
for l=1:L
 for g=1:G   
    temp1_z=0;
    for i=1:1
    temp1_z=temp1_z+ ln_gamma(i,g)*leng+log(det(B1));
    end
    for i=2:2
    temp1_z=temp1_z+ ln_gamma(i,g)*leng+log(det(B2));
    end
    for i=3:(N)
    temp1_z=temp1_z+ ln_gamma(i,g);
    end
     
    temp2_z=0;
    for i=1:(N)
    if i==1
    temp2_z=temp2_z+  gamma(i,g)*((mu((i-1)*leng+1:(i)*leng,:,l))'*B1*(mu((i-1)*leng+1:(i)*leng,:,l))+trace(B1*Sigma((i-1)*leng+1:(i)*leng,(i-1)*leng+1:(i)*leng,l)));
    elseif i==2 
    temp2_z=temp2_z+  gamma(i,g)*((mu((i-1)*leng+1:(i)*leng,:,l))'*B2*(mu((i-1)*leng+1:(i)*leng,:,l))+trace(B2*Sigma((i-1)*leng+1:(i)*leng,(i-1)*leng+1:(i)*leng,l)));
    else
    temp2_z=temp2_z+  gamma(i,g)*((mu((2*leng)+(i-2),l))'*(mu((2*leng)+(i-2),l))+trace(1*Sigma((2*leng)+(i-2),(2*leng)+(i-2),l))); 
    end
    end
    
    xi(l,g)= temp1_z-temp2_z;
  end
end

for l=1:L
 temp(l,:)= exp(xi(l,:)-max(max(xi(l,:)))); 
end
Z= diag(  1./sum(temp,2) ) *   temp;



if norm(squeeze(mu_previous) - squeeze(mu) ,'fro') < 1e-6
    break;
end

end


for l=1:L
   estimate(:,l)=mu(1:N,l); 
end


resultdai(total,iteration)=(norm(x - estimate ,'fro')/norm(x,'fro'))^2
resultdai(total,iteration)
% Z



cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\Dai\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\Dai\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,estimate);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\Dai\Z');
filename = sprintf('Z_%d_%d.csv', total,iteration );
csvwrite(filename,Z);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\Dai\ITERATION');
filename = sprintf('iteration_%d_%d.csv', total,iteration );
csvwrite(filename,i1);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\Dai\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=1./gamma;
csvwrite(filename,variance_interest);


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\Dai\FAULT_INDEX');
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
csvwrite(filename,fault_ind);


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\Code')
end




clear nonzeroW;
clear Phi;
clear Phi_delta;
clear Sigma;
clear mu_previous;
clear ln_gamma;
clear estimate;
clear temp1_gamma
clear temp3_gamma;
clear temp0_gamma;
clear diag_inv;
clear lambda;
clear temp0;
clear temp1;
clear temp2;
clear Z;
mean(resultdai);
end


