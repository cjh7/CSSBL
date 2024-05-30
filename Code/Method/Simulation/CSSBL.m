clear;
rng(100);
% cd('C:\Users\jihoon7\Desktop\gssbl_new\sim\propose\Result_Z2_60_60_modify_80_3');
for total=1:10
   
    G=2;
    M=8;
    N=40;
    leng=3;
    L1=60;
    L2=60;
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
        brho1= 0.2;
        brho2= 0.2;
    elseif total==7
        brho1= 0.4;
        brho2= 0.4;
    elseif total==8
        brho1= 0.5;
        brho2= 0.5;
    elseif total==9
        brho1= 0.7;
        brho2= 0.7;
    elseif total==10
        brho1= 0.8;
        brho2= 0.8;        
  end 
  Btrue1 = [1 brho1 brho1 ; brho1 1 brho1 ; brho1 brho1 1 ];
  R1 = chol(Btrue1);
  Btrue2 = [1 brho2 brho2 ; brho2 1 brho2 ; brho2 brho2 1 ];
  R2 = chol(Btrue2);
    
for iteration=1:20
rng(iteration*total)

% data generation + sparse+ spatial correlated+ group (different sparsity)

phi = randn(M,N);
phi = phi./(ones(M,1)*sqrt(sum(phi.^2)));

var=1;
var_small=0.01;

%%%%%%%%independent fault
total_index=[1:N];
fault_ind=randsample([7:N],6);
fault_ind=sort(fault_ind);
%%%%%%%%%%1st group

Wgen = zeros(N,L1);
indice1=[1,2,3];
for i=1:L1
    r = normrnd(0,var,[1,3]);
    nonzeroW=r*R1*(sqrt(var));
    Wgen(indice1,i)=nonzeroW;
end

indice2=[4,5,6];
for i=1:L1
    r = normrnd(0,var,[1,3]);
    nonzeroW=r*R2*(sqrt(var_small));
    Wgen(indice2,i)=nonzeroW;
end    

Wgen(fault_ind(1),:)= normrnd(0,var,[1,L1]);
Wgen(fault_ind(3),:)= normrnd(0,var,[1,L1]);
Wgen(fault_ind(5),:)= normrnd(0,var,[1,L1]);

for i=[setdiff(total_index,[indice1,indice2,fault_ind(1),fault_ind(3),fault_ind(5)])]
    Wgen(i,:)=normrnd(0,var_small,[1,L1]);
end      
x1=Wgen;

%%%%%%%%%%%2nd group

Wgen = zeros(N,L2);
indice1=[1,2,3];
for i=1:L2
    r = normrnd(0,var,[1,3]);
    nonzeroW=r*R1*(sqrt(var_small));
    Wgen(indice1,i)=nonzeroW;
end

indice2=[4,5,6];
for i=1:L2
    r = normrnd(0,var,[1,3]);
    nonzeroW=r*R2*(sqrt(var));
    Wgen(indice2,i)=nonzeroW;
end

Wgen(fault_ind(2),:)= normrnd(0,var,[1,L2]); 
Wgen(fault_ind(4),:)= normrnd(0,var,[1,L2]);
Wgen(fault_ind(6),:)= normrnd(0,var,[1,L2]);

for i=[setdiff(total_index,[indice1,indice2,fault_ind(2),fault_ind(4),fault_ind(6)])]
    Wgen(i,:)=normrnd(0,var_small,[1,L2]);
end           

x2=Wgen;



%%%%%%%%%%%%%combine both groups

x=zeros(N,L);
r = randperm(L);
x(:,r(1:L1))=x1;
x(:,r(L1+1:L1+L2))=x2;
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
gamma=ones((N-4),G);
mu=zeros(N,1,L);
B1 = eye(3);
B2 = eye(3);

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
     for i=3:(N-4)
         cov_independent(i) =gamma(i,:) * ( Z(l,:).') *1;
     end

    for ii=1:(N-4)
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
     for i=3:(N-4)
         cov_independent(i) =gamma(i,:) * ( Z(l,:).') *1;
     end

    for ii=1:(N-4)
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
temp0_gamma=zeros((N-4),G);
for i=1:2
  temp0_gamma(i,:) = 2*a-1 + sum(Z)*1*leng;
end
for i=3:(N-4)
  temp0_gamma(i,:) = 2*a-1 + sum(Z)*1*1;
end

temp1_gamma=zeros((N-4),G);
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
    for i=3:(N-4)
        temp2_gamma=0;
    for l=1:L
       temp2_gamma=temp2_gamma+Z(l,g)*((mu((2*leng)+(i-2),l))'*(mu((2*leng)+(i-2),l))+trace(1*Sigma((2*leng)+(i-2),(2*leng)+(i-2),l)));  
    end
     temp1_gamma(i,g)=temp2_gamma;
    end
end

temp3_gamma=(temp1_gamma)+2*b;

for i=1:(N-4) 
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
    for i=3:(N-4)
    temp1_z=temp1_z+ ln_gamma(i,g);
    end
     
    temp2_z=0;
    for i=1:(N-4)
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

%%%%%%%%%%%%%%% Estimate B1
for l=1:L
for g=1:G
    temp0_B1=0;
    for i=1:1
     temp0_B1=temp0_B1+gamma(i,g)*( Sigma((i-1)*leng+1:(i)*leng,(i-1)*leng+1:(i)*leng,l) + mu((i-1)*leng+1:(i)*leng,:,l)*mu((i-1)*leng+1:(i)*leng,:,l)');
    end
    temp1_B1(:,:,l,g)= temp0_B1;
end
end

temp2_B1=zeros(3);
for l=1:L
    for g=1:G
      temp2_B1=temp2_B1+ Z(l,g)*temp1_B1(:,:,l,g);
    end
end
temp3_B1=(temp2_B1);
temp4_B1=1*sum(sum(Z));
temp5_B1=temp3_B1/temp4_B1;
temp6_B1=trace(temp5_B1)/3;
temp7_B1=(temp5_B1(2,1)+temp5_B1(3,2))/2;
temp8_B1=temp7_B1/temp6_B1;
for i=1:leng
    for j=1:leng
temp5_B1(i,j)=temp8_B1;
    end
end
temp5_B1(1,1)=1;
temp5_B1(2,2)=1;
temp5_B1(3,3)=1;
d = eig(temp5_B1);
if d(1)<0
    temp5_B1=temp5_B1+2*eye(3);
end
B1=inv(temp5_B1);


%%%%%%%%%%%%%%% Estimate B2
for l=1:L
for g=1:G
    temp0_B2=0;
    for i=2:2
     temp0_B2=temp0_B2+gamma(i,g)*( Sigma((i-1)*leng+1:(i)*leng,(i-1)*leng+1:(i)*leng,l) + mu((i-1)*leng+1:(i)*leng,:,l)*mu((i-1)*leng+1:(i)*leng,:,l)');
    end
    temp1_B2(:,:,l,g)= temp0_B2;
end
end

temp2_B2=zeros(3);
for l=1:L
    for g=1:G
      temp2_B2=temp2_B2+ Z(l,g)*temp1_B2(:,:,l,g);
    end
end
temp3_B2=(temp2_B2);
temp4_B2=1*sum(sum(Z));
temp5_B2=temp3_B2/temp4_B2;
temp6_B2=trace(temp5_B2)/3;
temp7_B2=(temp5_B2(2,1)+temp5_B2(3,2))/2;
temp8_B2=temp7_B2/temp6_B2;
for i=1:leng
    for j=1:leng
temp5_B2(i,j)=temp8_B2;
    end
end
temp5_B2(1,1)=1;
temp5_B2(2,2)=1;
temp5_B2(3,3)=1;
d = eig(temp5_B2);
if d(1)<0
    temp5_B2=temp5_B2+2*eye(3);
end
B2=inv(temp5_B2);


if norm(squeeze(mu_previous) - squeeze(mu) ,'fro') < 1e-6
    break;
end

end


for l=1:L
   estimate(:,l)=mu(1:N,l); 
end


resultpr(total,iteration)=(norm(x - estimate ,'fro')/norm(x,'fro'))^2
resultpr(total,iteration)
Z
inv(B1)
inv(B2)
1./gamma




cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,estimate);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\Z');
filename = sprintf('Z_%d_%d.csv', total,iteration );
csvwrite(filename,Z);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\ITERATION');
filename = sprintf('iteration_%d_%d.csv', total,iteration );
csvwrite(filename,i1);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=1./gamma;
csvwrite(filename,variance_interest);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\R1');
filename = sprintf('r1_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(1:L1)));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\R2');
filename = sprintf('r2_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(L1+1:L)));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\B1');
filename = sprintf('B1_%d_%d.csv', total,iteration );
csvwrite(filename, inv(B1));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\B2');
filename = sprintf('B2_%d_%d.csv', total,iteration );
csvwrite(filename, inv(B2));


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\Propose\FAULT_INDEX');
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
csvwrite(filename,fault_ind);




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
mean(resultpr);
end
% mean(resultpr,2)
% csvwrite('proposed_simulation_1_5_9_modify80_3.csv',resultpr);

