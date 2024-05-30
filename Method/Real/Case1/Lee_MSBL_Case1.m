clear all;  


rng(100);
for total=1:5
    run_MSBL=1; early_stop_MSBL=0;
    run_LEE=1; early_stop_LEE=0;
    G=1;
    M=12;
    N=33;
    leng=3;
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


%======================= MSBL =======================   
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\Code')
rng(iteration);
    if run_MSBL == 1,
        tic;
        [X, stop,gamma,count]= MSBL(phi,y);
        Weight_MSBL =X;
        time_MSBL = toc;
        TIME_MSBL(iteration) = time_MSBL;
    
    if (stop==1) early_stop_MSBL = early_stop_MSBL +1; end;
      
        mse_MSBL(iteration) = (norm(x - Weight_MSBL,'fro')/norm(x,'fro'))^2;  
         
     fprintf(' MSBL(learn lambda): time = %5.2f;  Ave-MSE = %3.2f%%; Ave-Time = %4.3f\n',...
         time_MSBL, mean(mse_MSBL)*100,mean(TIME_MSBL));
    disp(early_stop_MSBL)
    
    
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\MSBL\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\MSBL\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,Weight_MSBL);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\MSBL\ITERATION');
filename = sprintf('iteration_%d_%d.csv', total,iteration );
csvwrite(filename,count);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\MSBL\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=gamma;
csvwrite(filename,variance_interest);


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\MSBL\FAULT_INDEX');
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
csvwrite(filename,fault_ind);
    
    
    
    
    
    end

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\Code')
%======================= Lee =======================       
rng(iteration);
    if run_LEE == 1,
        tic;
        [X, stop,est_sig2a ]= LEE(phi,y);
        Weight_LEE =X;
        time_LEE = toc;
        TIME_LEE(iteration) = time_LEE;
    
    if (stop==1) early_stop_LEE = early_stop_LEE +1; end;
      
        mse_LEE(iteration) = (norm(x - Weight_LEE,'fro')/norm(x,'fro'))^2;  
         
     fprintf(' LEE(learn lambda): time = %5.2f;  Ave-LEE = %3.2f%%; Ave-Time = %4.3f\n',...
         time_LEE,mean(mse_LEE)*100,mean(TIME_LEE));
    disp(early_stop_LEE)
    
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\LEE\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\LEE\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,Weight_LEE);


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\LEE\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=est_sig2a;
csvwrite(filename,variance_interest);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\REAL\before_submit_kaveh_final\LEE\FAULT_INDEX');
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
csvwrite(filename,fault_ind);    
    
    
    
end    
    
end


%     result(3,9,total) =mean(mse_MSBL)*100;
%     result(3,10,total) =mean(TIME_MSBL); 
%     result(3,11,total) =mean(early_stop_MSBL); 
%     result(3,13,total) =std(mse_MSBL)*100;
% 
%     result(6,9,total) =mean(mse_LEE)*100;
%     result(6,10,total) =mean(TIME_LEE); 
%     result(6,11,total) =mean(early_stop_LEE);
%     result(6,13,total) =std(mse_LEE)*100;  

clear Wgen;

end