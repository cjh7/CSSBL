clear all;  


rng(100);
for total=1:10
    run_MSBL=1; early_stop_MSBL=0;
    run_LEE=1; early_stop_LEE=0;
    run_SA_TSBL=0; early_stop_SA_TSBL=0;
    run_SA_MSBL=0; early_stop_SA_MSBL=0;
    run_TMSBL=0; early_stop_TMSBL=0;
    G=3;
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


%======================= MSBL =======================   
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\Code')
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
    
    
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\MSBL\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\MSBL\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,Weight_MSBL);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\MSBL\ITERATION');
filename = sprintf('iteration_%d_%d.csv', total,iteration );
csvwrite(filename,count);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\MSBL\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=gamma;
csvwrite(filename,variance_interest);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\MSBL\R1');
filename = sprintf('r1_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(1:L1)));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\MSBL\R2');
filename = sprintf('r2_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(L1+1:L1+L2)));


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\MSBL\FAULT_INDEX');
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
csvwrite(filename,fault_ind);
    
    
    
    
    
    end

    %======================= SA_MSBL =======================   
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\Code')
rng(iteration);
    if run_SA_MSBL == 1,
        tic;
        [X, stop,gamma,count]= SA_MSBL(phi,y);
        Weight_SA_MSBL =X;
        time_SA_MSBL = toc;
        TIME_SA_MSBL(iteration) = time_SA_MSBL;
    
    if (stop==1) early_stop_SA_MSBL = early_stop_SA_MSBL +1; end;
      
        mse_SA_MSBL(iteration) = (norm(x - Weight_SA_MSBL,'fro')/norm(x,'fro'))^2;  
         
     fprintf(' MSBL(learn lambda): time = %5.2f;  Ave-MSE = %3.2f%%; Ave-Time = %4.3f\n',...
         time_SA_MSBL, mean(mse_SA_MSBL)*100,mean(TIME_SA_MSBL));
    disp(early_stop_SA_MSBL)
    
    
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SA_MSBL\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SA_MSBL\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,Weight_MSBL);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SA_MSBL\ITERATION');
filename = sprintf('iteration_%d_%d.csv', total,iteration );
csvwrite(filename,count);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SA_MSBL\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=gamma;
csvwrite(filename,variance_interest);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SA_MSBL\R1');
filename = sprintf('r1_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(1:L1)));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SA_MSBL\R2');
filename = sprintf('r2_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(L1+1:L1+L2)));


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SA_MSBL\FAULT_INDEX');
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
csvwrite(filename,fault_ind);
    

    end

    %======================= TMSBL =======================   
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\Code')
rng(iteration);
    if run_TMSBL == 1,
        tic;
        [X, stop,gamma,count]= TMSBL(phi,y);
        Weight_TMSBL =X;
        time_TMSBL = toc;
        TIME_TMSBL(iteration) = time_TMSBL;
    
    if (stop==1) early_stop_TMSBL = early_stop_TMSBL +1; end;
      
        mse_TMSBL(iteration) = (norm(x - Weight_TMSBL,'fro')/norm(x,'fro'))^2;  
         
     fprintf(' TMSBL(learn lambda): time = %5.2f;  Ave-MSE = %3.2f%%; Ave-Time = %4.3f\n',...
         time_TMSBL, mean(mse_TMSBL)*100,mean(TIME_TMSBL));
    disp(early_stop_TMSBL)
    
    
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\TMSBL\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\TMSBL\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,Weight_TMSBL);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\TMSBL\ITERATION');
filename = sprintf('iteration_%d_%d.csv', total,iteration );
csvwrite(filename,count);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\TMSBL\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=gamma;
csvwrite(filename,variance_interest);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\TMSBL\R1');
filename = sprintf('r1_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(1:L1)));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\TMSBL\R2');
filename = sprintf('r2_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(L1+1:L1+L2)));


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\TMSBL\FAULT_INDEX');
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
csvwrite(filename,fault_ind);
    

    end
    
    

%======================= SA_TSBL =======================   
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\Code')
rng(iteration);
    if run_SA_TSBL == 1,
        tic;
        [X, stop,gamma1,count]= SA_TSBL(phi,y);
        Weight_SA_TSBL =X;
        time_SA_TSBL = toc;
        TIME_SA_TSBL(iteration) = time_SA_TSBL;
    
    if (stop==1) early_stop_SA_TSBL = early_stop_SA_TSBL +1; end;
      
        mse_SA_TSBL(iteration) = (norm(x - Weight_SA_TSBL,'fro')/norm(x,'fro'))^2;  
         
     fprintf(' MSBL(learn lambda): time = %5.2f;  Ave-MSE = %3.2f%%; Ave-Time = %4.3f\n',...
         time_SA_TSBL, mean(mse_SA_TSBL)*100,mean(TIME_SA_TSBL));
    disp(early_stop_SA_TSBL)
    
    
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SATSBL\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SATSBL\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,Weight_MSBL);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SATSBL\ITERATION');
filename = sprintf('iteration_%d_%d.csv', total,iteration );
csvwrite(filename,count);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SATSBL\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=gamma1;
csvwrite(filename,variance_interest);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SATSBL\R1');
filename = sprintf('r1_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(1:L1)));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SATSBL\R2');
filename = sprintf('r2_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(L1+1:L1+L2)));


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\SATSBL\FAULT_INDEX');
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
csvwrite(filename,fault_ind);
    
    
    
    
    
    end




cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\Code')
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
    
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\LEE\X');
filename = sprintf('x_%d_%d.csv', total,iteration );
csvwrite(filename,x);
cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\LEE\ESTIMATE');
filename = sprintf('estimate_%d_%d.csv', total,iteration );
csvwrite(filename,Weight_LEE);


cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\LEE\GAMMA');
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance_interest=est_sig2a;
csvwrite(filename,variance_interest);

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\LEE\R1');
filename = sprintf('r1_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(1:L1)));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\LEE\R2');
filename = sprintf('r2_%d_%d.csv', total,iteration );
csvwrite(filename,sort(r(L1+1:L1+L2)));

cd('C:\Users\User\Desktop\논문\TASE\TASE_CSSBL_FINAL\TASE_CSSBL_FINAL\GSSBL\GSSBL\Sim\result_before_submit_final\LEE\FAULT_INDEX');
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