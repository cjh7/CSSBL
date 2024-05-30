clear;
warning('off','all')
L1=60;
L2=60;
N=40;
AUC_RESULT_FINAL=zeros(5,20);
NMSE_RESULT_FINAL=zeros(5,20);
ITER_RESULT_FINAL=zeros(5,20);
ESTIMATE_RESULT=zeros(5,20);
for total=1:10
    disp(total)
   for iteration=1:20
disp(iteration)
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Simulation\CSSBL\ESTIMATE');       
filename = sprintf('estimate_%d_%d.csv', total,iteration );
estimate=csvread(filename);       
       
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Simulation\CSSBL\X');       
filename = sprintf('x_%d_%d.csv', total,iteration );
x=csvread(filename);   


cd('C:\Users\User\Desktop\새 폴더\Code\Result\Simulation\CSSBL\R1');       
filename = sprintf('r1_%d_%d.csv', total,iteration );
r1=csvread(filename);  

cd('C:\Users\User\Desktop\새 폴더\Code\Result\Simulation\CSSBL\R2');       
filename = sprintf('r2_%d_%d.csv', total,iteration );
r2=csvread(filename);  

cd('C:\Users\User\Desktop\새 폴더\Code\Result\Simulation\CSSBL\Z');       
filename = sprintf('z_%d_%d.csv', total,iteration );
z=csvread(filename);  
       
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Simulation\CSSBL\GAMMA');       
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance=csvread(filename);
       
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Simulation\CSSBL\FAULT_INDEX');       
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
fault_index=csvread(filename);


cd('C:\Users\User\Desktop\새 폴더\Code\Result\Simulation\CSSBL\ITERATION')
filename = sprintf('iteration_%d_%d.csv', total,iteration );
iter_result=csvread(filename);

iter=0;

% Label Information
label_1=zeros(N-4,1);
label_1(1)=1;
label_1(fault_index(1)-4)=1;
label_1(fault_index(3)-4)=1;
label_1(fault_index(5)-4)=1;

label_2=zeros(N-4,1);
label_2(2)=1;
label_2(fault_index(2)-4)=1;
label_2(fault_index(4)-4)=1;
label_2(fault_index(6)-4)=1;

% Variance Information
var_1=ones(N-4,1)*0.01;
var_1(1)=1;
var_1(fault_index(1)-4)=1;
var_1(fault_index(3)-4)=1;
var_1(fault_index(5)-4)=1;

var_2=ones(N-4,1)*0.01;
var_2(2)=1;
var_2(fault_index(2)-4)=1;
var_2(fault_index(4)-4)=1;
var_2(fault_index(6)-4)=1;

AUC_RESULT=0;
NMSE_RESULT=0;
for i=1:120

  estimation_variance=z(i,1)*variance(:,1)+z(i,2)*variance(:,2);
  
  if ismember(i,r1)
    var=var_1;
    label=label_1;
    NMSE=(norm(var -  estimation_variance ,'fro')/norm(var,'fro'))^2;
    mdl=fitglm( estimation_variance,label,'Distribution','binomial','Link','logit');
    scores=mdl.Fitted.Probability;
    [X,Y,T,AUC] = perfcurve(label,scores,'1');    
  else
    var=var_2;
    label=label_2;
    NMSE=(norm(var -  estimation_variance ,'fro')/norm(var,'fro'))^2;
    mdl=fitglm( estimation_variance,label,'Distribution','binomial','Link','logit');
    scores=mdl.Fitted.Probability;
    [X,Y,T,AUC] = perfcurve(label,scores,'1');         
  end
    AUC_RESULT=AUC_RESULT+AUC;
    NMSE_RESULT=NMSE_RESULT+NMSE;
end

if iter_result==5000
    iter=iter+1;
    ITER_RESULT_FINAL(total,iteration)=1;
end



AUC_RESULT_FINAL(total,iteration)= AUC_RESULT/120;
NMSE_RESULT_FINAL(total,iteration)=NMSE_RESULT/120;
% ESTIMATE_RESULT(total,iteration)=(norm(x - estimate ,'fro')/norm(x,'fro'))^2;

   end
end
 
%% https://plotly.com/matlab/roc-and-pr-curves/



mean(AUC_RESULT_FINAL,2)
mean(NMSE_RESULT_FINAL,2)
sum(ITER_RESULT_FINAL,2)
% mean(ESTIMATE_RESULT,2)

std(AUC_RESULT_FINAL,0,2)
std(NMSE_RESULT_FINAL,0,2)
