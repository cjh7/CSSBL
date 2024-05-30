clear;
warning('off','all')
L1=50;
L2=50;
N=33;
AUC_RESULT_FINAL=zeros(10,20);
NMSE_RESULT_FINAL=zeros(10,20);
ITER_RESULT_FINAL=zeros(10,20);
ESTIMATE_RESULT=zeros(10,20);
for total=1:10
    disp(total)
   for iteration=1:20
disp(iteration)
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case2\SSBL\ESTIMATE');       
filename = sprintf('estimate_%d_%d.csv', total,iteration );
estimate=csvread(filename);       
       
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case2\SSBL\X');       
filename = sprintf('x_%d_%d.csv', total,iteration );
x=csvread(filename);   
       

r1=(1:50);  
r2=(51:100);

cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case2\SSBL\Z');       
filename = sprintf('z_%d_%d.csv', total,iteration );
z=csvread(filename);  


cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case2\SSBL\GAMMA');       
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance=csvread(filename);
       
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case2\SSBL\FAULT_INDEX');       
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
fault_index=csvread(filename);

cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case2\SSBL\ITERATION')
filename = sprintf('iteration_%d_%d.csv', total,iteration );
iter_result=csvread(filename);

iter=0;

% Label Information
label_1=zeros(N-4,1);
label_1(1)=1;


label_2=zeros(N-4,1);
label_2(1)=1;
label_2(2)=1;
label_2(fault_index(1)-4)=1;
label_2(fault_index(2)-4)=1;

% Variance Information
var_1=ones(N-4,1)*0.01;
var_1(1)=1;


var_2=ones(N-4,1)*0.01;
var_2(1)=1;
var_2(2)=1;
var_2(fault_index(1)-4)=1;
var_2(fault_index(2)-4)=1;

AUC_RESULT=0;
NMSE_RESULT=0;
for i=1:100

  estimation_variance=variance(:,1);
  
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


AUC_RESULT_FINAL(total,iteration)= AUC_RESULT/100;
NMSE_RESULT_FINAL(total,iteration)=NMSE_RESULT/100;
% ESTIMATE_RESULT(total,iteration)=(norm(x - estimate ,'fro')/norm(x,'fro'))^2;

   end
end
 
%% https://plotly.com/matlab/roc-and-pr-curves/



mean(AUC_RESULT_FINAL,2)
mean(NMSE_RESULT_FINAL,2)
sum(ITER_RESULT_FINAL,2)
% mean(ESTIMATE_RESULT,2)
