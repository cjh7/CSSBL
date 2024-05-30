clear;
warning('off','all')
L1=50;
L2=50;
N=33;
AUC_RESULT_FINAL=zeros(10,20);
NMSE_RESULT_FINAL=zeros(10,20);
ESTIMATE_RESULT=zeros(10,20);
for total=1:5
    disp(total)
   for iteration=1:20
disp(iteration)

cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case1\LEE\ESTIMATE');       
filename = sprintf('estimate_%d_%d.csv', total,iteration );
estimate=csvread(filename);       
       
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case1\LEE\X');       
filename = sprintf('x_%d_%d.csv', total,iteration );
x=csvread(filename);   
       
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case1\LEE\GAMMA');       
filename = sprintf('variance_%d_%d.csv', total,iteration );
variance=csvread(filename);
       
cd('C:\Users\User\Desktop\새 폴더\Code\Result\Real\Case1\LEE\FAULT_INDEX');       
filename = sprintf('FAULT_INDEX_%d_%d.csv', total,iteration );
fault_index=csvread(filename);

r1=(1:50);  
r2=(51:100);
% Label Information
label_1=zeros(N,1);
label_1(1)=1;
label_1(2)=1;
label_1(3)=1;
label_1(4)=1;
label_1(5)=1;
label_1(6)=1;

label_2=zeros(N,1);
label_2(1)=1;
label_2(2)=1;
label_2(3)=1;
label_2(4)=1;
label_2(5)=1;
label_2(6)=1;
label_2(7)=1;
label_2(8)=1;
label_2(9)=1;
label_2(fault_index(1))=1;


% Variance Information
var_1=ones(N,1)*0.01;
var_1(1)=1;
var_1(2)=1;
var_1(3)=1;
var_1(4)=1;
var_1(5)=1;
var_1(6)=1;

var_2=ones(N,1)*0.01;
var_2(1)=1;
var_2(2)=1;
var_2(3)=1;
var_2(4)=1;
var_2(5)=1;
var_2(6)=1;
var_2(7)=1;
var_2(8)=1;
var_2(9)=1;
var_2(fault_index(1))=1;


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


AUC_RESULT_FINAL(total,iteration)= AUC_RESULT/100;
NMSE_RESULT_FINAL(total,iteration)=NMSE_RESULT/100;
% ESTIMATE_RESULT(total,iteration)=(norm(x - estimate ,'fro')/norm(x,'fro'))^2;

   end
end

%% https://plotly.com/matlab/roc-and-pr-curves/



mean(AUC_RESULT_FINAL,2)
mean(NMSE_RESULT_FINAL,2)
% mean(ESTIMATE_RESULT,2)

std(AUC_RESULT_FINAL,0,2)
std(NMSE_RESULT_FINAL,0,2)


% 
% cd('C:\Users\User\Desktop\학위논문심사');
% write.csvfilename = sprintf('kaveh_real_lee.csv' );
% csvwrite(filename,NMSE_RESULT_FINAL);