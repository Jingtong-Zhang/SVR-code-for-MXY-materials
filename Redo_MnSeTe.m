clear 
clc 
%% 载入Feature和Energy并进行归一化,注意，这里的ps1和ps2是归一化的系数，后面不要修改了
load('Feature_and_energy.mat');
Total_test_result=zeros(1008,2);
%% 虽然不知道会不会起到作用，但还是把p_test和t_test直接从Features拿到然后归一化吧
p_test=Features;
t_test=Energy;

p_train=Features(1:100,:);
t_train=Energy(1:100);

%输入样本归一化
[pn_train,ps1] = mapminmax(p_train');
pn_train = pn_train';
pn_test = mapminmax('apply',p_test',ps1);
pn_test = pn_test';
%输出样本归一化
[tn_train,ps2] = mapminmax(t_train');
tn_train = tn_train';
tn_test = mapminmax('apply',t_test',ps2);
tn_test = tn_test';
%注意，这里的归一化，是根据训练集数据进行的归一化，也许每次做的时候都需要重新归一化才行
%% SVR模型创建/训练
% 寻找最佳c参数/g参数——交叉验证方法
% SVM模型有两个非常重要的参数C与gamma。
% 其中 C是惩罚系数，即对误差的宽容度。
% c越高，说明越不能容忍出现误差,容易过拟合。C越小，容易欠拟合。C过大或过小，泛化能力变差
% gamma是选择RBF函数作为kernel后，该函数自带的一个参数。隐含地决定了数据映射到新的特征空间后的分布，
% gamma越大，支持向量越少，gamma值越小，支持向量越多。支持向量的个数影响训练与预测的速度。
[c,g] = meshgrid(-10:0.5:10,-10:0.5:10);
[m,n] = size(c);
cg = zeros(m,n);
eps = 10^(-4);
v = 5;
bestc = 0;
bestg = 0;
error = Inf;
for i = 1:m
    for j = 1:n
        cmd = ['-v ',num2str(v),' -t 2',' -c ',num2str(2^c(i,j)),' -g ',num2str(2^g(i,j) ),' -s 3 -p 0.1'];
        cg(i,j) = svmtrain(tn_train,pn_train,cmd);
        if cg(i,j) < error
            error = cg(i,j);
            bestc = 2^c(i,j);
            bestg = 2^g(i,j);
        end
        if abs(cg(i,j) - error) <= eps && bestc > 2^c(i,j)
            error = cg(i,j);
            bestc = 2^c(i,j);
            bestg = 2^g(i,j);
        end
    end
end
cmd = [' -t 2',' -c ',num2str(bestc),' -g ',num2str(bestg),' -s 3 -p 0.01'];
model = svmtrain(tn_train,pn_train,cmd);
%% SVR仿真预测
[Predict_2,error_2,dec_values_2] = svmpredict(tn_test,pn_test,model);
% 反归一化
predict_2 = mapminmax('reverse',Predict_2,ps2);
Total_test_result(:,1)=predict_2;
%% 计算误差
[len,~]=size(predict_2);
error = t_test - predict_2;
error = error';
MAE1=sum(abs(error./t_test'))/len;
MSE1=error*error'/len;
RMSE1=MSE1^(1/2);
R = corrcoef(t_test,predict_2);
r = R(1,2);
disp(['........支持向量回归误差计算................'])
disp(['平均绝对误差MAE为:',num2str(MAE1)])
disp(['均方误差为MSE:',num2str(MSE1)])
disp(['均方根误差RMSE为:',num2str(RMSE1)])
disp(['决定系数 R^2为:',num2str(r)])
%figure(1)
%plot(1:length(t_test),t_test,'r-*',1:length(t_test),predict_2,'b:o')
%不需要作图，就把这一行注释掉
%grid on
%legend('真实值','预测值')
%xlabel('样本编号')
%ylabel('值')
%string_2 = {'测试集预测结果对比'};
%title(string_2)
%% 排序与筛选
%下面的这个All_re是储存所有的Feature和Energy,储存的是归一化后的
ALL_re=zeros(1008,22);%19是Feature数+4
% SVR能量（归一化后）,DFT能量（归一化后），Features,编号，是否在训练集中
% 把SVR放在前面是因为如果sortrows的排序值是一样的话，接下来会按照第二列来进行排序的
ALL_re(:,3:20)=Features;
%ALL_re(:,1)=t_test;
%ALL_re(:,2)=predict_2; %如果按照SVR排序，这里应该是tn_test，前面的应该是Predict_2
ALL_re(:,2)=t_test;
ALL_re(:,1)=predict_2;

%上面这段改成用实际的值，而归一化的值
% pn_test的顺序不变的话，ALL_re的顺序也不会变，因此每次只需要更新ALL_re(:,1)就可以了
i=1;
while i<=1008
 ALL_re(i,21)=i; %如果材料换了，这里的19也要改一下
 i=i+1;
end
ALL_re(1:100,22)=1;
ALL_re(101:1008,22)=0; %代表前100个在训练集里，后面的不在。
B=sortrows(ALL_re);
% FR_check用来查看前20个的R2值
FR_check1=B(1:20,1);
FR_check2=B(1:20,2);
FR_R = corrcoef(FR_check2,FR_check1);
FR_R2 = FR_R(1,2);
[LLx,LLy]=size(pn_train);
%% 这里还是得重新做一遍归一化

t_train_new=t_train;
p_train_new=p_train;


ADD=10;%这里的ADD代表每次循环要加入多少个新的结构
i=1;
already_add=0;
while i<=1008
    if B(i,22)==0 %这里要修改成为size后的值，因为最后一列是储存‘是否输入TS’的
        t_train_new(LLx+1+already_add)=B(i,2);
        p_train_new(LLx+1+already_add,:)=B(i,3:20); %如果Feature的尺寸改变，这里也要修改
        Current=B(i,21);%这个代表编号
        ALL_re(Current,22)=1;
        already_add=already_add+1;
        if already_add>ADD
            i=999999;
        end
    end 
    i=i+1;
end

%% 新的训练集已经出来了，就是pn_train_new和tn_train_new，接下来要按照新的训练集来做了吧。
How_many_loop=30; %一共循环多少次
loop_now=1;
ALL_R2=zeros(2);
ALL_R2(1,1)=r;
ALL_R2(1,2)=FR_R2;
%这里的r是前面那一步的R2
while loop_now<=How_many_loop
    p_train=p_train_new;
    t_train=t_train_new;
    
    
%输入样本归一化
[pn_train,ps1] = mapminmax(p_train');
pn_train = pn_train';
pn_test = mapminmax('apply',p_test',ps1);
pn_test = pn_test';
%输出样本归一化
[tn_train,ps2] = mapminmax(t_train');
tn_train = tn_train';
tn_test = mapminmax('apply',t_test',ps2);
tn_test = tn_test';
    
    
%% Loop内部SVR模型创建/训练
% 寻找最佳c参数/g参数——交叉验证方法
% SVM模型有两个非常重要的参数C与gamma。
% 其中 C是惩罚系数，即对误差的宽容度。
% c越高，说明越不能容忍出现误差,容易过拟合。C越小，容易欠拟合。C过大或过小，泛化能力变差
% gamma是选择RBF函数作为kernel后，该函数自带的一个参数。隐含地决定了数据映射到新的特征空间后的分布，
% gamma越大，支持向量越少，gamma值越小，支持向量越多。支持向量的个数影响训练与预测的速度。
[c,g] = meshgrid(-10:0.5:10,-10:0.5:10);
[m,n] = size(c);
cg = zeros(m,n);
eps = 10^(-4);
v = 5;
bestc = 0;
bestg = 0;
error = Inf;
for i = 1:m
    for j = 1:n
        cmd = ['-v ',num2str(v),' -t 2',' -c ',num2str(2^c(i,j)),' -g ',num2str(2^g(i,j) ),' -s 3 -p 0.1'];
        cg(i,j) = svmtrain(tn_train,pn_train,cmd);
        if cg(i,j) < error
            error = cg(i,j);
            bestc = 2^c(i,j);
            bestg = 2^g(i,j);
        end
        if abs(cg(i,j) - error) <= eps && bestc > 2^c(i,j)
            error = cg(i,j);
            bestc = 2^c(i,j);
            bestg = 2^g(i,j);
        end
    end
end
cmd = [' -t 2',' -c ',num2str(bestc),' -g ',num2str(bestg),' -s 3 -p 0.01'];
model = svmtrain(tn_train,pn_train,cmd);
%% SVR仿真预测
[Predict_2,error_2,dec_values_2] = svmpredict(tn_test,pn_test,model);

% 反归一化
predict_2 = mapminmax('reverse',Predict_2,ps2);   
Total_test_result(:,loop_now+1)=predict_2;
[len,~]=size(predict_2);
error = t_test - predict_2;
error = error';
MAE1=sum(abs(error./t_test'))/len;
MSE1=error*error'/len;
RMSE1=MSE1^(1/2);
R = corrcoef(t_test,predict_2);
r = R(1,2);
%% Loop内添加新的训练集
    ALL_re(:,1)=predict_2;   %这里的预测和DFT的需要修改下
    B=sortrows(ALL_re);
    FR_check1=B(1:20,1);
    FR_check2=B(1:20,2);
    FR_R = corrcoef(FR_check2,FR_check1);
    FR_R2 = FR_R(1,2);
    t_train_new=t_train;
    p_train_new=p_train;
    
    
j=1;
already_add=0;
[LLx,LLy]=size(pn_train);
while j<=1008
    if B(j,22)==0
        t_train_new(LLx+1+already_add)=B(j,2); %这里需要注意DFT和SVR的顺序
        p_train_new(LLx+1+already_add,:)=B(j,3:20); %如果Feature的尺寸改变，这里也要修改
        Current=B(j,21);%这个代表编号
        ALL_re(Current,22)=1;
        already_add=already_add+1;
        if already_add>ADD
            j=999999;
        end 
    end 
    j=j+1;
end
    loop_now=loop_now+1;
    ALL_R2(loop_now,1)=r;
    ALL_R2(loop_now,2)=FR_R2;
end 
ALL_R2