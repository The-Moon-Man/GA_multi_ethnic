%% 智能算法第六例多种群遗传算法（简化修改）
clc
clear
close all
%% 定义遗传算法参数
NIND=40;    	%个体数目
NVAR=2;         %变量的维数
MAXGEN=10;      %最大遗传代数
PRECI=10;       %变量的二进制位数
GGAP=0.9;       %代沟
MP=10;          %多种群的种群数量
Y_max=0;                               	%最优值
px=0.7+(0.9-0.7)*rand(MP,1);            %交叉概率（0.7,0.9）
pm=0.001+(0.05-0.001)*rand(MP,1);      	%变异概率（0.001,0.05）
Scope_individual=[0 0;2 2];     %个体的自变量范围
%区域描述器
FieldD=[rep(PRECI,[1,NVAR]);Scope_individual;rep([1;0;1;1],[1,NVAR])];
Chrom=cell(1,MP);
ObjV=cell(1,MP);
for i=1:MP
    Chrom{i}=crtbp(NIND,NVAR*PRECI);  	%初始种群
end
for i=1:MP
    X=bs2rv(Chrom{i},FieldD);         	%计算初始种群的十进制转换
    ObjV{i}=Multi_fun(X);
end

%% 优化
gen=0;                                  %初始遗传代数
gen_keep=0;                             %初始保持代数
MaxObjV=zeros(MP,1);                    %记录精华种群
MaxChrom=zeros(MP,NVAR*PRECI);          %记录精华种群的编码
FitnV=cell(1,MP);
SelCh=cell(1,MP);
while gen_keep<=MAXGEN
    gen=gen+1;
    for i=1:MP
        FitnV{i}=ranking(-ObjV{i});                   	%分配适应度值
        SelCh{i}=select('sus',Chrom{i},FitnV{i},GGAP);	%选择
        SelCh{i}=recombin('xovsp',SelCh{i},px(i));   	%重组
        SelCh{i}=mut(SelCh{i},pm(i));                 	%变异
        X=bs2rv(SelCh{i},FieldD);                     	%计算初始种群的十进制转换
        ObjVSel=Multi_fun(X);                                               %计算目标函数值
        [Chrom{i},ObjV{i}]=reins(Chrom{i},SelCh{i},1,1,ObjV{i},ObjVSel);   %重新插入
    end
    [Chrom,ObjV]=Multi_immigration(Chrom,ObjV);                             %移民操作
    [MaxObjV,MaxChrom]=Artificial_selection(Chrom,ObjV,MaxObjV,MaxChrom);   %人工选择  
    
    %找出并记录每代最优值
    Y_current(gen)=max(MaxObjV);	%找到精华种群中最优的个体
    if Y_current(gen)>Y_max                                                
        Y_max=Y_current(gen);       %更新最优值
        gen_keep=0;
    else
        gen_keep=gen_keep+1;        %最优值保持次数加1
    end
end
%% 结果展示
%画进化图
figure(2);
plot(1:gen,Y_current);
grid on
xlabel('进化代数')
ylabel('最优解变化')
title('进化过程')
xlim([1,gen])
%最优解
[Y,I]=max(MaxObjV);
X=bs2rv(MaxChrom(I,:),FieldD);
fprintf(['最优解:\nX=',num2str(X),'\n最优值Y=',num2str(Y),'\n']);