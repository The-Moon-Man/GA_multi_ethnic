function [Chrom,ObjV]=Multi_immigration(Chrom,ObjV)
%%移民算子
    MP=length(Chrom);
    for i=1:MP
        [~,MaxI]=max(ObjV{i});
        Total_target=i+1;
        if Total_target>MP;Total_target=mod(Total_target,MP);end
        [~,MinI]=min(ObjV{Total_target});
        %%目标种群最劣个体替换成源种群最优个体
        Chrom{Total_target}(MinI,:)=Chrom{i}(MaxI,:);
        ObjV{Total_target}(MinI,:)=ObjV{i}(MaxI,:);
    end
end