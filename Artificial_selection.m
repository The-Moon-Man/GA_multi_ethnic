function [MaxObjV,MaxChrom]=Artificial_selection(Chrom,ObjV,MaxObjV,MaxChrom)
%%人工选择算子
    MP=length(Chrom);
    for i=1:MP
        [MaxY,MaxI]=max(ObjV{i});
        if MaxY>MaxObjV(i)
            MaxObjV(i)=MaxY;
            MaxChrom(i,:)=Chrom{i}(MaxI,:);
        end
    end
end