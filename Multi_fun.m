function ObjV=Multi_fun(X)
    [m,~]=size(X);
    ObjV=zeros(m,1);
    for i=1:m
        ObjV(i,1)=sin(pi*X(i,1)*X(i,2))+sin(X(i,2)^2);
    end
end