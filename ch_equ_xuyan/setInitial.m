function [ U0 ] = setInitial(Nelm,elm_size,x)

U0=zeros(Nelm*elm_size,1);
for ne=1:Nelm

    for i=1:elm_size
        xtemp=x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num=(ne-1)*elm_size+i;
        if xtemp<=0
            U0(num,1)=0.25*exp(xtemp);
        else
            U0(num,1)=0.25*exp(-1*xtemp);
        end
    end
end

end 
