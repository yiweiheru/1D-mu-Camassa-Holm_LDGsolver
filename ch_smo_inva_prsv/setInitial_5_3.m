function [ U0 ] = setInitial_5_3(Nelm,elm_size,x)

U0=zeros(Nelm*elm_size,1);

c=1;
a=30;
x0=-5;

for ne=1:Nelm

    for i=1:elm_size
        xtemp=x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num=(ne-1)*elm_size+i;
        if abs(xtemp-x0)<=a/2
            U0(num,1)=c/cosh(a/2)*cosh(xtemp-x0);
        else
            U0(num,1)=c/cosh(a/2)*cosh(a-(xtemp-x0));
        end
    end
end

end 
