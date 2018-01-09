function [ U_final,R_final ] = getFinalExactSol( Nelm,elm_size,x,Time )

U_final=zeros(Nelm*elm_size,1);
R_final=zeros(Nelm*elm_size,1);
for ne=1:Nelm
    for i=1:elm_size
        xtemp=x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num=(ne-1)*elm_size+i;
        if xtemp<=Time*0.25
            U_final(num,1)=0.25*exp(xtemp-0.25*Time);
            R_final(num,1)=0.25*exp(xtemp-0.25*Time);
        else
            U_final(num,1)=0.25*exp(-1*(xtemp-0.25*Time));
            R_final(num,1)=-0.25*exp(-1*(xtemp-0.25*Time));
        end
    end
end

end

