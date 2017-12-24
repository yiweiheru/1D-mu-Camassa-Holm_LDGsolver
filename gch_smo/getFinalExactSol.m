function [ U_final,R_final ] = getFinalExactSol( Nelm,elm_size,x,Time,Xexc,uexc,rexc,period,c )

U_final=zeros(Nelm*elm_size,1);
R_final=zeros(Nelm*elm_size,1);
for ne=1:Nelm
    for i=1:elm_size
        xtemp=x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
        num=(ne-1)*elm_size+i;
        if xtemp<=c*Time-period/2
        U_final(num,1)=interp1(Xexc,uexc,xtemp+period-c*Time);
        R_final(num,1)=interp1(Xexc,rexc,xtemp+period-c*Time);
        else
        U_final(num,1)=interp1(Xexc,uexc,xtemp-c*Time);
        R_final(num,1)=interp1(Xexc,rexc,xtemp-c*Time);
        end
    end
end

end

