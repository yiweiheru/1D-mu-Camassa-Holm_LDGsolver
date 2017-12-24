function [ LInf_err ] = LInf_error( Uh,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c )

elm_size=Ord+1;
[Ue,~]= getFinalExactSol( Nelm,elm_size,x,Time,Xexc,uexc,rexc,period,c );
U_diff=Ue-Uh;
u_diff=uhTransform( Nelm,elm_size,U_diff );

np = 11;
xi = linspace(-1,1,11);

un=zeros(elm_size,np);
for k = 1 : np
    un(:,k) = basis_1d(Ord,xi(k));
end

LInf_err = 0;
for ne = 1:Nelm
    for ik = 1:np
        utemp = abs(u_diff(ne,:)*un(:,ik));
        if utemp >= LInf_err
            LInf_err = utemp;
        end
    end
end

end

