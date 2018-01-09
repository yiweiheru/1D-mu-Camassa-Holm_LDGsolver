function [ Unew ] = RKn( x,Ord,Nelm,U,Amat,massMat_inv,Dumat,massMat,n_RK,dt )

if n_RK==1
    [ Ut_1 ] = getUt(x,Ord,Nelm,U,Amat,Dumat,massMat_inv);
    Unew=U+dt*Ut_1; 
end
if n_RK==2
    [ Ut_1 ] = getUt(x,Ord,Nelm,U,Amat,Dumat,massMat_inv);
    U1=U+dt*Ut_1; 
        
    [ Ut_2 ] = getUt(x,Ord,Nelm,U1,Amat,Dumat,massMat_inv);
    
    U2=U+1/2*U1+(dt/2)*Ut_2;

    Unew = U2;
end

if n_RK==3
    [ Ut_1 ] = getUt(x,Ord,Nelm,U,Amat,Dumat,massMat_inv);
    U1=U+dt*Ut_1;
    
    [ Ut_2 ] = getUt(x,Ord,Nelm,U1,Amat,Dumat,massMat_inv);
    U2=3/4*U+1/4*U1+(dt/4)*Ut_2;
    
    [ Ut_3 ] = getUt(x,Ord,Nelm,U2,Amat,Dumat,massMat_inv);
    U3=1/3*U+2/3*U2+2/3*dt*Ut_3;
    
    Unew=U3;
end
if n_RK==4
    [ Ut_1 ] = getUt(x,Ord,Nelm,U,Amat,Dumat,massMat_inv);
    U1=U+(dt/2)*Ut_1;
    
    [ Ut_2 ] = getUt(x,Ord,Nelm,U1,Amat,Dumat,massMat_inv);
    U2=U+(dt/2)*Ut_2;
    
    [ Ut_3 ] = getUt(x,Ord,Nelm,U2,Amat,Dumat,massMat_inv);
    U3=U+dt*Ut_3;
    
    [ Ut_4 ] = getUt(x,Ord,Nelm,U3,Amat,Dumat,massMat_inv);
    
    Unew=U+dt/6*(Ut_1+2*Ut_2+2*Ut_3+Ut_4);
end

end

