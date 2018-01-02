function [ Unew ] = RKn( x,Ord,Nelm,U,Amat,massMat_inv,Rmat,n_RK,dt )


if n_RK==3
    residue_1=getResidue(x,Ord,Nelm,U,Amat,massMat_inv,Rmat);
    Ut_1=Amat\residue_1;
    U1=U+dt*Ut_1;
    
    residue_2=getResidue(x,Ord,Nelm,U1,Amat,massMat_inv,Rmat);
    Ut_2=Amat\residue_2;
    U2=3/4*U+1/4*U1+(dt/4)*Ut_2;
    
    residue_3=getResidue(x,Ord,Nelm,U2,Amat,massMat_inv,Rmat);
    Ut_3=Amat\residue_3;
    U3=1/3*U+2/3*U2+2/3*dt*Ut_3;
    
    Unew=U3;
end


if n_RK==4
    residue_1=getResidue(x,Ord,Nelm,U,Amat,massMat_inv,Rmat);
    Ut_1=Amat\residue_1;
    U1=U+(dt/2)*Ut_1;
    
    residue_2=getResidue(x,Ord,Nelm,U1,Amat,massMat_inv,Rmat);
    Ut_2=Amat\residue_2;
    U2=U+(dt/2)*Ut_2;
    
    residue_3=getResidue(x,Ord,Nelm,U2,Amat,massMat_inv,Rmat);
    Ut_3=Amat\residue_3;
    U3=U+dt*Ut_3;
    
    residue_4=getResidue(x,Ord,Nelm,U3,Amat,massMat_inv,Rmat);
    Ut_4=Amat\residue_4;
    
    Unew=U+dt/6*(Ut_1+2*Ut_2+2*Ut_3+Ut_4);
end

end
