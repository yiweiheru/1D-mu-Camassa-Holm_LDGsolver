function [ Unew ] = RKn_A( Ord,x,Nelm,U,Amat,massMat_inv,n_RK,dt,Time )

switch n_RK
    
    case 1
        
        residue_1=AgetResidue(Ord,x,Nelm,U,massMat_inv,Time);
        Ut_1=Amat\residue_1;
        U1=U+dt*Ut_1;
        
        Unew =U1;
        
    case 2
        
        residue_1=AgetResidue(Ord,x,Nelm,U,massMat_inv,Time);
        Ut_1=Amat\residue_1;
        U1=U+dt*Ut_1;
        Time1=Time+dt;
        
        residue_2=AgetResidue(Ord,x,Nelm,U1,massMat_inv,Time1);
        Ut_2=Amat\residue_2;
        U2=U+1/2*U1+(dt/2)*Ut_2;

        Unew = U2;
        
    case 3
   
        residue_1=AgetResidue(Ord,x,Nelm,U,massMat_inv,Time);
        Ut_1=Amat\residue_1;
        U1=U+dt*Ut_1;
        Time1=Time+dt;
        
        residue_2=AgetResidue(Ord,x,Nelm,U1,massMat_inv,Time1);
        Ut_2=Amat\residue_2;
        U2=3/4*U+1/4*U1+(dt/4)*Ut_2;
        Time2=Time+1/2*dt;
        
        residue_3=AgetResidue(Ord,x,Nelm,U2,massMat_inv,Time2);
        Ut_3=Amat\residue_3;
        U3=1/3*U+2/3*U2+2/3*dt*Ut_3;
        
        Unew=U3;
        
    case 4
        
        residue_1=AgetResidue(Ord,x,Nelm,U,massMat_inv,Time);
        Ut_1=Amat\residue_1;
        U1=U+(dt/2)*Ut_1;
        Time1=Time+dt/2;
        
        residue_2=AgetResidue(Ord,x,Nelm,U1,massMat_inv,Time1);
        Ut_2=Amat\residue_2;
        U2=U+(dt/2)*Ut_2;
        Time2=Time+dt/2;
        
        residue_3=AgetResidue(Ord,x,Nelm,U2,massMat_inv,Time2);
        Ut_3=Amat\residue_3;
        U3=U+dt*Ut_3;
        Time3=Time+dt;
        
        residue_4=AgetResidue(Ord,x,Nelm,U3,massMat_inv,Time3);
        Ut_4=Amat\residue_4;
        
        Unew=U+dt/6*(Ut_1+2*Ut_2+2*Ut_3+Ut_4);
        
end

end


