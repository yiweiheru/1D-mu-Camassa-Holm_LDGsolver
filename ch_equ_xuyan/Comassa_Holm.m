clear
close
n_RK=3;
Tfinal=10;
ord_num=2;
ir_num=6;
period=50;

UStore=zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);
UexcStore=zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);

for Ord=1:ord_num
    for ir=1:ir_num
%         clearvars Nelm dx x elm_size U0 Amat massMat  R Time dt Tsteps dtFinal
        Nelm=10*2^(ir-1);
        dx=period/Nelm;
        x=-period/2:dx:period/2;
        cfl=0.1;
        dt=min(cfl*dx,Tfinal/10);
%         dt=Tfinal/(10*2^(ir-1));
        Tsteps=floor((Tfinal-0.00000000001)/dt)+1;
        
        elm_size=Ord+1;
        U0=setInitial(Nelm,elm_size,x);
        [ Amat,massMat_inv,A1 ] = getAmat( Ord,Nelm,x );
                
        Time=0;
        U=U0;
        for nt=1:Tsteps-1
            U=RKn( Ord,x,Nelm,U,Amat,massMat_inv,A1,n_RK,dt,Time ,period);
            Time=Time+dt;
        end
        dtFinal=Tfinal-Time;
        U=RKn( Ord,x,Nelm,U,Amat,massMat_inv,A1,n_RK,dtFinal,Time ,period);
        Time=Time+dtFinal;
        
        lengthU=elm_size*Nelm;
        UStore(1:lengthU,ir,Ord)=U;
%         R=getuh_Der( Ord,x,Nelm,U,massMat_inv,Time );
        [ Uexc,~ ] = getFinalExactSol( Nelm,elm_size,x,Time );
        UexcStore(1:lengthU,ir,Ord)=Uexc;
    end
end
% [ L2err ] = L2_error( Uh,Time,Ord,Nelm,x );
% [ L2error ] = L2err_discrete( UStore,UexcStore,ir_num,ord_num,period );
% [ LIerror ] = LIerr_discrete( UStore,UexcStore,ir_num,ord_num,period );
% for Ord=1:ord_num
%     disp(ErrorOrder(L2error(Ord,:)));
% %     disp(ErrorOrder(L2err(Ord,:)));
% end

figure(1)
plot_uh( U,Ord,Nelm,x,-period/2,period/2 )
[ U_final,R_final ] = getFinalExactSol( Nelm,elm_size,x,Time );
figure(2)
plot_uh( U_final,Ord,Nelm,x,-25,25 )
