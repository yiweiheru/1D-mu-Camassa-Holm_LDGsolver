close

n_RK=4;
Tfinal=1;
ord_num=3;
ir_num=5;
cfl=0.1;

L2_ErrorStore=zeros(ord_num,ir_num);
LInf_ErrorStore=zeros(ord_num,ir_num);
muH1_ErrorStore=zeros(ord_num,ir_num);
L2_OrderStore=zeros(ord_num,ir_num-1);
LInf_OrderStore=zeros(ord_num,ir_num-1);
muH1_OrderStore=zeros(ord_num,ir_num-1);

global orr
orr = 1; % flux of f(u):1 is central; 2 is LF.

for Ord=1:ord_num
    for ir=1:ir_num
        
        Nelm=10*2^(ir-1);
        dx=period/Nelm;

        dt=dx*cfl;
        Tsteps=floor((Tfinal-0.000000001)/dt)+1;
        
        x=-period/2:dx:period/2;
        elm_size=Ord+1;
        U0 = setInitial(Nelm,elm_size,x,Xexc,uexc);
        [ Amat,massMat_inv ] = getAmat( Ord,Nelm,x );
        U=U0;
        
        Time=0;
        for nt=1:Tsteps-1
            U=RKn( Ord,x,Nelm,U,Amat,massMat_inv,n_RK,dt);
            Time=Time+dt;
        end
        dtFinal = Tfinal-Time;
        U = RKn( Ord,x,Nelm,U,Amat,massMat_inv,n_RK,dtFinal );
        Time = Time+dtFinal;
        
        R = getuh_Der( Ord,x,Nelm,U,massMat_inv );

        L2_ErrorStore(Ord,ir)= L2_error( U,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c);
        LInf_ErrorStore(Ord,ir)= LInf_error( U,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c);
        muH1_ErrorStore(Ord,ir)= mu_H1_error( U,R,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c );
    
    end
    
    L2_OrderStore(Ord,:)=L2_OrderStore(Ord,:)+ErrorOrder(L2_ErrorStore(Ord,:));
    LInf_OrderStore(Ord,:)=LInf_OrderStore(Ord,:)+ErrorOrder(LInf_ErrorStore(Ord,:));
    muH1_OrderStore(Ord,:)=muH1_OrderStore(Ord,:)+ErrorOrder(muH1_ErrorStore(Ord,:));

end

format short
disp(L2_OrderStore)
disp(LInf_OrderStore)
disp(muH1_OrderStore)

format shortE
disp(L2_ErrorStore)
disp(LInf_ErrorStore)
disp(muH1_ErrorStore)
% 
% figure(1)
% plot_uh(U,Ord,Nelm,x,-period/2,period/2);
% grid on
% 
% figure(2)
% plot_uh(U0,Ord,Nelm,x,-period/2,period/2);
% hold on
% R0 = getuh_Der( Ord,x,Nelm,U0,massMat_inv );
% plot_uh(R0,Ord,Nelm,x,-period/2,period/2);
% % grid on
% [ U_final,R_final ] =getFinalExactSol( Nelm,elm_size,x,Time,Xexc,uexc,rexc ,period,c);
% figure(3)
% % hold on
% plot_uh(R_final,Ord,Nelm,x,-period/2,period/2);