
clearvars -except Xexc uexc rexc period c mu_u0
n_RK=4;
Tfinal=1;
ord_num=3;
ir_num=3;
cfl=0.1;


L2_ErrorStore = [ ];
muH1_ErrorStore=[ ];
T_point = [ ];



for Ord=ord_num:ord_num
    for ir=ir_num:ir_num
        
        Nelm=10*2^(ir-1);
        dx=period/Nelm;

        dt=dx*cfl;
        Tsteps=floor((Tfinal-0.000000001)/dt)+1;
        
        x=-period/2:dx:period/2;
        elm_size=Ord+1;
        U0 = setInitial(Nelm,elm_size,x,Xexc,uexc);
        [ Amat,massMat_inv ] = getAmat( Ord,Nelm,x );
        U=U0;
        
        Time = 0;
        itr = 0;
        for nt=1:Tsteps-1
            U=RKn( Ord,x,Nelm,U,Amat,massMat_inv,n_RK,dt);
            Time=Time+dt;
            if mod(Time,0.01) < dt
                itr = itr +1;
                R = getuh_Der( Ord,x,Nelm,U,massMat_inv );
                L2_ErrorStore(itr)= L2_error( U,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c);
                muH1_ErrorStore(itr)= mu_H1_error( U,R,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c );
                T_point(itr) = Time;
            end
        end

        dtFinal = Tfinal-Time;
        U = RKn( Ord,x,Nelm,U,Amat,massMat_inv,n_RK,dtFinal );
        Time = Time+dtFinal;
%         if mod(Time,0.01) < dt
%             itr = itr +1;
%             R = getuh_Der( Ord,x,Nelm,U,massMat_inv );
%              L2_ErrorStore(itr)= L2_error( U,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c);
%             muH1_ErrorStore(itr)= mu_H1_error( U,R,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c );
%             T_point(itr) = Time;
%         end
        
    end
end

hold on
global orr
if orr ==1
%     plot(T_point,muH1_ErrorStore,'r')
    plot(T_point,L2_ErrorStore,'r')
elseif orr == 2
%     plot(T_point,muH1_ErrorStore,'b')
    plot(T_point,L2_ErrorStore,'b')
    xlabel("t")
%     ylabel("H_\mu^1 error")
    ylabel("L^2 error")
	legend("Center","L-F")
end
