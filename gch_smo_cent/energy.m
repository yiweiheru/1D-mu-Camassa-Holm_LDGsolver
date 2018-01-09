
n_RK=4;
Tfinal=10;
ord_num=2;
ir_num=4;
cfl=0.1;

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
        Time=0;
        for nt=1:Tsteps
            if nt == Tsteps
               dt = Tfinal - Time;
            end
            U=RKn( Ord,x,Nelm,U,Amat,massMat_inv,n_RK,dt);
            R = getuh_Der( Ord,x,Nelm,U,massMat_inv );
            Ey(nt) = getEnergy( U,R,Ord,Nelm,x); 
            Time=Time+dt;
        end
    end
end
global orr
hold on
Ey_pt = 0.1*linspace(0,Tfinal,size(Ey,2));
if orr == 1
    plot(Ey_pt,Ey,'g')
else
    plot(Ey_pt,Ey,'k')
end
