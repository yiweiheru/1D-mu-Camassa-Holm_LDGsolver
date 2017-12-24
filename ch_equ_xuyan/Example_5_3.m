clear
close
n_RK=3;
Tfinal=0.1;
ord_num=2;
ir_num=6;
period=30;

Ord=5;
ir=6;
Nelm=10*2^(ir-1);
dx=period/Nelm;
x=0:dx:period;
cfl=0.1;
dt=cfl*dx;
Tsteps=floor((Tfinal-0.00000000001)/dt)+1;

elm_size=Ord+1;
U0=setInitial_5_3(Nelm,elm_size,x);
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

figure(1)
% plot_uh( U,Ord,Nelm,x,-period/2,period/2 )
