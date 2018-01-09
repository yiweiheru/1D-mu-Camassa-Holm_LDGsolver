close

global type_of_flux
% 1 is Alternating; 2 is Central.
type_of_flux = 1;

n_RK=2;
Tfinal=0.0001;
Ord=3;
ir=4;
cfl=0.1;

Nelm=10*2^(ir-1);
dx=period/Nelm;

dt=dx*cfl;
Tsteps=floor((Tfinal-0.000000001)/dt)+1;

x=-period/2:dx:period/2;
elm_size=Ord+1;
U0 = setInitial(Nelm,elm_size,x,Xexc,uexc);
[ Amat,Dumat,massMat,massMat_inv ] = getAmat( Ord,Nelm,x );
U=U0;
Time=0;
for nt=1:Tsteps
    if nt == Tsteps
        dt = Tfinal - Time;
    end
    U = RKn( x,Ord,Nelm,U,Amat,massMat_inv,Dumat,massMat,n_RK,dt );
    Time=Time+dt;
end
R = Dumat*U;


fig=0;
fig = fig+1;figure(fig)
plot_uh(U,Ord,Nelm,x,-period/2,period/2);
title("U")


% [ U_final,~ ] =getFinalExactSol( Nelm,elm_size,x,Time,Xexc,uexc,rexc ,period,c);
% fig = fig+1;figure(fig)
% plot_uh(U_final,Ord,Nelm,x,-period/2,period/2);
% title("Ufinal")
% 
% M = Amat*U;
% fig = fig+1;figure(fig)
% plot_uh(M,Ord,Nelm,x,-period/2,period/2);
% title("M1")
% % 
% R = Dumat*U;
% % fig = fig+1;figure(fig)
% % plot_uh(R,Ord,Nelm,x,-period/2,period/2);
% % title("R")
% [ Residue3 ] = residue3( x,Nelm,Ord,U,M );
% Q=massMat_inv*Residue3;
% % fig = fig+1;figure(fig)
% % plot_uh(Q,Ord,Nelm,x,-period/2,period/2);
% % title("Q")
% [ Residue4 ] = residue4( x,Nelm,Ord,Q,R,M );
% Mt = massMat_inv*Residue4;
% M = Mt*dt+M;
% 
% hold on
% figure(fig)
% plot_uh(M,Ord,Nelm,x,-period/2,period/2);
% 
% 
% fig = fig+1;figure(fig)
% plot_uh(Mt,Ord,Nelm,x,-period/2,period/2);
% title("Mt")
% 
