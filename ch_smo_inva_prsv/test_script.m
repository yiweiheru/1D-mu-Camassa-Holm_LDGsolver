close

n_RK=4;
Tfinal=0.1;
Ord=1;
ir=4;
cfl=0.1;

Nelm=10*2^(ir-1);
dx=period/Nelm;

dt=dx*cfl;
Tsteps=floor((Tfinal-0.000000001)/dt)+1;

x=-period/2:dx:period/2;
elm_size=Ord+1;
U0 = setInitial(Nelm,elm_size,x,Xexc,uexc);
[ Amat,massMat_inv,Rmat,massMat ] = getAmat( Ord,Nelm,x );
U=U0;

% uh=uhTransform(Nelm,elm_size,U);
% [ Residue1 ] = residue1( x,Nelm,Ord,uh);
% R=massMat_inv*Residue1;
% rh=uhTransform(Nelm,elm_size,R);
% 
% [ Residue2 ] = residue2( x,Nelm,Ord,uh,rh);
% M=massMat_inv*Residue2;
% mh=uhTransform(Nelm,elm_size,M);
% 
% [ Residue3 ] = residue3( x,Nelm,Ord,uh,mh );
% Q=massMat_inv*Residue3;
% qh=uhTransform(Nelm,elm_size,Q);
% 
% [ Residue4 ] = residue4( x,Nelm,Ord,qh,rh,mh );

M = Amat*U;
U = Amat\M;
% R = Rmat*U;
% uh=uhTransform(Nelm,elm_size,U);
% mh=uhTransform(Nelm,elm_size,M);
% rh=uhTransform(Nelm,elm_size,R);

% [ Residue3 ] = residue3( x,Nelm,Ord,uh,mh );
% Q=massMat_inv*Residue3;
% qh=uhTransform(Nelm,elm_size,Q);

% [ Residue4 ] = residue4( x,Nelm,Ord,qh,rh,mh );
% plot_uh(M,Ord,Nelm,x,-period/2,period/2);
% hold on
% Ut = Amat\Residue4;
% U = U + Ut*dt;
plot_uh(M,Ord,Nelm,x,-period/2,period/2);
hold on
plot_uh(U,Ord,Nelm,x,-period/2,period/2);