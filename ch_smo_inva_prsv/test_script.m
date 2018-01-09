clear 
Nelm=100;
Ord=2;
elm_size=Ord+1;
dx=1/Nelm;
x=-0.5:dx:0.5;

U0=setInitial(Nelm,elm_size,x);
uh=uhTransform(Nelm,elm_size,U0);

[ Amat,massMat_inv ] = getAmat( Ord,Nelm,x );


[ Residue1 ] = residue1( x,Nelm,Ord,uh );
R=massMat_inv*Residue1;
rh=uhTransform(Nelm,elm_size,R);

massMat_inv=full(massMat_inv);

figure(1)
plot_uh( rh,Ord,Nelm,x )