Ord = 3;
ir  = 5;
period = 1;
Nelm=10*2^(ir-1);
dx=period/Nelm;
x=-period/2:dx:period/2;
elm_size=Ord+1;
[ Amat,massMat_inv ] = getAmat( Ord,Nelm,x );
