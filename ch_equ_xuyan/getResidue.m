function [ Residue ] = getResidue( Ord,x,Nelm,U,massMat_inv,Time,period )

elm_size=Ord+1;

uh=uhTransform(Nelm,elm_size,U);

[ Residue1 ] = residue1( x,Nelm,Ord,uh,Time,period );
R=massMat_inv*Residue1;
% R= Limiter( R,Ord,Nelm,x );
rh=uhTransform(Nelm,elm_size,R);

[ Residue2 ] = residue2( x,Nelm,Ord,uh,rh,Time,period );
P=massMat_inv*Residue2;
% P= Limiter( P,Ord,Nelm,x );
ph=uhTransform(Nelm,elm_size,P);

[ Residue3 ] = residue3( x,Nelm,Ord,uh,rh,ph,Time,period );

Residue=Residue3;

end

