function [ Residue ] = AgetResidue( Ord,x,Nelm,U,massMat_inv,Time )

elm_size=Ord+1;

uh=uhTransform(Nelm,elm_size,U);

[ Residue1 ] = Aresidue1( x,Nelm,Ord,uh,Time );
R=massMat_inv*Residue1;
rh=uhTransform(Nelm,elm_size,R);

[ Residue2 ] = Aresidue2( x,Nelm,Ord,uh,rh,Time );
P=massMat_inv*Residue2;
ph=uhTransform(Nelm,elm_size,P);

[ Residue3 ] = Aresidue3( x,Nelm,Ord,uh,rh,ph,Time );

Residue=Residue3;

end


