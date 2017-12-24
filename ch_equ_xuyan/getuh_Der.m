function [ R ] = getuh_Der( Ord,x,Nelm,U,massMat_inv,Time )

elm_size=Ord+1;
uh=uhTransform(Nelm,elm_size,U);

[ Residue1 ] = residue1( x,Nelm,Ord,uh,Time );
R=massMat_inv*Residue1;

end

