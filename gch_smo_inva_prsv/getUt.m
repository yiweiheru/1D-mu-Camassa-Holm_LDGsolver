
function [ Ut ] = getUt(x,Ord,Nelm,U,Amat,Dumat,massMat_inv)

M = Amat*U;
R = Dumat*U;

[ Residue3 ] = residue3( x,Nelm,Ord,U,M );
Q=massMat_inv*Residue3;

[ Residue4 ] = residue4( x,Nelm,Ord,Q,R,M );
Mt = massMat_inv*Residue4;

Ut = Amat\Mt;