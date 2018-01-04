function [ Residue4 ] = getResidue(x,Ord,Nelm,U,Amat,massMat_inv,Rmat)

elm_size=Ord+1;

M = Amat*U;
R = Rmat*U;
uh=uhTransform(Nelm,elm_size,U);
mh=uhTransform(Nelm,elm_size,M);
rh=uhTransform(Nelm,elm_size,R);

[ Residue3 ] = residue3( x,Nelm,Ord,uh,mh );
Q=massMat_inv*Residue3;
qh=uhTransform(Nelm,elm_size,Q);

[ Residue4 ] = residue4( x,Nelm,Ord,qh,rh,mh );


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
end

