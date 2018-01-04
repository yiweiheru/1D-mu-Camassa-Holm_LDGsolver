function [ Residue3 ] = residue3( x,Nelm,Ord,uh,mh )

elm_size=Ord+1;

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un     = zeros(elm_size,npt_quad);
un_der = zeros(elm_size,npt_quad);
for k = 1 : npt_quad
    un(:,k)     = basis_1d(Ord,qpt(k));
    un_der(:,k) = basisDer_1d(Ord,qpt(k));
end

Residue3=zeros(Nelm*elm_size,1);
Intergral=zeros(Nelm*elm_size,1);

for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        u=0;m=0;
        for j=1:elm_size
            u=u+un(j,ik)*uh(ne,j);
            m=m+un(j,ik)*mh(ne,j);
        end
        
        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num)-Jaco*qwt(ik)*un(i,ik)*(u*m);
        end
    end
end

Residue3=Residue3-Intergral;

end


