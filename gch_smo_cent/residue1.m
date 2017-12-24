function [ Residue1 ] = residue1( x,Nelm,Ord,uh )

elm_size=Ord+1;

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un=zeros(elm_size,npt_quad);
un_der=zeros(elm_size,npt_quad);
uf=zeros(elm_size,2);

for k = 1 : npt_quad
    un(:,k)      =basis_1d(Ord,qpt(k));
    un_der(:,k)=basisDer_1d(Ord,qpt(k));
end

uf(:,1)=basis_1d(Ord,-1);
uf(:,2)=basis_1d(Ord,1);

Residue1=zeros(Nelm*elm_size,1);
Intergral=zeros(Nelm*elm_size,1);
Flux_L=zeros(Nelm*elm_size,1);
Flux_R=zeros(Nelm*elm_size,1);

%imply the periodic boundary condition to get the neighbourhood.
nei=zeros(Nelm,2);
for ne=1:Nelm
    if ne==1
        nei(ne,1)=Nelm;
        nei(ne,2)=2;
    elseif ne==Nelm
        nei(ne,1)=Nelm-1;
        nei(ne,2)=1;
    else
        nei(ne,1)=ne-1;
        nei(ne,2)=ne+1;
    end
end

for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        u=0;
        for j=1:elm_size
            u=u+un(j,ik)*uh(ne,j);
        end
        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num)+qwt(ik)*un_der(i,ik)*u;
        end
    end  
end

for ne=1:Nelm
    uLp=0;uRp=0;uLm=0;uRm=0;
    for j=1:elm_size
        uLp=uLp+uf(j,1)*uh(ne,j);            %At (j-1/2), uh takes plus
        uLm=uLm+uf(j,2)*uh(nei(ne,1),j);
        uRp=uRp+uf(j,1)*uh(nei(ne,2),j);   %At (j+1/2), uh takes plus
        uRm=uRm+uf(j,2)*uh(ne,j);
    end
    uhat_L=(uLp+uLm)/2;
    uhat_R=(uRp+uRm)/2;
    
    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Flux_L(num)=Flux_L(num)+uf(i,1)*uhat_L;
        Flux_R(num)=Flux_R(num)+uf(i,2)*uhat_R;
    end
end

Residue1=Residue1-(Intergral-Flux_R+Flux_L);

end


