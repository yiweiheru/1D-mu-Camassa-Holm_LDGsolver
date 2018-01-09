function [ Residue4 ] = residue4( x,Nelm,Ord,Q,R,M )

elm_size=Ord+1;
qh=uhTransform(Nelm,elm_size,Q);
mh=uhTransform(Nelm,elm_size,M);
rh=uhTransform(Nelm,elm_size,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un=zeros(elm_size,npt_quad);
un_der=zeros(elm_size,npt_quad);

for k = 1 : npt_quad
    un(:,k)      =basis_1d(Ord,qpt(k));
    un_der(:,k)=basisDer_1d(Ord,qpt(k));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ufL = basis_1d(Ord,-1); %left boundary of elment
ufR = basis_1d(Ord,1);  %right boundary of element

Residue4=zeros(Nelm*elm_size,1);
Intergral=zeros(Nelm*elm_size,1);
Flux_L=zeros(Nelm*elm_size,1);
Flux_R=zeros(Nelm*elm_size,1);

for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    % compute the integral part
    for ik=1:npt_quad
        q=0;m=0;r=0;
        for j=1:elm_size
            q=q+un(j,ik)*qh(ne,j);
            m=m+un(j,ik)*mh(ne,j);
            r=r+un(j,ik)*rh(ne,j);
        end

        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num) - qwt(ik)*un_der(i,ik)*q ...
                + Jaco*qwt(ik)*un(i,ik)*(r*m);
        end
    end
    % add the boundary part
    qLm=0;qLp=0;qRm=0;qRp=0;
    
    qLm = qLm + qh(nei(ne,1),:)*ufR;
    qLp = qLp + qh(ne,:)*ufL;
    qRm = qRm + qh(ne,:)*ufR;
    qRp = qRp + qh(nei(ne,2),:)*ufL;


    qhat_L = (1/2)*(qLm+qLp);
    qhat_R = (1/2)*(qRm+qRp);

    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Flux_L(num)=Flux_L(num)+ufL(i)*qhat_L;
        Flux_R(num)=Flux_R(num)+ufR(i)*qhat_R;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Residue4=Residue4-(Intergral+Flux_R-Flux_L);

end

