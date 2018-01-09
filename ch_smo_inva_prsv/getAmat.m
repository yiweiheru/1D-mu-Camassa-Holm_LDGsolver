function [ Amat,Dumat,massMat_inv ] = getAmat( Ord,Nelm,x )

% 1 represent Left, 2 represent Right

elm_size=Ord+1;

[ Mmat,Pmat,mu_Mmat,M_inv,FLPmat,FLMmat,FRPmat,FRMmat ] = Mat_sys( Ord,x,Nelm );

% compute the inverse of mass matrix
isp=zeros(Nelm*elm_size*elm_size+1,1);
jsp=zeros(Nelm*elm_size*elm_size+1,1);
sdat=zeros(Nelm*elm_size*elm_size+1,1);

it=0;
for ne=1:Nelm
    Jacob=(x(ne+1)-x(ne))/2;
    for i=1:elm_size
        for j=1:elm_size
            it=it+1;
            isp(it)=(ne-1)*elm_size+i;
            jsp(it)=(ne-1)*elm_size+j;
            sdat(it)=M_inv(i,j)/Jacob;
        end
    end
end
it=it+1;
isp(it)=Nelm*elm_size;
jsp(it)=Nelm*elm_size;
sdat(it)=0;

massMat_inv=sparse(isp,jsp,sdat);

%compute the mass matrix & prime matrix
isp1=zeros(Nelm*elm_size*elm_size+1,1);
jsp1=zeros(Nelm*elm_size*elm_size+1,1);
sdat1=zeros(Nelm*elm_size*elm_size+1,1);

isp2=zeros(Nelm*elm_size*elm_size+1,1);
jsp2=zeros(Nelm*elm_size*elm_size+1,1);
sdat2=zeros(Nelm*elm_size*elm_size+1,1);


it=0;jt=0;
for ne=1:Nelm
    for i=1:elm_size
        for j=1:elm_size
            it=it+1;
            isp1(it)=(ne-1)*elm_size+i;
            jsp1(it)=(ne-1)*elm_size+j;
            sdat1(it)=Mmat(i,j,ne);
            
            jt=jt+1;
            isp2(jt)=(ne-1)*elm_size+i;
            jsp2(jt)=(ne-1)*elm_size+j;
            sdat2(jt)=Pmat(i,j,ne);
        end
    end
end

it=it+1;
isp1(it)=Nelm*elm_size;
jsp1(it)=Nelm*elm_size;
sdat1(it)=0;

jt=jt+1;
isp2(jt)=Nelm*elm_size;
jsp2(jt)=Nelm*elm_size;
sdat2(jt)=0;

massMat=sparse(isp1,jsp1,sdat1);
PriMat=sparse(isp2,jsp2,sdat2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global type_of_flux
% 1 is Alternating; 2 is Central.

if type_of_flux == 1
    isp3=zeros(2*Nelm*elm_size*elm_size+1,1);
    jsp3=zeros(2*Nelm*elm_size*elm_size+1,1);
    sdat3=zeros(2*Nelm*elm_size*elm_size+1,1);
    
    isp4=zeros(2*Nelm*elm_size*elm_size+1,1);
    jsp4=zeros(2*Nelm*elm_size*elm_size+1,1);
    sdat4=zeros(2*Nelm*elm_size*elm_size+1,1);
    it=0;jt=0;
    for ne=1:Nelm
        for i=1:elm_size
            for j=1:elm_size
                %uh_hat take plus
                it=it+1;
                isp3(it)  = (ne-1)*elm_size+i;
                jsp3(it)  = (ne-1)*elm_size+j;
                sdat3(it) = FLPmat(i,j,1,ne);
                it=it+1;
                isp3(it)  = (ne-1)*elm_size+i;
                jsp3(it)  = (nei(ne,2)-1)*elm_size+j;
                sdat3(it) = (-1) * FRPmat(i,j,2,ne); %take care of the flux at (j+1/2) by -1
                
                % rh_hat take minus
                % take care of the flux on (j+1/2)by -1
                jt=jt+1;
                isp4(jt)  = (ne-1)*elm_size+i;
                jsp4(jt)  = (nei(ne,1)-1)*elm_size+j;
                sdat4(jt) = FLMmat(i,j,1,ne);
                jt=jt+1;
                isp4(jt)  = (ne-1)*elm_size+i;
                jsp4(jt)  = (ne-1)*elm_size+j;
                sdat4(jt) = (-1) * FRMmat(i,j,2,ne); 
                
            end
        end
    end
    
    
    it=it+1;
    isp3(it)=Nelm*elm_size;
    jsp3(it)=Nelm*elm_size;
    sdat3(it)=0;
    jt=jt+1;
    isp4(jt)=Nelm*elm_size;
    jsp4(jt)=Nelm*elm_size;
    sdat4(jt)=0;
    
    Pu = PriMat+sparse(isp3,jsp3,sdat3);
    Pr = PriMat+sparse(isp4,jsp4,sdat4);
    
    % matrix between U and R, R = Dumat*U;
    Dumat = (-1) * massMat_inv * Pu;
    Drmat = (-1) * massMat_inv * Pr;
    
    % matrix between U and M, M = Amat*U;
    Amat  = massMat_inv*massMat - Drmat*Dumat;
else
    isp3=zeros(2*2*Nelm*elm_size*elm_size+1,1);
    jsp3=zeros(2*2*Nelm*elm_size*elm_size+1,1);
    sdat3=zeros(2*1*Nelm*elm_size*elm_size+1,1);
    
    it=0;
    for ne=1:Nelm
        for i=1:elm_size
            for j=1:elm_size
                %uh_hat take plus
                it=it+1;
                isp3(it)  = (ne-1)*elm_size+i;
                jsp3(it)  = (ne-1)*elm_size+j;
                sdat3(it) = (0.5) * FLPmat(i,j,1,ne);
                it=it+1;
                isp3(it)  = (ne-1)*elm_size+i;
                jsp3(it)  = (nei(ne,1)-1)*elm_size+j;
                sdat3(it) = (0.5) * FLMmat(i,j,1,ne);
                %take care of the flux at (j+1/2) by -1
                it=it+1;
                isp3(it)  = (ne-1)*elm_size+i;
                jsp3(it)  = (nei(ne,2)-1)*elm_size+j;
                sdat3(it) = (-0.5) * FRPmat(i,j,2,ne); 
                it=it+1;
                isp3(it)  = (ne-1)*elm_size+i;
                jsp3(it)  = (ne-1)*elm_size+j;
                sdat3(it) = (-0.5) * FRMmat(i,j,2,ne); 
            
                
            end
        end
    end
    
    
    it=it+1;
    isp3(it)=Nelm*elm_size;
    jsp3(it)=Nelm*elm_size;
    sdat3(it)=0;

    Pu = PriMat+sparse(isp3,jsp3,sdat3);

    % matrix between U and R, R = Dumat*U;
    Dumat = (-1) * massMat_inv * Pu;
    Drmat = Dumat;
    
    % matrix between U and M, M = Amat*U;
    Amat  = massMat_inv*massMat - Drmat*Dumat;
end


end

