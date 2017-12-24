function [ l2error ] = L2err_discrete( UStore,UexcStore,ir_num,ord_num,period )


l2error=zeros(ord_num,ir_num);
for Ord=1:ord_num
    for ir=1:ir_num
        Nelm=10*2^(ir-1);
        dx=period/Nelm;
        x=-period/2:dx:period/2;
        num_pt=20;
        xk=zeros(Nelm,(num_pt+1));
        for ne=1:Nelm
            xk(ne,:)=linspace(x(ne),x(ne+1),num_pt+1);
        end
        elm_size=Ord+1;
        sizeU=elm_size*Nelm;
        U=UStore(1:sizeU,ir,Ord);
        Uexc=UexcStore(1:sizeU,ir,Ord);
        uh=uhTransform(Nelm,Ord+1,U);
        uhexc=uhTransform(Nelm,Ord+1,Uexc);
        evalU=zeros(Nelm,(num_pt+1));
        evalUexc=zeros(Nelm,(num_pt+1));
        for ne=1:Nelm
            for i=1:num_pt+1
                evalU(ne,i)=evalue_uh( Ord,Nelm,x,uh,xk(ne,i) );
                evalUexc(ne,i)=evalue_uh( Ord,Nelm,x,uhexc,xk(ne,i) );
            end
        end
        Du=evalU-evalUexc;
        val=0;
        NP=0;
        for ne=1:Nelm
            for i=1:num_pt+1
%                 if (xk(ne,i)<=-5 && xk(ne,i)>=-25) || (xk(ne,i)>=5 && xk(ne,i)<=25)
                if ~(xk(ne,i)>=-4 && xk(ne,i)<=4)
                    val=val+(Du(ne,i))^2;
                    NP=NP+1;
                end
            end
        end
        l2error(Ord,ir)=sqrt(val/NP);
        
    end
end


end

