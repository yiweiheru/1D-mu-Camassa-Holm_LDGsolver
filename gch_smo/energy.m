close all
clear

run getExc.m

n_RK=4;
Tfinal=10;
ord_num=2;
ir_num=4;
cfl=0.1;
fig_Ey_at_time = 0:0.1:(Tfinal);
Ey  = zeros(size(fig_Ey_at_time,2),2);

% 1:center; 2:LF
global orr
orr = 2;

fig = 0;
nEy = 0;
for Ord=ord_num:ord_num
    for ir=ir_num:ir_num        
        Nelm=10*2^(ir-1);
        dx=period/Nelm;

        dt=dx*cfl;
        Tsteps=floor((Tfinal-0.000000001)/dt)+1;
        
        x=-period/2:dx:period/2;
        elm_size=Ord+1;
        
        [ Amat,massMat_inv ] = getAmat( Ord,Nelm,x );
        U0 = setInitial(Nelm,elm_size,x,Xexc,uexc);
        U=U0;

        nEy = nEy+1;
        R = getuh_Der( Ord,x,Nelm,U,massMat_inv );
        Ey(nEy,:) = getEnergy( U,R,Ord,Nelm,x);
%         
        Time=0;
        for nt=1:Tsteps
            if nt == Tsteps
               dt = Tfinal - Time;
            end
            Time=Time+dt;
            U=RKn( Ord,x,Nelm,U,Amat,massMat_inv,n_RK,dt);
            for p=1:size(fig_Ey_at_time,2)
                if fig_Ey_at_time(p)-dt/2 < Time && Time <= fig_Ey_at_time(p)+dt/2
                    nEy = nEy+1;
                    R = getuh_Der( Ord,x,Nelm,U,massMat_inv );
                    Ey(nEy,:) = getEnergy( U,R,Ord,Nelm,x);
                end
            end
            
        end

        
  
    end
end

len = size(fig_Ey_at_time,2);
fig = fig+1;
figure(fig)
diff_Ey0 = abs(Ey(2:len,1)'-Ey(1,1)*ones(1,len-1));
semilogy(fig_Ey_at_time(2:len),diff_Ey0,'-.','LineWidth',1.5)
hold on
diff_Ey1 = abs(Ey(2:len,2)'-Ey(1,2)*ones(1,len-1));
semilogy(fig_Ey_at_time(2:len),diff_Ey1,'-','LineWidth',1.5)
xlabel('t')
legend('E_0','E_1')
title('semilogy of E(t)-E(0)')
grid on


% hold on
% Ey_pt = 0.1*linspace(0,Tfinal,size(Ey,2));
% if orr == 1
%     plot(Ey_pt,Ey(),'r')
% else
%     plot(Ey_pt,Ey,'b')
% end
% 
