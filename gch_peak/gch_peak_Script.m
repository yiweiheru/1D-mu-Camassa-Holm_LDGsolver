clear
close all

Tfinal = 15;
fig_u_at_time  = [0,1,5,10,Tfinal];
fig_Ey_at_time = 0:0.5:(Tfinal);

ord_num = 2;
ir_num = 4;

n_RK  = 4;
cfl = 0.1;
period = 1; 

peak_type = 2;

P0 = zeros(peak_type,1);
Q0 = zeros(peak_type,1);
switch peak_type
    case 1
        P0(1) = 0.1;      Q0(1) = 0;
    case 2
        P0(1) = 0.1;      Q0(1) = 0.2;
        P0(2) = 0.08;     Q0(2) = 0.1;
    case 3
        P0(1) = 0.1;      Q0(1) = 0.2;
        P0(2) = 0.08;     Q0(2) = 0.1;
        P0(3) = 0.12;     Q0(3) = 0.05;
end

Eyc  = zeros(size(fig_Ey_at_time,2),2);
Eyd  = zeros(size(fig_Ey_at_time,2),2);

fig = 0;
nEy = 0;
for Ord = ord_num:ord_num
    for ir = ir_num:ir_num
        
        Nelm = 10*2^(ir-1)+1;
        elm_size = Ord+1;
        
        dx = period/Nelm;
        x = 0:dx:period;
        
        Time = 0;
        dt = cfl * dx;
        Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
        dt_final = Tfinal - (Tsteps-1) * dt;
        %--------------------------------------------------------------------------   
        [ Amat,~,~,~,massMat_inv,~ ] = getAmat( Ord,Nelm,x,'R','L' );
        U0 = setInitial(Nelm,elm_size,x,peak_type,period,P0,Q0);
        R0 = getuh_Der( Ord,x,Nelm,U0,massMat_inv );
        Uc = U0;
        Ud = U0;
        %--------------------------------------------------------------------------
        nEy = nEy+1;
        Eyc(nEy,:) = getEnergy( U0,R0,Ord,Nelm,x);
        Eyd(nEy,:) = Eyc(nEy,:);
        %------------------------------------------------------------------        
        fig = fig + 1;
        figure(fig)
        plot_uh( U0,Ord,Nelm,x ,"exact")
        grid on
        xlabel('x')
        ylabel('u')
        Title_str = strcat('t=',num2str(Time));
        title(Title_str);
        %--------------------------------------------------------------------------
        for nt = 1:Tsteps
            if nt == Tsteps
                dt = dt_final;
            end
            %--------------------------------------------------------------
            Uc = RKn( Ord,x,Nelm,Uc,Amat,massMat_inv,n_RK,dt,Time,'Csv' );
            Ud = RKn( Ord,x,Nelm,Ud,Amat,massMat_inv,n_RK,dt,Time,'Dsp' );
            Time = Time+dt;
            %--------------------------------------------------------------
            for p=1:size(fig_Ey_at_time,2)
                if fig_Ey_at_time(p)-dt/2 < Time && Time <= fig_Ey_at_time(p)+dt/2
                    nEy = nEy+1;
                    Rc = getuh_Der( Ord,x,Nelm,Uc,massMat_inv );
                    Rd = getuh_Der( Ord,x,Nelm,Ud,massMat_inv );
                    Eyc(nEy,:) = getEnergy( Uc,Rc,Ord,Nelm,x);
                    Eyd(nEy,:) = getEnergy( Ud,Rd,Ord,Nelm,x);
                end
            end
            %--------------------------------------------------------------
            for p=1:size(fig_u_at_time,2)
                if fig_u_at_time(p)-dt/2 < Time && Time <= fig_u_at_time(p)+dt/2
                    fig = fig+1;
                    figure(fig)
                    plot_uh( Uc,Ord,Nelm,x ,"numerical" )
                    hold on
                    plot_uh( Ud,Ord,Nelm,x ,"other" )
                    [ Uexc,~,~ ] = multi_pkns_solu( P0,Q0,Nelm,elm_size, x ,period,Time,peak_type );
                    plot_uh( Uexc,Ord,Nelm,x ,"exact" )
                    grid on
                    xlabel('x')
                    ylabel('u')
                    Title_str = strcat('t=',num2str(Time));
                    title(Title_str);
                    legend('LDG(Csv)','LDG(Dsp)','Exact','Location','best')
                end
            end
        end
    end
end
% fig = fig+1;
% figure(fig)
% plot(fig_Ey_at_time,Ey(:,1)','-.','LineWidth',1.5)
% hold on
% plot(fig_Ey_at_time,Ey(:,2)','-','LineWidth',1.5)
% xlabel('t')
% legend('E_0','E_1')
% ylim([0,0.4])


len = size(fig_Ey_at_time,2);
diff_Eyc0 = abs(Eyc(2:len,1)'-Eyc(1,1)*ones(1,len-1));
diff_Eyd0 = abs(Eyd(2:len,1)'-Eyd(1,1)*ones(1,len-1));
diff_Eyc1 = abs(Eyc(2:len,2)'-Eyc(1,2)*ones(1,len-1));
diff_Eyd1 = abs(Eyd(2:len,2)'-Eyd(1,2)*ones(1,len-1));

fig = fig+1;
figure(fig)
semilogy(fig_Ey_at_time(2:len),diff_Eyc0,'-^','LineWidth',1.5)
hold on
semilogy(fig_Ey_at_time(2:len),diff_Eyd0,'-x','LineWidth',1.5)
semilogy(fig_Ey_at_time(2:len),diff_Eyc1,'-v','LineWidth',1.5)
semilogy(fig_Ey_at_time(2:len),diff_Eyd1,'-o','LineWidth',1.5)
grid on

xlabel('t')
legend('E_0(csv)','E_0(dsp)','E_1(csv)','E_1(dsp)', 'Location','NorthEastOutside')
title('|E(t)-E(0)|')

