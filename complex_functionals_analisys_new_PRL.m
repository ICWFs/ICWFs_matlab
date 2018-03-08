clear all
close all

ind_f = 1;
renormalize = true;

xe_axis = load('eix_e.txt');
xn_axis = load('eix_n.txt');
pn_axis = load('eix_pn.txt');
dx_e = abs(xe_axis(2)-xe_axis(1));
dx_n = abs(xn_axis(2)-xn_axis(1));
dp_n = abs(pn_axis(2)-pn_axis(1));

time = load('time.txt');
time = time/1E-15;
t_final = length(time);
Dim_ele = length(xe_axis);
Dim_nuc = length(xn_axis);

time = time(1:t_final);



%% BOPES & PES
PES = load('PES.txt');
% BOPES = load('BOPES.txt');
% % BOEIGST_GR = load('BO_eigst_gr.txt');
% % BOEIGST_EX1 = load('BO_eigst_ex1.txt');
% % BOEIGST_EX2 = load('BO_eigst_ex2.txt');
% % BOEIGST_GR = reshape(BOEIGST_GR,Dim_nuc,Dim_ele);
% % BOEIGST_EX1 = reshape(BOEIGST_EX1,Dim_nuc,Dim_ele);
% % BOEIGST_EX2 = reshape(BOEIGST_EX2,Dim_nuc,Dim_ele);
% 
% % figure(ind_f)
% % plot(xn_axis,BOPES(:,1),'-r',xn_axis,BOPES(:,2),'-b',xn_axis,BOPES(:,3),'-g')
% % ylim([-0.31 0.3])
% % xlim([-8 8])
% % ylabel('Population','FontSize',30)
% % xlabel('time(fs)','FontSize',30)
% % set(gca,'FontSize',30)
% % ind_f = ind_f + 1;

figure(ind_f);
mesh(xe_axis,xn_axis,PES)
ind_f = ind_f + 1;



%% TRAJECTORIES vs TIME:
str = strcat('xe_IPW',int2str(1),'.txt');
energy_e = fopen(str,'r');
xe_exact_time = fread(energy_e,'double');
fclose(energy_e);
xe_exact_aux = xe_exact_time;
N_traj = size(xe_exact_aux,1);

xe_exact = zeros(t_final,N_traj);
xe_IPW = zeros(t_final,N_traj);
xe_hermitian = zeros(t_final,N_traj);
ve_exact = zeros(t_final,N_traj);
ve_IPW = zeros(t_final,N_traj);
ve_hermitian = zeros(t_final,N_traj);
xn_exact = zeros(t_final,N_traj);
xn_IPW = zeros(t_final,N_traj);
xn_hermitian = zeros(t_final,N_traj);
vn_exact = zeros(t_final,N_traj);
vn_IPW = zeros(t_final,N_traj);
vn_hermitian = zeros(t_final,N_traj);
for t=1:t_final
    str = strcat('xe_exact',int2str(t),'.txt');
    energy_e = fopen(str,'r');
    xe_exact_time = fread(energy_e,'double');
    fclose(energy_e);
    xe_exact(t,:) = xe_exact_time;
    
    str = strcat('xe_IPW',int2str(t),'.txt');
    xe_all_cond = fopen(str,'r');
    xe_cond_time = fread(xe_all_cond,'double');
    fclose(xe_all_cond);
    xe_IPW(t,:) = xe_cond_time;   
    
    str = strcat('xe_hermitian',int2str(t),'.txt');
    xe_all_cond = fopen(str,'r');
    xe_cond_time = fread(xe_all_cond,'double');
    fclose(xe_all_cond);
    xe_hermitian(t,:) = xe_cond_time;       
 
    
    
    
    str = strcat('xn_exact',int2str(t),'.txt');
    xn_all_exact = fopen(str,'r');
    xn_exact_time = fread(xn_all_exact,'double');
    fclose(xn_all_exact);
    xn_exact(t,:) = xn_exact_time;

    str = strcat('xn_IPW',int2str(t),'.txt');
    xn_all_cond = fopen(str,'r');
    xn_cond_time = fread(xn_all_cond,'double');
    fclose(xn_all_cond);
    xn_IPW(t,:) = xn_cond_time;
    
    str = strcat('xn_hermitian',int2str(t),'.txt');
    xe_all_cond = fopen(str,'r');
    xe_cond_time = fread(xe_all_cond,'double');
    fclose(xe_all_cond);
    xn_hermitian(t,:) = xe_cond_time;  
    
    
    
    
    
    str = strcat('ve_exact',int2str(t),'.txt');
    ve_all_exact = fopen(str,'r');
    ve_exact_time = fread(ve_all_exact,'double');
    fclose(ve_all_exact);
    ve_exact(t,:) = ve_exact_time;

    str = strcat('ve_IPW',int2str(t),'.txt');
    ve_all_cond = fopen(str,'r');
    ve_cond_time = fread(ve_all_cond,'double');
    fclose(ve_all_cond);
    ve_IPW(t,:) = ve_cond_time;
    
    str = strcat('ve_hermitian',int2str(t),'.txt');
    xe_all_cond = fopen(str,'r');
    xe_cond_time = fread(xe_all_cond,'double');
    fclose(xe_all_cond);
    ve_hermitian(t,:) = xe_cond_time;      
    

    
    

    str = strcat('vn_exact',int2str(t),'.txt');
    vn_all_exact = fopen(str,'r');
    vn_exact_time = fread(vn_all_exact,'double');
    fclose(vn_all_exact);
    vn_exact(t,:) = vn_exact_time;

    str = strcat('vn_IPW',int2str(t),'.txt');
    vn_all_cond = fopen(str,'r');
    vn_cond_time = fread(vn_all_cond,'double');
    fclose(vn_all_cond);
    vn_IPW(t,:) = vn_cond_time;
    
    str = strcat('vn_hermitian',int2str(t),'.txt');
    xe_all_cond = fopen(str,'r');
    xe_cond_time = fread(xe_all_cond,'double');
    fclose(xe_all_cond);
    vn_hermitian(t,:) = xe_cond_time;      

end
interv = 1;

figure(ind_f);
hold on
h1 = plot(time,xe_exact(:,1:interv:end),'-k','Linewidth',4);
h2 = plot(time,xe_IPW(:,1:interv:end),'-b','Linewidth',2);
h3 = plot(time,xe_hermitian(:,1:interv:end),'-r','Linewidth',2);
legend([h1(1);h2(1);h3(1)],'exact','IPW','hermitian');
title('Electronic Trajectories')
xlabel('time(fs)','FontSize',40)
ylabel('r(a_0)','FontSize',40)
set(gca,'FontSize',40)
ind_f = ind_f + 1;

figure(ind_f);
hold on
h1 = plot(time,xn_exact(:,1:interv:end),'-k','Linewidth',4);
h2 = plot(time,xn_IPW(:,1:interv:end),'-b','Linewidth',2);
h3 = plot(time,xn_hermitian(:,1:interv:end),'-r','Linewidth',2);
legend([h1(1);h2(1);h3(1)],'exact','IPW','hermitian');
title('Nuclear Trajectories')
xlabel('time(fs)','FontSize',40)
ylabel('R(a_0)','FontSize',40)
set(gca,'FontSize',40)
ind_f = ind_f + 1;

figure(ind_f);
hold on
h1 = plot(time,ve_exact(:,1:interv:end),'-k','MarkerSize',10);
h2 = plot(time,ve_IPW(:,1:interv:end),'-b','MarkerSize',10);
h3 = plot(time,ve_hermitian(:,1:interv:end),'-r','MarkerSize',10);
legend([h1(1);h2(1);h3(1)],'exact','IPW','hermitian');
title('Bohmian Trajectories Projected on the Electronic Subspace')
xlabel('time(fs)')
ylabel('ve(a_0)')
ind_f = ind_f + 1;

figure(ind_f);
hold on
h1 = plot(time,vn_exact(:,1:interv:end),'-k','MarkerSize',10);
h2 = plot(time,vn_IPW(:,1:interv:end),'-b','MarkerSize',10);
h3 = plot(time,vn_hermitian(:,1:interv:end),'-r','MarkerSize',10);
legend([h1(1);h2(1);h3(1)],'exact','IPW','hermitian');
title('Bohmian Trajectories Projected on the Nuclear Subspace')
xlabel('time(fs)')
ylabel('vn(a_0)')
ind_f = ind_f + 1;

figure(ind_f)
hold on
h1 = plot(xe_exact(:,1:interv:end),xn_exact(:,1:interv:end),'-k');
h2 = plot(xe_IPW(:,1:interv:end),xn_IPW(:,1:interv:end),'-b');
h3 = plot(xe_hermitian(:,1:interv:end),xn_hermitian(:,1:interv:end),'-r');
legend([h1(1);h2(1);h3(1)],'exact','IPW','hermitian');
contour(xe_axis,xn_axis,PES,50)
ind_f = ind_f + 1;


pause

%% EXACT WAVEFUNCTION EVOLUTION
[xe_grid,xn_grid] = meshgrid(xe_axis,xn_axis);
xe_mean = zeros(t_final,1);
xn_mean = zeros(t_final,1);
red_e_exact = zeros(t_final,Dim_ele);
red_n_exact = zeros(t_final,Dim_nuc);
gr_comp_exact = zeros(Dim_nuc,t_final);
ex1_comp_exact = zeros(Dim_nuc,t_final);
decoh_exact = zeros(t_final,1);
dist_P_exact = zeros(Dim_nuc,t_final);

for t=1:t_final%t_final-20
    str = strcat('exact_wavepacket',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    full_wavepacket = fread(wavepacket_file,'double');
    fclose(wavepacket_file);
    full_wavepacket = reshape(full_wavepacket,Dim_nuc,Dim_ele);
    full_wavepacket = full_wavepacket/(sum(full_wavepacket(:))*dx_e*dx_n);
    
    str = strcat('gr_comp_exact',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    gr_comp_exact(:,t) = fread(wavepacket_file,'double');
    fclose(wavepacket_file);    

    str = strcat('ex1_comp_exact',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    ex1_comp_exact(:,t) = fread(wavepacket_file,'double');
    fclose(wavepacket_file);    
    
%     str = strcat('decoh_exact',int2str(t),'.txt');
%     wavepacket_file = fopen(str,'r');
%     decoh_exact(t) = fread(wavepacket_file,'double');
%     fclose(wavepacket_file);   
   
    str = strcat('dist_P',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    dist_P_exact(:,t) = fread(wavepacket_file,'double');
    fclose(wavepacket_file);       
    
    xe_mean(t) = sum(sum(full_wavepacket.*xe_grid))*dx_e*dx_n;
    xn_mean(t) = sum(sum(full_wavepacket.*xn_grid))*dx_e*dx_n;   
    
    red_e_exact(t,:) = sum(full_wavepacket,1)*dx_n;
    red_n_exact(t,:) = sum(full_wavepacket,2)*dx_e;       

    decoh_exact(t) = gr_comp_exact(:,t).'*(ex1_comp_exact(:,t))*dx_n;    %.*red_n_exact(t,:).'
    
%     figure(ind_f);
%     pcolor(xe_axis,xn_axis,full_wavepacket); 
%     shading interp
%     ylabel('R(a_0)')%,'fontsize',30)    
%     xlabel('r(a_0)')%,'fontsize',30)
%     str = num2str(time(t));
%     title(strcat('time: ',str,'fs'))%,'FontSize',40)    
%     %xlim([-20 20])
%     

end
ind_f = ind_f + 2;
% % v = VideoWriter('exact_dynamics.mp4');
% open(v)
% writeVideo(v,M)
% close(v)
% clear M

%% REDUCED DENSITIES AND NUCLEAR MOMENTUM DISTRIBUTION
red_n_IPW = zeros(t_final,Dim_nuc);
red_e_IPW = zeros(t_final,Dim_ele);
red_n_hermitian = zeros(t_final,Dim_nuc);
red_e_hermitian = zeros(t_final,Dim_ele);
red_n_IPW_recon = zeros(t_final,Dim_nuc);
red_e_IPW_recon = zeros(t_final,Dim_ele);

mean_R_IPW = zeros(t_final,1);
mean_R_recon = zeros(t_final,1);
mean_R_hermitian = zeros(t_final,1);
mean_r_IPW = zeros(t_final,1);
mean_r_recon = zeros(t_final,1);
mean_r_hermitian = zeros(t_final,1);

gr_comp_IPW = zeros(Dim_nuc,t_final);
ex1_comp_IPW = zeros(Dim_nuc,t_final);
gr_comp_hermitian = zeros(Dim_nuc,t_final);
ex1_comp_hermitian = zeros(Dim_nuc,t_final);
gr_comp_e_hermitian = zeros(Dim_nuc,t_final);
ex1_comp_e_hermitian = zeros(Dim_nuc,t_final);

decoh_IPW = zeros(t_final,1);
decoh_hermitian = zeros(t_final,1);
decoh_e_hermitian = zeros(t_final,1);
%decoh_n_hermitian = zeros(t_final,1);

momentum_n_IPW = zeros(Dim_nuc,t_final);
momentum_n_hermitian = zeros(Dim_nuc,t_final);

index = 0;
for t=1:t_final
 
    % REDUCED NUCLEAR DENSITY        
    str = strcat('red_n_IPW',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    red_n_IPW_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);    
    
    str = strcat('red_n_IPW_recon',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    red_n_IPW_recon_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);  
    
    str = strcat('red_n_hermitian',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    red_n_hermitian_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);      
    
    red_n_hermitian(t,:) = red_n_hermitian_aux/(sum(red_n_hermitian_aux)*dx_n);
    red_n_IPW(t,:) = red_n_IPW_aux/(sum(red_n_IPW_aux)*dx_n);    
    red_n_IPW_recon(t,:) = red_n_IPW_recon_aux/(sum(red_n_IPW_recon_aux)*dx_n);        

    mean_R_IPW(t) = red_n_IPW(t,:)*xn_axis*dx_n;
    mean_R_recon(t) = red_n_IPW_recon(t,:)*xn_axis*dx_n;
    mean_R_hermitian(t) = red_n_hermitian(t,:)*xn_axis*dx_n;
    
%    if (mod(t,5) == 0 || t == 1)
%        figure(ind_f);
%        plot(xn_axis,red_n_exact(t,:),'-k','linewidth',4) 
%        hold on       
%        plot(xn_axis,red_n_IPW(t,:),'-b','markersize',7)      
%        plot(xn_axis,red_n_IPW_recon(t,:),'ob','markersize',7)
%        plot(xn_axis,red_n_hermitian(t,:),'-r','markersize',7)
%        legend('exact','IPW','recon','hermitian')
%        hold off
%        str = num2str(time(t));
%        title(strcat('time: ',str,'fs'))%,'FontSize',40)
%        xlabel('R(a_0)')%,'FontSize',40)
%        ylabel('reduced n dens.')%,'FontSize',40)
% %       set(gca,'FontSize',40)
%        hold off
%        index = index + 1;
%    end
    
    
    % REDUCED ELECTRONIC DENSITY                
    str = strcat('red_e_IPW',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    red_e_IPW_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);        
        
    str = strcat('red_e_IPW_recon',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    red_e_IPW_recon_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file); 
    
    str = strcat('red_e_hermitian',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    red_e_hermitian_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);     

    
    red_e_hermitian(t,:) = red_e_hermitian_aux/(sum(red_e_hermitian_aux)*dx_e);        
    red_e_IPW(t,:) = red_e_IPW_aux/(sum(red_e_IPW_aux)*dx_e);
    red_e_IPW_recon(t,:) = red_e_IPW_recon_aux/(sum(red_e_IPW_recon_aux)*dx_e);    

    mean_r_IPW(t) = red_e_IPW(t,:)*xe_axis*dx_e;
    mean_r_recon(t) = red_e_IPW_recon(t,:)*xe_axis*dx_e;
    mean_r_hermitian(t) = red_e_hermitian(t,:)*xe_axis*dx_e;

    
%     if (mod(t,5) == 0 || t == 1)
%         figure(ind_f+1);
%         title('REDUCED ELECTRON')
%         plot(xe_axis,red_e_exact(t,:),'-k','linewidth',4)
%         hold on           
%         plot(xe_axis,red_e_IPW(t,:),'-b','markersize',7)     
%         plot(xe_axis,red_e_IPW_recon(t,:),'ob','markersize',7)
%         plot(xe_axis,red_e_hermitian(t,:),'-r','markersize',7)        
%         legend('exact','IPW','recon','hermitian')
%         hold off
%         str = num2str(time(t));
%         title(strcat('time: ',str,'fs'))%,'FontSize',40)
%         xlabel('r(a_0)')%,'FontSize',40)
%         ylabel('reduced e density')%,'FontSize',40)
%         pause(0.3)
%     end
    
    
    str = strcat('gr_comp_IPW',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    gr_comp_IPW_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file); 

    str = strcat('ex1_comp_IPW',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    ex1_comp_IPW_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);    

    gr_comp_IPW(:,t) = gr_comp_IPW_aux./((sum(gr_comp_IPW_aux,1)+sum(ex1_comp_IPW_aux,1))*dx_n);    
    ex1_comp_IPW(:,t) = ex1_comp_IPW_aux./((sum(gr_comp_IPW_aux,1)+sum(ex1_comp_IPW_aux,1))*dx_n);    
    
    decoh_IPW(t) = gr_comp_IPW(:,t).'*(ex1_comp_IPW(:,t))*dx_n;  %.*red_n_IPW(t,:).'
    
    
    str = strcat('comp_gr_n_hermitian',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    gr_comp_hermitian_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);    

    str = strcat('comp_ex1_n_hermitian',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    ex1_comp_hermitian_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);    
    
    gr_comp_hermitian(:,t) = gr_comp_hermitian_aux./((sum(gr_comp_hermitian_aux,1)+sum(ex1_comp_hermitian_aux,1))*dx_n);    
    ex1_comp_hermitian(:,t) = ex1_comp_hermitian_aux./((sum(gr_comp_hermitian_aux,1)+sum(ex1_comp_hermitian_aux,1))*dx_n);  
    
    decoh_hermitian(t) = gr_comp_hermitian(:,t).'*(ex1_comp_hermitian(:,t))*dx_n;   %.*red_n_hermitian(t,:).'


    
    str = strcat('comp_gr_e_hermitian',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    gr_comp_e_hermitian_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);    

    str = strcat('comp_ex1_e_hermitian',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    ex1_comp_e_hermitian_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);    
    
    gr_comp_e_hermitian(:,t) = gr_comp_e_hermitian_aux./((sum(gr_comp_e_hermitian_aux,1)+sum(ex1_comp_e_hermitian_aux,1))*dx_n);    
    ex1_comp_e_hermitian(:,t) = ex1_comp_e_hermitian_aux./((sum(gr_comp_e_hermitian_aux,1)+sum(ex1_comp_e_hermitian_aux,1))*dx_n);  

    decoh_e_hermitian(t) = gr_comp_e_hermitian(:,t).'*(ex1_comp_e_hermitian(:,t))*dx_n;   %.*red_n_hermitian(t,:).'
    
    
    
    
    str = strcat('momentum_n_IPW',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    momentum_n_IPW_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);        
    
    momentum_n_IPW(:,t) = momentum_n_IPW_aux/(sum(momentum_n_IPW_aux)*dp_n);
    
    str = strcat('momentum_n_hermitian',int2str(t),'.txt');
    wavepacket_file = fopen(str,'r');
    momentum_n_hermitian_aux = fread(wavepacket_file,'double');
    fclose(wavepacket_file);  
    
    momentum_n_hermitian(:,t) = momentum_n_hermitian_aux/(sum(momentum_n_hermitian_aux)*dp_n);

    
%    if (mod(t,5) == 0 || t == 1)
%        figure(ind_f+2);
%        plot(xn_axis,gr_comp_exact(:,t),'-r','linewidth',4) 
%        hold on       
%        plot(xn_axis,ex1_comp_exact(:,t),'-b','linewidth',4) 
%        plot(xn_axis,gr_comp_IPW(:,t),'-or','linewidth',4) 
%        plot(xn_axis,ex1_comp_IPW(:,t),'-ob','linewidth',4) 
%        plot(xn_axis,gr_comp_hermitian(:,t),'.r','linewidth',4) 
%        plot(xn_axis,ex1_comp_hermitian(:,t),'.b','linewidth',4) 
%        plot(xn_axis,gr_comp_e_hermitian(:,t),':r','linewidth',4) 
%        plot(xn_axis,ex1_comp_e_hermitian(:,t),':b','linewidth',4) 
%        
% %        legend('exact','IPW','recon','hermitian')
%        hold off
%        str = num2str(time(t));
%        title(strcat('time: ',str,'fs'))%,'FontSize',40)
%        xlabel('R(a_0)')%,'FontSize',40)
%        ylabel('adiabatic components')%,'FontSize',40)
% %       set(gca,'FontSize',40)
%        hold off
%        index = index + 1;
%    end    
   
   if (mod(t,5) == 0 || t == 1)
       figure(ind_f+3);
       plot(pn_axis,dist_P_exact(:,t),'-k','linewidth',4) 
       hold on       
       plot(pn_axis,momentum_n_IPW(:,t),'-ob','linewidth',4) 
       plot(pn_axis,momentum_n_hermitian(:,t),'-r','linewidth',4) 
       legend('exact','IPW','hermitian')
       hold off
       str = num2str(time(t));
       title(strcat('time: ',str,'fs'))%,'FontSize',40)
       xlabel('P(a_0)')%,'FontSize',40)
       ylabel('nuclear momentum distribution.')%,'FontSize',40)
%       set(gca,'FontSize',40)
       hold off
       index = index + 1;
   end      
     pause(0.1)
end
ind_f = ind_f + 4;





%% MEAN ELECTRONIC AND NUCLEAR POSITIONS
figure(ind_f)
hold on
plot(time,xe_mean,'-k','linewidth',4)
plot(time,mean_r_IPW,'-b')%,'linewidth',4)
plot(time,mean_r_recon,'ob')%,'linewidth',4)
plot(time,mean_r_hermitian,'-r')%,'linewidth',4)
legend('exact','IPW','recon','hermitian')
% xlim([0 30])
ylabel('mean e-position')%,'FontSize',40)
xlabel('time(fs)')%,'FontSize',40)
%set(gca,'FontSize',40)
box on
ind_f = ind_f + 1;



figure(ind_f)
hold on
plot(time,xn_mean,'-k','linewidth',4)
plot(time,mean_R_IPW,'-b')%,'linewidth',4)
plot(time,mean_R_recon,'ob')%,'linewidth',4)
plot(time,mean_R_hermitian,'-r')%,'linewidth',4)
legend('exact','IPW','recon','hermitian')
% xlim([0 30])
ylabel('mean n-position')%,'FontSize',40)
xlabel('time(fs)')%,'FontSize',40)
%set(gca,'FontSize',40)
box on
ind_f = ind_f + 1;



figure(ind_f)
plot(time,decoh_exact,'-k','linewidth',5)
hold on
plot(time,decoh_IPW,'ob','markersize',15,'linewidth',2)
% plot(time,decoh_hermitian,'-r',time,decoh_e_hermitian,'--r')
xlabel('time (fs)','fontsize',45)
ylabel('decoherence','fontsize',45)
set(gca,'FontSize',40)
xlim([0 32])
ylim([0 0.15])
box on
set(gca,'linewidth',6)
legend('exact','interacting CWF')
ind_f = ind_f + 1;

figure(ind_f)
plot(time,sum(gr_comp_exact,1)*dx_n,'-r',time,sum(ex1_comp_exact,1)*dx_n,'-b')
hold on
plot(time,sum(gr_comp_IPW,1)*dx_n,'-or',time,sum(ex1_comp_IPW,1)*dx_n,'-ob')
plot(time,sum(gr_comp_hermitian,1)*dx_n,'-.r',time,sum(ex1_comp_hermitian,1)*dx_n,'-.b')
plot(time,sum(gr_comp_e_hermitian,1)*dx_n,':r',time,sum(ex1_comp_e_hermitian,1)*dx_n,':b')

%% ENERGY
% energy = zeros(t_final,1);
% for t=1:t_final
%     str = strcat('total_energy',int2str(t),'.txt');
%     prob_gr_exact = fopen(str,'r');
%     prob_gr_time = fread(prob_gr_exact,'double');
%     fclose(prob_gr_exact);
%     energy(t) = prob_gr_time;
%     
%     str = strcat('potential_energy',int2str(t),'.txt');
%     prob_gr_exact = fopen(str,'r');
%     prob_gr_time = fread(prob_gr_exact,'double');
%     fclose(prob_gr_exact);
%     p_energy(t) = prob_gr_time; 
%     
%     str = strcat('kinetic_energy_e',int2str(t),'.txt');
%     prob_gr_exact = fopen(str,'r');
%     prob_gr_time = fread(prob_gr_exact,'double');
%     fclose(prob_gr_exact);
%     k_energy_e(t) = prob_gr_time;     
%     
%     str = strcat('e_q_potential_energy',int2str(t),'.txt');
%     prob_gr_exact = fopen(str,'r');
%     prob_gr_time = fread(prob_gr_exact,'double');
%     fclose(prob_gr_exact);
%     e_q_p_energy(t) = prob_gr_time; 
%     
%     
%     
% 
%     str = strcat('nuclear_k_energy_cond',int2str(t),'.txt');
%     prob_gr_exact = fopen(str,'r');
%     prob_gr_time = fread(prob_gr_exact,'double');
%     fclose(prob_gr_exact);
%     nuclear_k_energy_cond(t) = prob_gr_time;
% 
%     str = strcat('nuclear_p_energy_cond',int2str(t),'.txt');
%     prob_gr_exact = fopen(str,'r');
%     prob_gr_time = fread(prob_gr_exact,'double');
%     fclose(prob_gr_exact);
%     nuclear_p_energy_cond(t) = prob_gr_time;
%     
%     str = strcat('electronic_k_energy_cond',int2str(t),'.txt');
%     prob_gr_exact = fopen(str,'r');
%     prob_gr_time = fread(prob_gr_exact,'double');
%     fclose(prob_gr_exact);
%     electronic_k_energy_cond(t) = prob_gr_time;    
% 
%     str = strcat('electronic_p_energy_cond',int2str(t),'.txt');
%     prob_gr_exact = fopen(str,'r');
%     prob_gr_time = fread(prob_gr_exact,'double');
%     fclose(prob_gr_exact);
%     electronic_p_energy_cond(t) = prob_gr_time;    
% 
% end
% figure(ind_f)
% hold on
% 
% plot(time,energy,'-k','linewidth',4)
% plot(time,k_energy_e.','--k','linewidth',4)
% plot(time,energy-p_energy.'-k_energy_e.','-.k','linewidth',4)
% plot(time,p_energy.','--k','linewidth',4)
% plot(time,e_q_p_energy,'*k','linewidth',4)
% 
% 
% plot(time,electronic_k_energy_cond+0.5*electronic_p_energy_cond + ...
%           nuclear_k_energy_cond+0.5*nuclear_p_energy_cond,'-b','linewidth',2)
% plot(time,electronic_k_energy_cond,'--b','linewidth',2)
% plot(time,nuclear_k_energy_cond,'-.b','linewidth',2)
% plot(time,electronic_p_energy_cond,'--b','linewidth',2)
% plot(time,nuclear_p_energy_cond,'-.b','linewidth',2)
% % plot(time,electronic_k_energy_cond+nuclear_k_energy_cond,'-.b','linewidth',2)
% % plot(time,0.5*electronic_p_energy_cond+0.5*nuclear_p_energy_cond,'--b','linewidth',2)
% 
% ind_f = ind_f + 1;
% xlabel('time (fs)','fontsize',36)
% ylabel('energy (Ha)','fontsize',36)
% set(gca,'FontSize',25)
% % legend('tot energy exact','k energy e exact','k energy n exact','p energy exact','q pot e energy exact',...
% %        'tot energy ehre','k energy e ehre','k energy n ehre','p energy ehre',...
% %        'tot energy cond','k energy e cond','k energy n cond','p energy e cond','p energy n cond',...
% %        'tot energy qc','k energy e qc','k energy n qc','p energy e qc','p energy n qc')
% 
% pause