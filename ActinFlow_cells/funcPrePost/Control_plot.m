%----------CELL ALONG TIME ----------------------
hfig=figure(4); hold on;

pos = get(hfig,'position');
set(hfig,'position',pos.*[0.1 0.1 2 2])

tt = timedata/60;

%----------PLOT DENSITY  RHO ALONG TIME (PER NODE)----------------------
subplot(2,2,1,'fontsize',16);
xlabel('t [min]'); ylabel('Density [-]'); hold on

    plot(tt,cellDensdata(1,:),'k',tt,cellDensdataM(1,:),'b',tt,cellDensdataAdh(1,:),'r',tt,cellDensdata(end,:),'k--',tt,cellDensdataM(end,:),'b--',tt,cellDensdataAdh(end,:),'r--','LineWidth',2); hold on;
    lgd = legend('Actin fiber', 'Attached Myosin','Adhesion','Location','northeast');
    lgd.FontSize = 8;

    


%----------PLOT VELOCITY ALONG TIME (PER NODE)---------------------

subplot(2,2,2,'fontsize',16);
xlabel('t [min]'); ylabel('Velocity [\mum/s]'); 
yyaxis left
ax=gca;
ax.YColor='k';
asint = find(cellinfo{1}.V_cell>0,1);

yyaxis right
ax=gca;
ax.YColor='b';

semilogy(tt,vpol(1,:),'m--',tt,velocitydata(1,:),'b-.',tt,wdata(1,:),'c:','LineWidth',2); 

lgd = legend('Velocity-Polymerisation', 'Velocity-Lab frame','Velocity-Cell frame','Location','northeast');
lgd.FontSize = 8;

%----------PLOT TENSION ALONG TIME (PER NODE)---------------------
subplot(2,2,3,'fontsize',16);
xlabel('t [min]'); ylabel('Tension [kPa]'); hold on
yyaxis left
ax=gca;
ax.YColor='k';
plot(tt,memTension,'k','LineWidth',2)
yyaxis right
ax=gca;
ax.YColor='b';
ylabel('Cell radius [\mum]');
plot(tt,cell_length_info/2.0,'b','linewidth',2)

%----------PLOT CELL VELOCITY ALONG TIME ---------------------
subplot(2,2,4,'fontsize',16);
xlabel('t [min]'); ylabel('Cell Velocity [\mum/s]'); hold on
plot(tt,cellinfo{1}.V_cell(:),'k','LineWidth',2); 


%----------CELL ALONG SPACE ----------------------
hfig2=figure(5); hold on;

pos = get(hfig2,'position');
set(hfig2,'position',pos.*[0.1 0.1 3.2 0.9])

xx = cellinfo{1}.meshparam.X;
%----------PLOT DENSITY ALONG SPACE (PER NODE)----------------------
subplot(1,3,1,'fontsize',16);
xlabel('x [\mum]'); ylabel('Density [-]'); hold on

    plot(xx,cellDensdata(:,end),'k',xx,cellDensdataM(:,end),'b',xx,cellDensdataAdh(:,end),'r', 'LineWidth',2);
    lgd = legend('Actin fiber', 'Attached Myosin','Adhesion','Location','northeast');
    lgd.FontSize = 8;

%----------PLOT VELOCITY ALONG SPACE (PER NODE)---------------------
subplot(1,3,2,'fontsize',16);
xlabel('x [\mum]'); ylabel('Velocity [\mum/s]'); hold on
plot(xx,velocitydata(:,end),'k',xx,wdata(:,end),'b','LineWidth',2); 
lgd = legend('Velocity-Lab frame','Velocity-Cell frame','Location','northeast');
lgd.FontSize = 8;
%----------PLOT TENSION ALONG SPACE (PER NODE)---------------------
subplot(1,3,3,'fontsize',16);
xlabel('x [\mum]'); ylabel('Tension [s]'); hold on
plot(xx,tension_actin(:,end),'k','LineWidth',2); 

%-----------PLOT ADHESION DENSITY along time and space -----------------------------------

hfig3 = figure(6); hold on;
pos = get(hfig3,'position');
set(hfig3,'position',pos.*[0.1 0.1 3.2 0.9])

xx = cellinfo{1}.meshparam.X;


subplot(1,2,1,'fontsize',16);
xlabel('x [\mum]'); ylabel('Density [-]'); hold on

    plot(xx,cellDensdataAdh(:,end),'r','LineWidth',2);
    lgd = legend('Adhesion','Location','northeast');
    lgd.FontSize = 8;

subplot(1,2,2, 'fontsize',16);
xlabel('t [min]'); ylabel('Density [-]'); hold on
    plot(tt,cellDensdataAdh(1,:),'k',tt,cellDensdataAdh(end,:),'r--','LineWidth',2);



%----------------SAVING FIGURES AND PARAMETERS------------------------------

Cname = string('Cell_migration');

saveas(hfig,char(Cname + string('_time')),'svg');
saveas(hfig2,char(Cname + string('_space')),'svg');

%save(char(Cname+string('.mat')),'timedata', 'cellxcoordsdata', 'cellDensdata', 'cellDensdataM', 'RSDensdata', 'cellDensdataG', 'cellDensdatam', 'wdata', 'parameter', 'cellinfo','cell_length_info','velocitydata','vpol','tension_actin','chargeDensdata','memTension')
save(char(Cname),'timedata', 'cellxcoordsdata', 'cellDensdata', 'cellDensdataM', 'cellDensdataAdh', 'wdata', 'parameter', 'cellinfo','cell_length_info','velocitydata','vpol','tension_actin','memTension')
