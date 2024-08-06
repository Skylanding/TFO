clear all
set(0,'defaultAxesFontSize',20,'DefaultLineLineWidth', 2)
    sensitive_4 = zeros(1,20);
    attenuation_t = zeros(1,20);
    attenuation_M = zeros(1,20);
    SNR_4 = zeros(1,20);
    SD = zeros(1,20);
for i=1:50
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
    %
    % In this example, we simulate a 4-layer brain model using MCXLAB.
    % We will investigate the differences between the solutions with and 
    % witout boundary reflections (both external and internal) and show
    % you how to display and analyze the resulting data.
    %
    % This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear cfg;

    %% preparing the input data
    % set seed to make the simulation repeatible
    cfg.seed=hex2dec('623F9A9E'); 
    
    precision = 1;      %0.1cm
    diminmm = 100;      %10cm
    cfg.nphoton=1e7;

    %number of pixels
    dim = diminmm * precision;

    % define a 4 layer structure
    cfg.vol=ones(dim,dim,60);
    cfg.vol(:,:,3:end)=2;
    cfg.vol=uint8(cfg.vol);
    % cfg.unitinmm=1/precision;  %this statement specifies the resultion in mm

    % define the source position
    cfg.srcpos=[10*precision,50*precision,1*precision]; %a 1 by 3 vector, the position of the source in grid unit
    cfg.srcdir=[0 0 1];           % a 1 by 3 vector, specifying the incident vector; if srcdir
                                  % contains a 4th element, it specifies the focal length of
                                  % the source (only valid for focuable src, such as planar, disk,
                                  % fourier, gaussian, pattern, slit, etc); if the focal length
                                  % is nan, all photons will be launched isotropically regardless
                                  % of the srcdir direction.
    cfg.unitinmm=1/precision;  %this statement specifies the resultion in mm

    % use the brain optical properties defined at
    % http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/CollinsAtlasMesh
    % format: [mua(1/mm) mus(1/mm) g n]

    cfg.prop=[0 0 1 1            % medium 0: the environment
       0.02 1.3  0.8  1.4     % medium 1: Maternal abdominal wall
       0.01 1.2    0.8  1.4];     % medium 2: Maternal uterus

    % time-domain simulation parameters
    cfg.tstart=0;
    cfg.tend=100e-9;
    cfg.tstep=0.5e-10;

    % GPU thread configuration
    cfg.autopilot=1;
    cfg.gpuid=1;

    %% running simulation without boundary reflection
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    cfg.isreflect=0; % disable reflection at exterior boundary
    tic;
    f1=mcxlab(cfg);
    toc;

    %% running simulation with boundary reflection enabled
    fprintf('running simulation ... this takes about 50 seconds on a GTX 470\n');
    cfg.isreflect=1; % enable reflection at exterior boundary
    cfg.isrefint=1;  % enable reflection at interior boundary too
    cfg.issavedet=1; % enable recording partial pathlength of detected photons
    cfg.detpos=[(10+4*i)*precision 50*precision 1*precision 1*precision];
    tic;
    [f2,det2]=mcxlab(cfg);
    toc;

    %% plot the results
%     figure;
%     contourf(log10(squeeze(sum(f1.data(:,51*precision,:,:),4))'),1:0.5:8);
%     hold on
%     plot([0 diminmm],[3 3],'--r');
%     title('flux with no reflection');
%     set(gca,'clim',[1 8]);
%     colorbar
% 
%     figure;
%     contourf(log10(squeeze(sum(f2.data(:,51*precision,:,:),4))'),1:0.5:8);
%     hold on
%     plot([0 diminmm],[3 3],'--r');
%     title('flux with reflection at boundaries');
%     set(gca,'clim',[1 8]);
%     colorbar
% 
%     figure;
%     xi=1e9*((cfg.tstart:cfg.tstep:cfg.tend-cfg.tstep)+0.5*cfg.tstep);
%     semilogy(xi,squeeze(f2.data(51,41,1,:)));
%     xlabel('time (ns)')
%     ylabel('flux (1/mm^2/s)')
%     set(gca,'yscale','log');
%     title('time point spread functions (TPSF)');
% % 
%     figure;
%     [hs1,c1]=hist(det2.ppath(find(det2.ppath(:,1)),1),200);
%     hold on;
%     [hs2,c2]=hist(det2.ppath(find(det2.ppath(:,2)),2),200);
%     [hs3,c3]=hist(det2.ppath(find(det2.ppath(:,3)),3),200);
%     hold on;
%     [hs4,c4]=hist(det2.ppath(find(det2.ppath(:,4)),4),200);
%     bar(c1,hs1,'edgecolor','none','facecolor','r');
%     bar(c2,hs2,'edgecolor','none','facecolor','b');
%     bar(c3,hs3,'edgecolor','none','facecolor','g');
%     bar(c4,hs4,'edgecolor','none','facecolor','k');
%     legend('tissue 1','tissue 2','tissue 3','tissue 3');
%     xlabel('partial pathlength (mm)');
%     title(sprintf('detected %d photons',size(det2.ppath,1)));

    SD(i)= 4*i*0.1;

%     sensitive_4(i)=length(find(det2.ppath(:,4)))/size(det2.ppath,1);
    attenuation_t(i)=size(det2.ppath,1)/3e5/1.4;
%     attenuation_M(i)=(size(det2.ppath,1)-length(find(det2.ppath(:,4))))/cfg.nphoton;
%     SNR_4(i)=length(find(det2.ppath(:,4)))/cfg.nphoton;
end

% figure;
% plot(SD,sensitive_4,'-ob','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b');
% xlabel('Source-detector distance (cm)');
% ylabel('Fetal signal sensitivity');
% xticks(0:1:8)
% grid on;
% grid minor

figure;
semilogy(SD,attenuation_t,'-or','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');
xlabel('Source-detector distance (cm)');
ylabel('Attenuation Ratio');
xticks(0:1:8)
grid on;


% figure;
% semilogy(SD,SNR_4,'-or','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');
% xlabel('Source-detector distance (cm)');
% ylabel('Attenuation Ratio');
% xticks(0:1:8)
% grid on;

% fig=figure;
% set(fig,'defaultAxesColorOrder',[[0 0 1];[1 0 0]]);
% yyaxis right
% plot(SD,sensitive_4,'-r','LineWidth',2,'MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r');
% hold on
% yyaxis left
% plot(SD,attenuation_t,'-b','LineWidth',2,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b');
% hold on
% plot(SD,attenuation_M,':b','LineWidth',2,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b');
% hold on
% plot(SD,SNR_4,'--g','LineWidth',2,'MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','g');
% set(findall(0,'YAxisLocation','right'),'Yscale','linear');
% set(findall(0,'YAxisLocation','left'),'Yscale','log');
% xlabel('Source-detector distance (cm)');
% xticks(0:1:8)
% grid on
% grid minor