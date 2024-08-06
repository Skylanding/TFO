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
clear all;

%% preparing the input data
% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E'); 

cfg.nphoton=5e6;

% define a 4 layer structure
cfg.vol=ones(100,100,50);
cfg.vol(:,:,20:end)=2;
cfg.vol=uint8(cfg.vol);

% define the source position
cfg.srcpos=[50,50,0]+1;
cfg.srcdir=[0 0 1];

% use the brain optical properties defined at
% http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/CollinsAtlasMesh
% format: [mua(1/mm) mus(1/mm) g n]
% paper
% cfg.prop=[0 0 1 1            % medium 0: the environment
%    0.02 6.5  0.8  1.4     % medium 1: skin & skull
%    0.01 0.7/0.2    0.8  1.4];   % medium 2: white matter

% cfg.srctype='pencil';

g=0.8;
cfg.prop=[0 0 1 1            % medium 0: the environment
   1.3 0.02  g  1.4     % medium 1: skin & skull
   0.7 0.01   g  1.4];   % medium 2: white matter

% time-domain simulation parameters
cfg.tstart=0;
cfg.tend=40e-9;
cfg.tstep=0.5e-10;

% GPU thread configuration
cfg.autopilot=1;
cfg.gpuid=1;

%% running simulation without boundary reflection
% fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
% cfg.isreflect=0; % disable reflection at exterior boundary
% tic;
% f1=mcxlab(cfg);
% toc;

%% running simulation with boundary reflection enabled
fprintf('running simulation ... this takes about 50 seconds on a GTX 470\n');
cfg.isreflect=1; % enable reflection at exterior boundary
cfg.isrefint=1;  % enable reflection at interior boundary too
cfg.issavedet=1; % enable recording partial pathlength of detected photons
cfg.detpos=[50 60 1 2];
tic;
[f2,det2]=mcxlab(cfg);
toc;

%% plot the results
% figure;
% contourf(log10(squeeze(sum(f1.data(:,51,:,:),4))'),1:0.5:8);
% hold on
% plot([0 100],[21 21],'--r');
% title('flux with no reflection');
% set(gca,'clim',[1 8]);
% 
figure
contourf(log10(squeeze(sum(f2.data(:,51,:,:),4))'),1:0.5:8);
hold on
plot([0 100],[26 26],'--r');
title('flux with reflection at boundaries');
set(gca,'clim',[1 8]);


% % figure
% xi=1e9*((cfg.tstart:cfg.tstep:cfg.tend-cfg.tstep)+0.5*cfg.tstep);
% plot(xi,squeeze(f2.data(50,90,1,:))./cfg.nphoton);
% xlabel('time (ns)')
% ylabel('flux (1/mm^2/s)')
% xlim([0 6])
% % set(gca,'yscale','log');
% title('time-resolved reflectance (TFO)');
% grid on;
% hold on;

% figure;
% hold on;
% [hs1,c1]=hist(det2.ppath(find(det2.ppath(:,1)),1),200);
% [hs2,c2]=hist(det2.ppath(find(det2.ppath(:,2)),2),200);
% bar(c1,hs1,'edgecolor','none','facecolor','r');
% bar(c2,hs2,'edgecolor','none','facecolor','b');
% legend('tissue 1','tissue 2');
% xlabel('partial pathlength (mm)');
% title(sprintf('detected %d photons',size(det2.ppath,1)));

% pi=1:1:20;
% for i=1:20
% photon(i)=sum(squeeze(f2.data(50,50+i,1,:))./cfg.nphoton)/2000;
% end
% semilogy(pi,photon);
% xlabel('distance (mm)')
% ylabel('flux (1/mm^2/s)')
% title('reflectance');
% grid on;
% hold on;
