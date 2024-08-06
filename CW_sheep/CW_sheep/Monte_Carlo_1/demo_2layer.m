for i = 1 : 10
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


    precision = 1;      %0.5cm
    diminmm = 100;      %10cm
    cfg.nphoton=5e7;

    %number of pixels
    dim = diminmm * precision;

    % define a 4 layer structure
    cfg.vol=ones(dim,dim,100);
    cfg.vol(:,:,3*precision:end)=2;
    cfg.vol=uint8(cfg.vol);
    % cfg.unitinmm=1/precision;  %this statement specifies the resultion in mm

    % define the source position
    cfg.srcpos=[50*precision,50*precision,0]+1;
    cfg.srcdir=[0 0 1];
    cfg.unitinmm=1/precision;  %this statement specifies the resultion in mm

    % use the brain optical properties defined at
    % http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/CollinsAtlasMesh
    % format: [mua(1/mm) mus(1/mm) g n]

    cfg.prop=[0 0 1 1            % medium 0: the environment
       0.02 6.5  0.8  1.4     % medium 1: skin & skull
       0.01 6    0.8  1.4];   % medium 2: white matter

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
    cfg.detpos=[31*precision 51*precision 1*precision 2*precision];
    tic;
    [f2,det2]=mcxlab(cfg);
    toc;

    %% plot the results
    figure;
    contourf(log10(squeeze(sum(f1.data(:,51*precision,:,:),4))'),1:0.5:8);
    hold on
    plot([0 diminmm],[3 3],'--r');
    title('flux with no reflection');
    set(gca,'clim',[1 8]);
    colorbar

    figure;
    contourf(log10(squeeze(sum(f2.data(:,51*precision,:,:),4))'),1:0.5:8);
    hold on
    plot([0 diminmm],[3 3],'--r');
    title('flux with reflection at boundaries');
    set(gca,'clim',[1 8]);
    colorbar

    figure;
    xi=1e9*((cfg.tstart:cfg.tstep:cfg.tend-cfg.tstep)+0.5*cfg.tstep);
    semilogy(xi,squeeze(f2.data(51,41,1,:)));
    xlabel('time (ns)')
    ylabel('flux (1/mm^2/s)')
    set(gca,'yscale','log');
    title('time point spread functions (TPSF)');

    figure;
    [hs1,c1]=hist(det2.ppath(find(det2.ppath(:,1)),1),200);
    hold on;
    [hs2,c2]=hist(det2.ppath(find(det2.ppath(:,2)),2),200);
    bar(c1,hs1,'edgecolor','none','facecolor','r');
    bar(c2,hs2,'edgecolor','none','facecolor','b');
    legend('tissue 1','tissue 2');
    xlabel('partial pathlength (mm)');
    title(sprintf('detected %d photons',size(det2.ppath,1)));
end