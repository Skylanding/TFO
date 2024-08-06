clc;
clear all;

set(0,'defaultAxesFontSize',20,'DefaultLineLineWidth', 1)
    sensitive_4 = zeros(1,20);
    attenuation_t = zeros(1,20);
    attenuation_M = zeros(1,20);
    SNR_4 = zeros(1,20);
    SD = zeros(1,20);
    number_deep = zeros(1,20);
    number_M = zeros(1,20);
    number_T = zeros(1,20);

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

    % which layer is the first layer, and what is the max num of layers
    first_layer = 1;
    maxlayer = 8;

    % fetal depth   %%(Unit: 1.8cm)
    depth = 18;
    derm_depth = depth - (1+7);   % cass1: fetal - amniotic - uterus
%     derm_depth = depth - (1+16);   % case2: fetal - amniotic - uterus
%     derm_depth = depth - (1+11.5);   % case3: fetal - amniotic - uterus

    dep = derm_depth;

    precision = 2;      %0.5cm (0.5mm??)
    diminmm = 100;      %10cm
    cfg.nphoton = 120e6;

    %number of pixels
    dim = diminmm * precision;

%     % model and simulation parameters configuration (layers thickness) 
%     % (1:maternal dermal, 2:maternal subdermal, 3:maternal uterus, 4:amniotic fluid, 5:fetal scalp, 6:fetal arterial, 7:fetal skull, 8:fetal brain)
%     airlay = 1 * precision;
%     lay = zeros(1,maxlayer);
%     lay(first_layer) = 2 * precision;
%     lay(2) = 13 * precision;
%     lay(3) = 16 * precision;
%     lay(4) = 1 * precision;
%     lay(5) = 2 * precision;
%     lay(6) = 1 * precision;
%     lay(7) = 2 * precision; 
%     lay(8) = 100 * precision;

    % defining the volume
    [xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);
    cfg.vol=uint8(zeros(dim,dim,dim));

    % more realistic with 4 layers (1:air/fat, 2:maternal dermal, 3:maternal uterus, 4:amniotic fluid)
%     % case 1 extreme thin
%     airlay_thick = 1; 
%     derm_thick = dep; %now fix (1.6cm depth, depth =17;)
%     uter_thick = 7; 
%     amio_thick = 1; 

    % case 2 extreme thick
    airlay_thick = 1; 
    derm_thick = dep; %now fix (1.6cm depth, depth =17;)
    uter_thick = 7; 
    amio_thick = 1; 

%     % case 3 extreme thick
%     airlay_thick = 1; 
%     derm_thick = dep; %now fix (1.6cm depth, depth =17;)
%     uter_thick = 11.5; 
%     amio_thick = 1; 

    interest_layer = 4;

    % fat tissue
    xfat = diminmm/2;
    yfat = diminmm/2;
    rfat = 160; %%mm
    zfat = rfat + airlay_thick; % in mm
    center_fat_sphere = [xfat*precision yfat*precision zfat*precision];
    dist = (xi-center_fat_sphere(1)).^2+(yi-center_fat_sphere(2)).^2+(zi-center_fat_sphere(3)).^2;
    cfg.vol(dist<(rfat*precision).^2)=1;

    % uterus
    xuterus = diminmm/2;
    yuterus = diminmm/2;
    ruterus = 130; %mm
    zuterus = ruterus + airlay_thick + derm_thick;
    center_uterus_sphere = [xuterus*precision yuterus*precision zuterus*precision];
    dist = (xi-center_uterus_sphere(1)).^2+(yi-center_uterus_sphere(2)).^2+(zi-center_uterus_sphere(3)).^2;
    cfg.vol(dist<(ruterus*precision).^2)=2;

    % amiotic fluid
    xaf = diminmm/2;
    yaf = diminmm/2;
    raf = 120; %mm
    zaf = raf + airlay_thick + derm_thick + uter_thick;
    center_af_sphere = [xaf*precision yaf*precision zaf*precision];
    dist = (xi-center_af_sphere(1)).^2+(yi-center_af_sphere(2)).^2+(zi-center_af_sphere(3)).^2;
    cfg.vol(dist<(raf*precision).^2)=3;

    % fetus 
    xfetus = diminmm/2;
    yfetus = diminmm/2;
    rfetus = 70; %mm
    zfetus = rfetus + airlay_thick + derm_thick + uter_thick + amio_thick;
    center_fetus_sphere = [xfetus*precision yfetus*precision zfetus*precision];
    dist = (xi-center_fetus_sphere(1)).^2+(yi-center_fetus_sphere(2)).^2+(zi-center_fetus_sphere(3)).^2;
    cfg.vol(dist<(rfetus*precision).^2)=4;

    % fetus body
    xfetusb = diminmm/2;
    yfetusb = diminmm/2 + diminmm/4;
    rfetusb = 120; %mm
    zfetusb = rfetusb + airlay_thick + derm_thick + uter_thick + amio_thick + rfetus + rfetus/2;
    center_fetusb_sphere = [xfetusb*precision yfetusb*precision zfetusb*precision];
    dist = (xi-center_fetusb_sphere(1)).^2+(yi-center_fetusb_sphere(2)).^2+(zi-center_fetusb_sphere(3)).^2;
    cfg.vol(dist<(rfetusb*precision).^2)=4;

    % (1:maternal dermal, 2:maternal subdermal, 3:maternal uterus, 4:amniotic fluid, 5:fetal scalp, 6:fetal arterial, 7:fetal skull, 8:fetal brain)
    lay(first_layer) = 1 * precision;
    lay(2) = (dep - 1) * precision;
    lay(3) = 7 * precision;
    lay(4) = 1 * precision;
    lay(5) = 2 * precision;
    lay(6) = 1 * precision;
    lay(7) = 2 * precision; 
    lay(8) = (rfetus-2-1-2) * precision;

    % rounding
    cfg.vol=uint8(cfg.vol);
    cfg.unitinmm=1/precision;  %this statement specifies the resultion in mm

    datadim =[dim dim dim];

    % with 3+1 layers
    air_prop = [0 0 1 1];
    
    % wavelength 850nm
    matderm_prop    =[0.0125  17.7  0.9  1.4];
    matsubderm_prop =[0.0088  11.1  0.9  1.4];
    uter_prop       =[0.01    8.15  0.9  1.4] ;
    amni_prop       =[0.0042  0.1   0.9  1.334];
    fetscalp_prop   =[0.0157  6.23  0.9  1.3];
    fetart_prop     =[0.0155  30    0.9  1.3];
    fetskull_prop   =[0.0215  9.1   0.9  1.3];
    fet_brain       =[0.0132  9.8   0.9  1.3];

    % wavelength 735nm
%     matderm_prop    =[0.017  23   0.9 1.4];
%     matsubderm_prop =[0.0085 12  0.9  1.4];
%     uter_prop       =[0.016  10.8  0.9  1.4] ;
%     amni_prop       =[0.0025 0.1   0.9  1.334];
%     fetscalp_prop   =[0.0157 6.81  0.9  1.3];
%     fetart_prop     =[0.0175 35   0.9  1.3];
%     fetskull_prop   =[0.021  10.9   0.9  1.3];
%     fet_brain       =[0.0187 12.2   0.9  1.3];

    tot_derm_prop = (matderm_prop*lay(1)+matsubderm_prop*lay(2))/(lay(1)+lay(2));
    tot_fet_prop = (fet_brain*lay(8) + fetskull_prop*lay(7)+ fetart_prop*lay(6) + fetscalp_prop*lay(5))/(lay(8)+lay(7)+lay(6)+lay(5));

%     tot_derm_prop(1,2)=tot_derm_prop(1,2)*20;
    cfg.prop =[air_prop;           % air
               tot_derm_prop;        % maternal derm + subderm
               uter_prop;            % uterus
               amni_prop;            % amniotic fluid
               tot_fet_prop];        % fetal sclap
     % [mua,mus,g,n] in 1/mm

    %%%%%%%%%%%%%%%%%%%%%%%

    %% time it runs
    cfg.tstart=0;
    cfg.tend=5e-9;
    cfg.tstep=5e-10;

%      %% other properties
%      cfg.issrcfrom0 = 1;
%      cfg.gpuid = 1;
%      cfg.autopilot = 1;
%      cfg.isnormalized = 1;
%      cfg.isreflect = 0;

    % GPU thread configuration
    cfg.autopilot=1;
    cfg.gpuid=1;

%     figure;
%     crosssectionloc = diminmm/2 ; %mm
%     mol =squeeze(cfg.vol(crosssectionloc*precision,:,:));
%     imagesc(0.1*(0.001:1/precision:dim/precision),0.1*(0.001:1/precision:dim/precision),mol);
%     xlabel('z [cm]');
%     ylabel('x [cm]');
%     set(gca,'XTickLabels',[0:2:10])
%     set(gca,'YTickLabels',[10:-2:0])
%     title(sprintf('domain cross section at x=%.2fmm',crosssectionloc)); 


% for i=1:1:32
i=6;    %SD=12mm
    % define the source position
    cfg.srctype='pencil';
    xsrc = dim/2 - 2*i;
    ysrc = dim/2;  % move the source to the left for 1mm  (SD=2mm, 4mm, 6mm ,....)
    tmp = find(cfg.vol(xsrc, ysrc,:));
    zsrc = tmp(1);
    cfg.srcpos=[xsrc, ysrc, zsrc]; %a 1 by 3 vector, the position of the source in grid unit
    cfg.srcdir=[0 0 1];           % a 1 by 3 vector, specifying the incident vector; if srcdir
                                  % contains a 4th element, it specifies the focal length of
                                  % the source (only valid for focuable src, such as planar, disk,
                                  % fourier, gaussian, pattern, slit, etc); if the focal length
                                  % is nan, all photons will be launched isotropically regardless
                                  % of the srcdir direction.

    %% running simulation with boundary reflection enabled
    fprintf('running simulation ... this takes about 50 seconds on a GTX 470\n');
    cfg.isreflect=1; % enable reflection at exterior boundary
    cfg.isrefint=1;  % enable reflection at interior boundary too
    cfg.issavedet=1; % enable recording partial pathlength of detected photons

    xdet = dim/2 + 2*i;
    ydet = dim/2;  % move the source to the left for 1mm  (SD=2mm, 4mm, 6mm ,....)
    tmp=find(cfg.vol(round(xdet),round(ydet),:));
    zdet = tmp(1);
    cfg.detpos=[xdet ydet zdet 0.1*10*precision];
    tic;
    [f2,det2]=mcxlab(cfg);

% % %  picture the CW
%{
% %     calculate the fluence distribution with the given config
% %     calculate run time
%     [fluence,detpt,vol,~,~] = mcxlab(cfg);
% %     toc;
%     cwfluence(1,:,:,:)=sum(fluence.data,4) ;  % fluence rate sum time gate
%     cwfluence_avg = squeeze(sum(cwfluence,1));
% 
%     crosssectionloc = diminmm/2 ; %mm
%     mol =squeeze(cfg.vol(crosssectionloc*precision,:,:));
%     mol1=mol;
%     BW1 = edge(mol1,'Roberts');
%     [row,col]=find(BW1==1);
% 
%     figure;
%     Img=imagesc(0.1*(0.001:1/precision:dim/precision),0.1*(0.001:1/precision:dim/precision),squeeze(log(cwfluence_avg(:,cfg.srcpos(2),:))));
% %     title('fluence at y=50mm');
%     xlabel('z [cm]');
%     ylabel('x [cm]');
%     fet_progress=colorbar;
% %     set(gca,'XTickLabels',[0:2:10])
% %     set(gca,'YTickLabels',[10:-2:0])
%     hold on
%     for k=1:length(col)
%     plot(0.1.*(col(k)-1)./precision,0.1.*row(k)./precision,'marker','.','MarkerSize',4,'Color','black')
%     hold on
%     end
%}
    %% plot the results
%{
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

%}
    SD(i)= (xdet-xsrc)*1/precision*0.1;

%% weight mu_a*pathlength
    
    mua_v=transpose(cfg.prop(2:5,1));
    mua_layer = repmat(mua_v, size(det2.ppath,1),1);
    layerwise_weight = det2.ppath*cfg.unitinmm.*mua_layer;
    BBL_temp = exp(-layerwise_weight);
    BBL_temp(BBL_temp==1) = 0;
    number_deep(i) = sum(exp(-sum(layerwise_weight(find(det2.ppath(:,4)),:)')))/cfg.nphoton;
    number_M(i) = sum(exp(-sum(layerwise_weight(:,1:3)')))/cfg.nphoton;
    number_T(i) = sum(exp(-sum(layerwise_weight')))/cfg.nphoton;

% end



fig=figure;
set(fig,'defaultAxesColorOrder',[[0 0 0];[0 0 1]]);
plot(SD(2:end),number_T(2:end),'-b','LineWidth',2,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b');
hold on
plot(SD(2:end),number_deep(2:end),'--b','LineWidth',2,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b');
set(findall(0,'YAxisLocation','left'),'Yscale','log');
xlabel('Source-detector distance (cm)');
xticks(0:1:8)
grid on
grid minor



