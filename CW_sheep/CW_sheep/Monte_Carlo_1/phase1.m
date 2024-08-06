%%% project phase 1 (detectors location)

% Detector distances (SD): 1.5, 3, 4.5, 7, 10 cm
% Fetal depths: 1 cm

% 5.5 by 5.5 square 
% smbb850ds gaussian dist
% 10% pump count photons

% Wavelength: (850nm) (740nm if enough time?)

% Target metrics:
% - dependent vars: Fetal signal sensitivity, Partial path-length through tissues < L >
% - independent: SD distance, fetal depth
% 
% Target Plot: Sensitivity plots vs SD distance and Fetal Depths (save the raw values so I can Plot them manually)


clc;
clear all;

%justmodel = 1;    %only model
justmodel = 0;

% which layer is the first layer, and what is the max num of layers
first_layer = 1;
maxlayer = 8;

% number of the simulations to be averaged
num_sim = 1;
%num_sim = 10;

% fetal depth
depth = 10;

derm_depth = depth - (1+7);   %fetal - amniotic - uterus

% number of different fetal depth
maxd = 1;             

% number of Wavelengths (735nm = 1 850nm = 2)
maxWV = 1;          
%maxWV = 2; 
            
for WV = 1:maxWV

    for d=1:maxd

        dep = derm_depth(d);

        % info = mcxlab('gpuinfo')

        % size and precision and number of photons
        % define the simulation using TFO curved
        % precision = 1;
        % maxdiminmm = 200;
        % cfg.nphoton = 5e6;

        precision = 2;      %0.5cm
        diminmm = 200;      %20cm
        cfg.nphoton = 5e6;

        % precision = 3;
        % maxdiminmm = 100;
        % cfg.nphoton = 50e6*24;

        % precision = 4;
        % maxdiminmm = 60;
        % cfg.nphoton = 50e6*64;

        %number of pixels
        dim = diminmm * precision;

        % model and simulation parameters configuration (layers thickness) 
        % (1:maternal dermal, 2:maternal subdermal, 3:maternal uterus, 4:amniotic fluid, 5:fetal scalp, 6:fetal arterial, 7:fetal skull, 8:fetal brain)
        airlay = 1 * precision;
        lay = zeros(1,maxlayer);
        lay(first_layer) = 2 * precision;
        lay(2) = 13 * precision;
        lay(3) = 7 * precision;
        lay(4) = 1 * precision;
        lay(5) = 2 * precision;
        lay(6) = 1 * precision;
        lay(7) = 2 * precision; 
        lay(8) = 100 * precision;

        % defining the volume
        [xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);
        cfg.vol=uint8(zeros(dim,dim,dim));

        % cubic layer
        % interest_layer = 8;
        % stacklayer = 0;
        % stacklayer = stacklayer + airlay;
        % cfg.vol(:,:,stacklayer+1:stacklayer+lay(first_layer))=1;
        % stacklayer = stacklayer + lay(first_layer);
        % cfg.vol(:,:,stacklayer+1:stacklayer+lay(2))=2;
        % stacklayer = stacklayer + lay(2);
        % cfg.vol(:,:,stacklayer+1:stacklayer+lay(3))=3;
        % stacklayer = stacklayer + lay(3);
        % cfg.vol(:,:,stacklayer+1:stacklayer+lay(4))=4;
        % stacklayer = stacklayer + lay(4);
        % cfg.vol(:,:,stacklayer+1:stacklayer+lay(5))=5;
        % stacklayer = stacklayer + lay(5);
        % cfg.vol(:,:,stacklayer+1:stacklayer+lay(6))=6;
        % stacklayer = stacklayer + lay(6);
        % cfg.vol(:,:,stacklayer+1:stacklayer+lay(7))=7;
        % stacklayer = stacklayer + lay(7);
        % cfg.vol(:,:,stacklayer+1:dim)=8;

        % concenteric spheres
        % interest_layer = 8;
        % xfat = diminmm/2;
        % yfat = diminmm/2;
        % zfat = diminmm; % in mm
        % center_fat_sphere = [xfat*precision yfat*precision zfat*precision];
        % dist=(xi-center_fat_sphere(1)).^2+(yi-center_fat_sphere(2)).^2+(zi-center_fat_sphere(3)).^2;
        % stacklayer=0;
        % stacklayer = stacklayer + airlay;
        % cfg.vol(dist<(dim-(stacklayer)).^2)=1;
        % stacklayer = stacklayer + lay(first_layer);
        % cfg.vol(dist<(dim-(stacklayer)).^2)=2;
        % stacklayer = stacklayer + lay(2);
        % cfg.vol(dist<(dim-(stacklayer)).^2)=3;
        % stacklayer = stacklayer + lay(3);
        % cfg.vol(dist<(dim-(stacklayer)).^2)=4;
        % stacklayer = stacklayer + lay(4);
        % cfg.vol(dist<(dim-(stacklayer)).^2)=5;
        % stacklayer = stacklayer + lay(5);
        % cfg.vol(dist<(dim-(stacklayer)).^2)=6;
        % stacklayer = stacklayer + lay(6);
        % cfg.vol(dist<(dim-(stacklayer)).^2)=7;
        % stacklayer = stacklayer + lay(7);
        % cfg.vol(dist<(dim-(stacklayer)).^2)=8;


        % more realistic with 4 layers (1:air/fat, 2:maternal dermal, 3:maternal uterus, 4:amniotic fluid)
        airlay_thick = 1; 
        derm_thick = dep; %now fix
        uter_thick = 7; 
        amio_thick = 1; 

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

        % medium properties 
        % FOR LMBDA = 735NM
        % in Dan's paper the data is in cm^-1
        % with 8+1 layers
        % cfg.prop=[0 0 1 1;        % air
        % 0.017  23   0.9 1.4;      % maternal dermal
        % 0.0085 12   0.9 1.4;      % maternal subdermal
        % 0.016  10.8 0.9 1.4;      % uterus
        % 0.0025 0.1  0.9 1.334;    % amniotic fluid
        % 0.0157 6.81 0.9 1.3;      % fetal scalp
        % 0.0175 35   0.9 1.3;      % fetal arterial         
        % 0.021  10.9 0.9 1.3;      % fetal skull
        % 0.0187 12.2 0.9 1.3];     % fetal brain
        % [mua,mus,g,n] in 1/mm

        % FOR LMBDA = 850NM
        % with 8+1 layers
        % cfg.prop=[0 0 1 1;        % air
        % 0.0125 17.7 0.9 1.4;      % maternal dermal
        % 0.0088 11.1 0.9 1.4;      % maternal subdermal
        % 0.01   8.15 0.9 1.4;      % uterus
        % 0.0042 0.1  0.9 1.334;    % amniotic fluid
        % 0.0157 6.23 0.9 1.3;      % fetal scalp
        % 0.0155 30   0.9 1.3;      % fetal arterial         
        % 0.0215 9.1  0.9 1.3;      % fetal skull
        % 0.0132 9.8  0.9 1.3];     % fetal brain
        % [mua,mus,g,n] in 1/mm

        % with 3+1 layers
        air_prop = [0 0 1 1];
        if WV ==1 
            matderm_prop    =[0.017   23    0.9  1.4];
            matsubderm_prop =[0.0085  12    0.9  1.4];
            uter_prop       =[0.016   10.8  0.9  1.4] ;
            amni_prop       =[0.0025  0.1   0.9  1.334];
            fetscalp_prop   =[0.0157  6.81  0.9  1.3];
            fetart_prop     =[0.0175  35    0.9  1.3];
            fetskull_prop   =[0.021   10.9  0.9  1.3];
            fet_brain       =[0.0187  12.2  0.9  1.3];
        elseif WV ==2
            matderm_prop    =[0.0125  17.7  0.9  1.4];
            matsubderm_prop =[0.0088  11.1  0.9  1.4];
            uter_prop       =[0.01    8.15  0.9  1.4] ;
            amni_prop       =[0.0042  0.1   0.9  1.334];
            fetscalp_prop   =[0.0157  6.23  0.9  1.3];
            fetart_prop     =[0.0155  30    0.9  1.3];
            fetskull_prop   =[0.0215  9.1   0.9  1.3];
            fet_brain       =[0.0132  9.8   0.9  1.3]; 
        end

        tot_derm_prop = (matderm_prop*lay(1)+matsubderm_prop*lay(2))/(lay(1)+lay(2));
        tot_fet_prop = (fet_brain*lay(8) + fetskull_prop*lay(7)+ fetart_prop*lay(6) + fetscalp_prop*lay(5))/(lay(8)+lay(7)+lay(6)+lay(5));

        cfg.prop =[air_prop;           % air
                 tot_derm_prop;        % maternal derm + subderm
                 uter_prop;            % uterus
                 amni_prop;            % amniotic fluid
                 tot_fet_prop];        % fetal sclap
        % [mua,mus,g,n] in 1/mm


        % source position

        % src type gaussian
        % waistrad = 5.5 ;      %waist radius of beam mm
        % cfg.srctype='gaussian';
        % cfg.srcparam1=[round(waistrad*precision) 0 0 0];
        % cfg.srcparam2=[0 0 0 0];

        % src type pencil beam
        cfg.srctype='pencil';
        xsrc = (3*diminmm) /4 ; %in mm (changing detectors and source together)
        ysrc = (2*diminmm) /4 + 5;
        tmp = find(cfg.vol(xsrc*precision, ysrc*precision,:));
        srcdistsurface = tmp(1);
        
        % if sphere
        srcdirvec = [center_fat_sphere(1)-xsrc*precision center_fat_sphere(2)-ysrc*precision center_fat_sphere(3)-srcdistsurface];
        cfg.srcdir = srcdirvec/norm(srcdirvec);

        % if box
        % cfg.srcdir=[0 0 1];

        % src directly on the surface
        cfg.srcpos = [xsrc*precision ysrc*precision srcdistsurface];

        % src with a mm distance form the surface
        % distfromsurf = 1;
        % cfg.srcpos=[xsrc*precision ysrc*precision srcdistsurface]-cfg.srcdir*distfromsurf*precision;

        cfg.srcpos = round(cfg.srcpos);

        %% other properties
        cfg.issrcfrom0 = 1;
        cfg.gpuid = 1;
        cfg.autopilot = 1;
        cfg.isnormalized = 1;
        cfg.isreflect = 0;

        %save trajectory
        cfg.maxjumpdebug = 1;
        % cfg.maxjumpdebug = cfg.nphoton;
        cfg.maxdetphoton = cfg.nphoton;
        cfg.issaveref = 1;

        %% time it runs
        cfg.tstart=0;
        cfg.tend=5e-9;
        cfg.tstep=5e-10;
        % cfg.tend=1e-6;
        % cfg.tstep=1e-6/10;

        num_gates = (cfg.tend-cfg.tstart)/cfg.tstep;
        timegates = [cfg.tstart:cfg.tstep:cfg.tend];

        % detector(s) position
        % detector separation
        % mindistdet = 5;
        % spacedet = 5;
        % maxdisdet = diminmm/2;
        % disdet = (mindistdet:spacedet:maxdisdet);
        disdet =[15, 30, 45, 70, 100];  %mm
        
        % change radius in longer disances
        detrad = 2; %mm
        phaseshift = (-2*pi)/4;
        % phaseshift = 0;
        theta = 0+phaseshift:pi:pi+phaseshift-0.0001;
        cfg.detpos=[];

        for j=1:length(disdet)
            xdet = cfg.srcpos(1)/precision + disdet(j)*cos(theta);
            ydet = cfg.srcpos(2)/precision + disdet(j)*sin(theta);

                if (round(xdet*precision)<= diminmm*precision) && (round(xdet*precision)>=1) && (round(ydet*precision)<=diminmm*precision) && (round(ydet*precision)>=1)
                    tmp=find(cfg.vol(round(xdet*precision),round(ydet*precision),:));
                    detdistsurface=tmp(1);
                    cfg.detpos = [cfg.detpos; round(xdet*precision) round(ydet*precision) detdistsurface detrad*precision];
                end
        end

        % creating random seed
        s='0123456789ABCDEF';
        numRands = length(s);
        sLength = 8;
        randString = s(ceil(rand(1,sLength)*numRands));

        % simulating for num_sim times
        ratio = zeros(num_sim,length(disdet));
        ntotphoton = zeros(num_sim,length(disdet));
        nintlayphoton = zeros(num_sim,length(disdet));
        
        meanpartpathintlayphoton = zeros(num_sim,length(disdet));
        stdpartpathdintlayphoton = zeros(num_sim,length(disdet));
        maxpartpathdintlayphoton = zeros(num_sim,length(disdet));
        minpartpathdintlayphoton = zeros(num_sim,length(disdet));
        
        %reached layer4
        meanpartpathintlayphoton2 = zeros(num_sim,length(disdet),4);
        stdpartpathdintlayphoton2 = zeros(num_sim,length(disdet),4);
        maxpartpathdintlayphoton2 = zeros(num_sim,length(disdet),4);
        minpartpathdintlayphoton2 = zeros(num_sim,length(disdet),4);
        
        dists = zeros(num_sim,length(disdet));
        cwfluence = zeros(num_sim,dim,dim,dim);
        cwdref = zeros(num_sim,dim,dim,dim);

        if ~justmodel
            progress_bar = waitbar(0,'Please wait...');
            for n=1:num_sim
                waitbar((n/num_sim)*(d/maxd)*(WV/maxWV),progress_bar,sprintf('Simulation number = %d, Light WV=%d, depth d=%d',n,WV,d)); 
                clear fluence
                clear detpt

                % generating random seed everytime
                % it seems it was overflowing with 32 bits, so divide by 2
                randString = s(ceil(rand(1,sLength)*numRands));
                cfg.seed=hex2dec(randString)/2; 
                % cfg.seed=hex2dec('623F9A9E'); 

                % calculate the fluence distribution with the given config
                % calculate run time
                [fluence,detpt,vol,~,~] = mcxlab(cfg);

                % calculating the ratio
                %%%% in detp.data columns are detected photons
                %%%% first row is id of detector which detected the photon
                %%%% partial path in layers start from row 2 -> interest_layer+1
                %%%% so partial path in fetal layer is in row interest_layer+1

                for j=1:length(disdet)
                    if ismember(j,detpt.data(1,:))  %if photon detected
                        indices_photon_detj = find(detpt.data(1,:)==j);
                        ntotphoton(n,j) = length(find(detpt.data(first_layer+1,indices_photon_detj)));
                        find_interst_lay_phot = find(detpt.data(interest_layer+1,indices_photon_detj));
                        nintlayphoton(n,j) = length(find_interst_lay_phot);
                        
                        if find_interst_lay_phot
                            sumpartialpaths = sum(detpt.data(2:interest_layer+1,find_interst_lay_phot),1) ;
                            meanpartpathintlayphoton(n,j) = mean(sumpartialpaths)/precision;
                            stdpartpathdintlayphoton(n,j) = std(sumpartialpaths)/precision;
                            maxpartpathdintlayphoton(n,j) = max(sumpartialpaths)/precision;
                            minpartpathdintlayphoton(n,j) = min(sumpartialpaths)/precision;
                        end
                         
                        if ntotphoton(n,j)
                            partialpaths2 = detpt.data(2:interest_layer+1,indices_photon_detj);
                            meanpartpathintlayphoton2(n,j,:) = squeeze(sum(partialpaths2, 2)/ntotphoton(n,j)/precision);
                            stdpartpathdintlayphoton2(n,j,:) = squeeze(std(partialpaths2, 0 , 2)/precision);
                            maxpartpathdintlayphoton2(n,j,:) = squeeze(max(partialpaths2, [], 2)/precision);
                            minpartpathdintlayphoton2(n,j,:) = squeeze(min(partialpaths2, [], 2)/precision);
                        end

                        ratio(n,j) =  nintlayphoton(n,j)/ntotphoton(n,j);
                        dists(n,j) =  disdet(j);

                    else
                        fprintf('NO PHOTON DETECTED in detector %d and depth %d and simulation number %d and WaveLength %d!!!!!\n',j,d,n,WV);
                    end
                end

                %% integrate time-axis (4th dimension) to get CW solutions
                cwfluence(n,:,:,:)=sum(fluence.data,4) ;  % fluence rate sum time gate
                cwdref(n,:,:,:)=sum(fluence.dref,4) ;     % diffuse reflectance sum time gate

            end
            close(progress_bar);    
        end

        %% sum and average all results
        ratio_avg = sum(ratio,1)/num_sim;
        ntotphoton_avg = sum(ntotphoton,1)/num_sim;
        nintlayphoton_avg = sum(nintlayphoton,1)/num_sim;
        
        meanpartpathintlayphoton_avg = sum(meanpartpathintlayphoton,1)/num_sim;
        stdpartpathdintlayphoton_avg = sum(stdpartpathdintlayphoton,1)/num_sim;
        maxpartpathdintlayphoton_avg = sum(maxpartpathdintlayphoton,1)/num_sim;
        minpartpathdintlayphoton_avg = sum(minpartpathdintlayphoton,1)/num_sim;

        meanpartpathintlayphoton2_avg = squeeze(sum(meanpartpathintlayphoton2,1)/num_sim);
        stdpartpathdintlayphoton2_avg = squeeze(sum(stdpartpathdintlayphoton2,1)/num_sim);
        maxpartpathdintlayphoton2_avg = squeeze(sum(maxpartpathdintlayphoton2,1)/num_sim);
        minpartpathdintlayphoton2_avg = squeeze(sum(minpartpathdintlayphoton2,1)/num_sim);

        dists_avg = sum(dists,1)/num_sim;

        cwfluence_avg = squeeze(sum(cwfluence,1)/num_sim);
        cwdref_avg = squeeze(sum(cwdref,1)/num_sim);

        %if justmodel
            % plot the volume
            figure(6);
            subplot(231);
            mcxpreview(cfg);
            title(sprintf('domain preview %2.2e photons precision = %.2fmm',cfg.nphoton,1/precision));
            figure(7);
            crosssectionloc = diminmm/2 ; %mm
            mol =squeeze(cfg.vol(crosssectionloc*precision,:,:));
            imagesc(0.001:1/precision:dim/precision,0.001:1/precision:dim/precision,mol);  
            title(sprintf('domain cross section at x=%.2fmm',crosssectionloc)); 
        %end 
 
        %plot the results 
        figure(5)
        subplot(232);
        imagesc(0.001:1/precision:dim/precision,0.001:1/precision:dim/precision,squeeze(log(cwfluence_avg(:,cfg.srcpos(2),:))));
        title('fluence at y=50mm');
        xlabel('z');
        ylabel('x');
        fet_progress=colorbar;
        lim=caxis;
        caxis([-11 22])
        ylabel(fet_progress,'log(W/mm^2)')
 
        subplot(233);
        imagesc(0.001:1/precision:dim/precision,0.001:1/precision:dim/precision,squeeze(log(cwdref_avg(:,:,1*precision))));
        title('diffuse refle. at z=1');
        xlabel('y');
        ylabel('x');       
        fet_progress=colorbar;
        lim=caxis;
        caxis([-11 22])
        ylabel(fet_progress,'log(W/mm^2)')

 
 
        subplot(234);
        imagesc(0.001:1/precision:dim/precision,0.001:1/precision:dim/precision,squeeze(log(cwdref_avg(:,:,1*precision+1))));
        title('diffuse refle. at z=1.5');
        xlabel('y');
        ylabel('x');  
        fet_progress=colorbar;
        lim=caxis;
        caxis([-11 22])
        ylabel(fet_progress,'log(W/mm^2)')

 
        subplot(235)
        imagesc(0.001:1/precision:dim/precision,0.001:1/precision:dim/precision,squeeze(log((cwdref_avg(:,:,1*precision)+cwdref_avg(:,:,1*precision+1)))));
        title('diffuse refle. at z=1 + diffuse refle. at z=1.5');
        xlabel('y');
        ylabel('x');  
        fet_progress=colorbar;
        lim=caxis;
        caxis([-11 22])
        ylabel(fet_progress,'log(W/mm^2)')

 
        subplot(236)
        imagesc(0.001:1/precision:dim/precision,0.001:1/precision:dim/precision,squeeze(log(sum(cwdref_avg,3))));
        title('diffuse refle. mapped to surface');
        xlabel('y');
        ylabel('x');  
        fet_progress=colorbar;
        lim=caxis;
        caxis([-11 22])
        ylabel(fet_progress,'log(W/mm^2)')

 
        % plot ratio of photons reach fetal layer to total photons detected
        figure(7)
        plot(dists_avg,ratio_avg);
        xlabel('distance from source(mm) at 45 degree offset')
        title('Ratio of photons reached the fetal layer to the total photons detected at the surface');

        %% RESULT STRUCT
        if WV ==1
            if d==1

                light735.depth10.ratio_avg = ratio_avg;
                light735.depth10.ntotphoton_avg = ntotphoton_avg;
                light735.depth10.nintlayphoton_avg = nintlayphoton_avg;
                light735.depth10.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light735.depth10.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light735.depth10.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light735.depth10.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light735.depth10.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light735.depth10.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light735.depth10.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light735.depth10.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light735.depth10.dists_avg = dists_avg;

                light735.depth10.ratio = ratio;
                light735.depth10.ntotphoton = ntotphoton;
                light735.depth10.nintlayphoton = nintlayphoton;
                light735.depth10.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light735.depth10.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light735.depth10.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light735.depth10.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light735.depth10.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light735.depth10.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light735.depth10.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light735.depth10.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light735.depth10.dists = dists;
                light735.depth10.photons = cfg.nphoton;
                light735.depth10.layerprop = cfg.prop;
                light735.depth10.volunitinmm = 1/precision;
                light735.depth10.diminmm = diminmm;
                light735.depth10.volume = cfg.vol;
                light735.depth10.sourcepos = cfg.srcpos/2;
                light735.depth10.time = timegates;

            elseif d==2

                light735.depth20.ratio_avg = ratio_avg;
                light735.depth20.ntotphoton_avg = ntotphoton_avg;
                light735.depth20.nintlayphoton_avg = nintlayphoton_avg;
                light735.depth20.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light735.depth20.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light735.depth20.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light735.depth20.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light735.depth20.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light735.depth20.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light735.depth20.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light735.depth20.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light735.depth20.dists_avg = dists_avg;

                light735.depth20.ratio = ratio;
                light735.depth20.ntotphoton = ntotphoton;
                light735.depth20.nintlayphoton = nintlayphoton;
                light735.depth20.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light735.depth20.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light735.depth20.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light735.depth20.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light735.depth20.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light735.depth20.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light735.depth20.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light735.depth20.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light735.depth20.dists = dists;
                light735.depth20.photons = cfg.nphoton;
                light735.depth20.layerprop = cfg.prop;
                light735.depth20.volunitinmm = 1/precision;
                light735.depth20.diminmm = diminmm;
                light735.depth20.volume = cfg.vol;
                light735.depth20.sourcepos = cfg.srcpos/2;
                light735.depth20.time = timegates;
                
            elseif d==3

                light735.depth30.ratio_avg = ratio_avg;
                light735.depth30.ntotphoton_avg = ntotphoton_avg;
                light735.depth30.nintlayphoton_avg = nintlayphoton_avg;
                light735.depth30.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light735.depth30.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light735.depth30.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light735.depth30.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light735.depth30.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light735.depth30.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light735.depth30.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light735.depth30.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light735.depth30.dists_avg = dists_avg;

                light735.depth30.ratio = ratio;
                light735.depth30.ntotphoton = ntotphoton;
                light735.depth30.nintlayphoton = nintlayphoton;
                light735.depth30.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light735.depth30.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light735.depth30.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light735.depth30.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light735.depth30.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light735.depth30.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light735.depth30.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light735.depth30.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light735.depth30.dists = dists;
                light735.depth30.photons = cfg.nphoton;
                light735.depth30.layerprop = cfg.prop;
                light735.depth30.volunitinmm = 1/precision;
                light735.depth30.diminmm = diminmm;
                light735.depth30.volume = cfg.vol;
                light735.depth30.sourcepos = cfg.srcpos/2;
                light735.depth30.time = timegates;
                
            elseif d==4

                light735.depth40.ratio_avg = ratio_avg;
                light735.depth40.ntotphoton_avg = ntotphoton_avg;
                light735.depth40.nintlayphoton_avg = nintlayphoton_avg;
                light735.depth40.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light735.depth40.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light735.depth40.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light735.depth40.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light735.depth40.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light735.depth40.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light735.depth40.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light735.depth40.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light735.depth40.dists_avg = dists_avg;

                light735.depth40.ratio = ratio;
                light735.depth40.ntotphoton = ntotphoton;
                light735.depth40.nintlayphoton = nintlayphoton;
                light735.depth40.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light735.depth40.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light735.depth40.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light735.depth40.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light735.depth40.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light735.depth40.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light735.depth40.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light735.depth40.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light735.depth40.dists = dists;
                light735.depth40.photons = cfg.nphoton;
                light735.depth40.layerprop = cfg.prop;
                light735.depth40.volunitinmm = 1/precision;
                light735.depth40.diminmm = diminmm;
                light735.depth40.volume = cfg.vol;
                light735.depth40.sourcepos = cfg.srcpos/2;
                light735.depth40.time = timegates;
                
            elseif d==5

                light735.depth50.ratio_avg = ratio_avg;
                light735.depth50.ntotphoton_avg = ntotphoton_avg;
                light735.depth50.nintlayphoton_avg = nintlayphoton_avg;
                light735.depth50.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light735.depth50.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light735.depth50.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light735.depth50.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light735.depth50.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light735.depth50.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light735.depth50.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light735.depth50.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light735.depth50.dists_avg = dists_avg;

                light735.depth50.ratio = ratio;
                light735.depth50.ntotphoton = ntotphoton;
                light735.depth50.nintlayphoton = nintlayphoton;
                light735.depth50.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light735.depth50.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light735.depth50.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light735.depth50.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light735.depth50.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light735.depth50.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light735.depth50.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light735.depth50.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light735.depth50.dists = dists;
                light735.depth50.photons = cfg.nphoton;
                light735.depth50.layerprop = cfg.prop;
                light735.depth50.volunitinmm = 1/precision;
                light735.depth50.diminmm = diminmm;
                light735.depth50.volume = cfg.vol;
                light735.depth50.sourcepos = cfg.srcpos/2;
                light735.depth50.time = timegates;
            end
        elseif WV ==2
            if d==1

                light850.depth10.ratio_avg = ratio_avg;
                light850.depth10.ntotphoton_avg = ntotphoton_avg;
                light850.depth10.nintlayphoton_avg = nintlayphoton_avg;
                light850.depth10.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light850.depth10.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light850.depth10.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light850.depth10.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light850.depth10.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light850.depth10.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light850.depth10.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light850.depth10.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light850.depth10.dists_avg = dists_avg;

                light850.depth10.ratio = ratio;
                light850.depth10.ntotphoton = ntotphoton;
                light850.depth10.nintlayphoton = nintlayphoton;
                light850.depth10.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light850.depth10.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light850.depth10.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light850.depth10.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light850.depth10.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light850.depth10.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light850.depth10.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light850.depth10.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light850.depth10.dists = dists;
                light850.depth10.photons = cfg.nphoton;
                light850.depth10.layerprop = cfg.prop;
                light850.depth10.volunitinmm = 1/precision;
                light850.depth10.diminmm = diminmm;
                light850.depth10.volume = cfg.vol;
                light850.depth10.sourcepos = cfg.srcpos/2;
                light850.depth10.time = timegates;
                
            elseif d==2

                light850.depth20.ratio_avg = ratio_avg;
                light850.depth20.ntotphoton_avg = ntotphoton_avg;
                light850.depth20.nintlayphoton_avg = nintlayphoton_avg;
                light850.depth20.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light850.depth20.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light850.depth20.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light850.depth20.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light850.depth20.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light850.depth20.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light850.depth20.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light850.depth20.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light850.depth20.dists_avg = dists_avg;

                light850.depth20.ratio = ratio;
                light850.depth20.ntotphoton = ntotphoton;
                light850.depth20.nintlayphoton = nintlayphoton;
                light850.depth20.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light850.depth20.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light850.depth20.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light850.depth20.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light850.depth20.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light850.depth20.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light850.depth20.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light850.depth20.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light850.depth20.dists = dists;
                light850.depth20.photons = cfg.nphoton;
                light850.depth20.layerprop = cfg.prop;
                light850.depth20.volunitinmm = 1/precision;
                light850.depth20.diminmm = diminmm;
                light850.depth20.volume = cfg.vol;
                light850.depth20.sourcepos = cfg.srcpos/2;
                light850.depth20.time = timegates;
                
            elseif d==3

                light850.depth30.ratio_avg = ratio_avg;
                light850.depth30.ntotphoton_avg = ntotphoton_avg;
                light850.depth30.nintlayphoton_avg = nintlayphoton_avg;
                light850.depth30.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light850.depth30.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light850.depth30.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light850.depth30.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light850.depth30.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light850.depth30.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light850.depth30.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light850.depth30.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light850.depth30.dists_avg = dists_avg;

                light850.depth30.ratio = ratio;
                light850.depth30.ntotphoton = ntotphoton;
                light850.depth30.nintlayphoton = nintlayphoton;
                light850.depth30.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light850.depth30.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light850.depth30.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light850.depth30.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light850.depth30.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light850.depth30.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light850.depth30.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light850.depth30.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light850.depth30.dists = dists;
                light850.depth30.photons = cfg.nphoton;
                light850.depth30.layerprop = cfg.prop;
                light850.depth30.volunitinmm = 1/precision;
                light850.depth30.diminmm = diminmm;
                light850.depth30.volume = cfg.vol;
                light850.depth30.sourcepos = cfg.srcpos/2;
                light850.depth30.time = timegates;
                
            elseif d==4

                light850.depth40.ratio_avg = ratio_avg;
                light850.depth40.ntotphoton_avg = ntotphoton_avg;
                light850.depth40.nintlayphoton_avg = nintlayphoton_avg;
                light850.depth40.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light850.depth40.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light850.depth40.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light850.depth40.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light850.depth40.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light850.depth40.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light850.depth40.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light850.depth40.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light850.depth40.dists_avg = dists_avg;

                light850.depth40.ratio = ratio;
                light850.depth40.ntotphoton = ntotphoton;
                light850.depth40.nintlayphoton = nintlayphoton;
                light850.depth40.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light850.depth40.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light850.depth40.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light850.depth40.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light850.depth40.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light850.depth40.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light850.depth40.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light850.depth40.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light850.depth40.dists = dists;
                light850.depth40.photons = cfg.nphoton;
                light850.depth40.layerprop = cfg.prop;
                light850.depth40.volunitinmm = 1/precision;
                light850.depth40.diminmm = diminmm;
                light850.depth40.volume = cfg.vol;
                light850.depth40.sourcepos = cfg.srcpos/2;
                light850.depth40.time = timegates;
                
            elseif d==5

                light850.depth50.ratio_avg = ratio_avg;
                light850.depth50.ntotphoton_avg = ntotphoton_avg;
                light850.depth50.nintlayphoton_avg = nintlayphoton_avg;
                light850.depth50.meanpartpathintlayphoton_avg = meanpartpathintlayphoton_avg;
                light850.depth50.stdpartpathdintlayphoton_avg = stdpartpathdintlayphoton_avg;
                light850.depth50.maxpartpathdintlayphoton_avg = maxpartpathdintlayphoton_avg;
                light850.depth50.minpartpathdintlayphoton_avg = minpartpathdintlayphoton_avg;
                light850.depth50.meanpartpathintlayphoton2_avg = meanpartpathintlayphoton2_avg;
                light850.depth50.stdpartpathdintlayphoton2_avg = stdpartpathdintlayphoton2_avg;
                light850.depth50.maxpartpathdintlayphoton2_avg = maxpartpathdintlayphoton2_avg;
                light850.depth50.minpartpathdintlayphoton2_avg = minpartpathdintlayphoton2_avg;
                
                light850.depth50.dists_avg = dists_avg;

                light850.depth50.ratio = ratio;
                light850.depth50.ntotphoton = ntotphoton;
                light850.depth50.nintlayphoton = nintlayphoton;
                light850.depth50.meanpartpathintlayphoton = meanpartpathintlayphoton;
                light850.depth50.stdpartpathdintlayphoton = stdpartpathdintlayphoton;
                light850.depth50.maxpartpathdintlayphoton = maxpartpathdintlayphoton;
                light850.depth50.minpartpathdintlayphoton = minpartpathdintlayphoton;
                light850.depth50.meanpartpathintlayphoton2 = meanpartpathintlayphoton2;
                light850.depth50.stdpartpathdintlayphoton2 = stdpartpathdintlayphoton2;
                light850.depth50.maxpartpathdintlayphoton2 = maxpartpathdintlayphoton2;
                light850.depth50.minpartpathdintlayphoton2 = minpartpathdintlayphoton2;
                
                light850.depth50.dists = dists;
                light850.depth50.photons = cfg.nphoton;
                light850.depth50.layerprop = cfg.prop;
                light850.depth50.volunitinmm = 1/precision;
                light850.depth50.diminmm = diminmm;
                light850.depth50.volume = cfg.vol;
                light850.depth50.sourcepos = cfg.srcpos/2;
                light850.depth50.time = timegates;
                
            end
        end

    end
end

if ~justmodel
    done = waitbar(1,'DONE !!');
end