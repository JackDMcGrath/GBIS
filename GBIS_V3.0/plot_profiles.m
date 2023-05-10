
gbis_folder='/nfs/a285/homes/eejdm/InSAR/combined_inversion/UTM_data/GBIS_output/AF_models';
gbis_inp='Little/little_noshear';
inversion_type='invert_GPS111_C';

gbismat=sprintf('%s/%s/%s/%s.mat',gbis_folder,gbis_inp,inversion_type,inversion_type);


figure
sgtitle(sprintf('Profile Comparison for GBIS %s/%s',gbis_inp,inversion_type))

prof_width=50;
prof_nums=[1 2 3];
binsize=5;
bins=-20:binsize:120;
bin_center=(bins(1)+0.5*binsize):binsize:(bins(end)-0.5*binsize);


nprof=size(prof_nums,2);
corners=sprintf('/nfs/a285/homes/eejdm/InSAR/combined_inversion/UTM_data/mcmc_profiles/%.0fkm_prof_corners.txt',prof_width);
fault=readmatrix(sprintf('/nfs/a285/homes/eejdm/InSAR/combined_inversion/UTM_data/mcmc_profiles/%.0fkm_outputs/faulttrace_shifted2.txt',prof_width));

corners=load(corners); % KM -> UTM Zoned

load(gbismat)

g=readmatrix('3comp_insarTC_fullres.txt');
fprintf('Full Res\n')
gps.ll=g(:,1:2);



[gps.utm(:,[1 2])]=lonlat2utm(gps.ll(:,[1 2])); % METERS -> UTM Zoned
gps.utm=gps.utm/1000; % KM -> UTM Zoned

par_ix=[];
par_ll=[];
for ii=1:size(invResults.model.parName,2)
    if strcmpi(invResults.model.parName{ii}(end),'X') || strcmpi(invResults.model.parName{ii}(end),'Y')
        par_ix=[par_ix,ii];
        par_ll=[par_ll,invResults.model.optimal(ii)];
    end
end

par_ll=reshape(par_ll,numel(par_ll)/2,2)/1000; %now in KM -> Local

par_ll=local2llh(par_ll',geo.referencePoint); % Convert GBIS co-ords to Lonlat CORRECT
[par_ll]=lonlat2utm(par_ll'); % Convert to UTM -> METERS -> UTM Zoned

for ii=1:size(par_ix,2)
    invResults.model.optimal(par_ix(ii))=par_ll(ii);
end

xlimit=[min(gps.utm(:,1)) max(gps.utm(:,1))];
ylimit=[min(gps.utm(:,2)) max(gps.utm(:,2))];

xbox=xlimit([1 1 2 2 1]);
ybox=ylimit([1 2 2 1 1]);

[xi,~]=polyxpoly(fault(:,1),fault(:,2),xbox,ybox);

fault1=find(fault(:,1)>min(xi));
fault1=fault1(1);

fault2=find(fault(:,1)<max(xi));
fault2=fault2(end);

fault=fault([fault1 fault2],:);

clear g fault1 fault2 xi yi xbox ybox xlimit ylimit

[dist,~] = point_to_line(gps.utm, fault(1,:), fault(2,:));

modGPS_original = forwardGPSModel(gps.utm*1000,invpar,invResults,modelInput)*1000; % Calculate modelled displacements (output in mm) Input GPS in meters

mod.ve=modGPS_original(1,:);
mod.vn=modGPS_original(2,:);

[mod]=GPS2FP(mod,90-55);

modGPS=[mod.vx;mod.vy;modGPS_original(3,:)];
fprintf('Converted ENU to FP\n')

in=cell(1,nprof);

for ii=1:nprof
    
    prof_id=prof_nums(ii);
    
    in{ii}=inpolygon(gps.utm(:,1),gps.utm(:,2),corners(:,prof_id*2-1),corners(:,prof_id*2));
    
end


for jj=1:size(bin_center,2)
    dist_ix=(dist>=bins(jj))+(dist<bins(jj+1));
    dist_ix=double(dist_ix==2); % Find all points in the right distance bin
    dist_ix(dist_ix==0)=10;
    for ii=1:nprof
        bin_ix=(dist_ix==in{ii});
        bindata{ii}(1,jj)=mean(modGPS(1,bin_ix));
        bindata{ii}(2,jj)=std(modGPS(1,bin_ix));
        bindata{ii}(3,jj)=mean(modGPS(2,bin_ix));
        bindata{ii}(4,jj)=std(modGPS(2,bin_ix));
        bindata{ii}(5,jj)=mean(modGPS(3,bin_ix));
        bindata{ii}(6,jj)=std(modGPS(3,bin_ix));
    end
end


colors=['r','g','b','k'];

if size(prof_nums,2)>1
    col=size(prof_nums,2)+1;
else
    col=size(prof_nums,2);
end

paray=[-40 0; floor(min(modGPS(1,:))) ceil(max(modGPS(1,:)))];
perpy=[-15 0; floor(min(modGPS(2,:))) ceil(max(modGPS(2,:)))];
verty=[-5 10; floor(min(modGPS(3,:))) ceil(max(modGPS(3,:)))];
paray=[min(paray(:,1)) max(paray(:,2))];
perpy=[min(perpy(:,1)) max(perpy(:,2))];
verty=[min(verty(:,1)) max(verty(:,2))];


modelFun=@(p,x) (p(1)/pi).*atan((x+p(2))/p(3))+p(4); % P1= sliprate, p2= horizontal offset, p3=locking, p4=vertical offset
starting=[-14,-13,12.2,-6];

Dsort=bins(1):1:bins(end);

ax=gobjects(3,col);

for ii=1:nprof
    
    coefEsts{(ii-1)*2+1}=nlinfit(dist(in{ii})',modGPS(1,in{ii}),modelFun,starting);
    coefEsts{(ii-1)*2+2}=nlinfit(dist(in{ii})',modGPS(2,in{ii}),modelFun,starting);
    
    maxdX=modelFun(coefEsts{ii*2-1},Dsort(1:end-1))-modelFun(coefEsts{ii*2},Dsort(2:end));
    maxdX=find(abs(maxdX)==max(abs(maxdX)));
    maxdY=modelFun(coefEsts{ii*2-1},Dsort(1:end-1))-modelFun(coefEsts{ii*2},Dsort(2:end));
    maxdY=find(abs(maxdY)==max(abs(maxdY)));
    
    maxshearloc((ii-1)*2+[1:2])=[Dsort(maxdX),Dsort(maxdY)];
    
    
    ax(ii*3-2)=subplot(3,col,ii)
    plot(dist(in{ii}),modGPS(1,in{ii}),'.','Color',[0.4 0.4 0.4])
    hold on
    plot(sort(dist),modelFun(coefEsts{(ii-1)*2+1},sort(dist)),'--','Color',colors(ii),'LineWidth',1.5)
    title(['Profile ',num2str(prof_nums(ii))])
    ylabel('Para Rate (mm/yr)')
    xline(0,'--');
    legend('Data','Simple Arctan Fit','Fault Trace');
    ylim(paray)
    xlim([-20 120])
    
    ax(ii*3-1)=subplot(3,col,ii+col)
    plot(dist(in{ii}),modGPS(2,in{ii}),'.','Color',[0.4 0.4 0.4])
    hold on
    plot(sort(dist),modelFun(coefEsts{(ii-1)*2+2},sort(dist)),'--','Color',colors(ii),'LineWidth',1.5)
    ylabel('Perp Rate (mm/yr)')
    xline(0,'--');
    legend('Data','Simple Arctan Fit','Fault Trace');
    ylim(perpy)
    xlim([-20 120])
    
    ax(ii*3)=subplot(3,col,ii+2*col)
    plot(dist(in{ii}),modGPS(3,in{ii}),'.','Color',[0.4 0.4 0.4])
    hold on
    errorbar(bin_center,bindata{ii}(5,:),bindata{ii}(6,:),'k','LineWidth',3,'Color',colors(ii))
    ylabel('Uplift Rate (mm/yr)')
    xlabel('Distance from Fault (km)')
    xline(0,'--');
    legend('Data','Mean Trend','Fault Trace');
    ylim(verty)
    xlim([-20 120])
end

if nprof>1
    legend_lab=[];
    legend_para=[];
    legend_perp=[];
    legend_ver=[];
    
    for ii=1:size(prof_nums,2)
        legend_lab=[legend_lab,{sprintf('Profile %.0f',prof_nums(ii))}];
        legend_para=[legend_para,{sprintf('Profile %.0f',prof_nums(ii))},{sprintf('Max Shear: %.1f km',-coefEsts{ii*2-1}(2))}];
        legend_perp=[legend_perp,{sprintf('Profile %.0f',prof_nums(ii))},{sprintf('Max Shear: %.1f km',-coefEsts{ii*2}(2))}];
        legend_ver=[legend_ver,{sprintf('Profile %.0f',prof_nums(ii))},{sprintf('Max Mean Uplift: %.1f mm/yr',max(bindata{ii}(5,:)))},{sprintf('Max Uplift: %.1f mm/yr',max((modGPS(3,in{ii}))))}];
    end
    
    subplot(3,col,col)
    title('All Profiles')
    hold on
    for ii=1:nprof
        plot(sort(dist),modelFun(coefEsts{(ii-1)*2+1},sort(dist)),'--','Color',colors(ii),'LineWidth',2)
        xline(-coefEsts{ii*2-1}(2),':','Color',colors(ii),'LineWidth',2);
    end
    ylabel('Para Rate (mm/yr)')
    xline(0,'k--');
    ylim(paray)
    xlim([-20 120])
    legend(legend_para);
    
    subplot(3,col,col*2)
    hold on
    for ii=1:nprof
        plot(sort(dist),modelFun(coefEsts{(ii-1)*2+2},sort(dist)),'--','Color',colors(ii),'LineWidth',2)
        xline(-coefEsts{ii*2}(2),':','Color',colors(ii),'LineWidth',2);
    end
    ylabel('Perp Rate (mm/yr)')
    xline(0,'k--');
    ylim(perpy)
    xlim([-20 120])
    legend(legend_perp);
    
    subplot(3,col,col*3)
    hold on
    for ii=1:nprof
        errorbar(bin_center,bindata{ii}(5,:),bindata{ii}(6,:),'Color',colors(ii),'LineWidth',3)
        xline(bin_center(find(bindata{ii}(5,:)==max((bindata{ii}(5,:))))),':','Color',colors(ii),'LineWidth',2);
        d=dist(in{ii});
        xline(d(find(modGPS(3,in{ii})==max((modGPS(3,in{ii}))))),'.-','Color',colors(ii),'LineWidth',2);
    end
    ylabel('Uplift Rate (mm/yr)')
    xlabel('Distance from Fault (km)')
    xline(0,'k--');
    ylim(verty)
    xlim([-20 120])
    legend(legend_ver);
end

