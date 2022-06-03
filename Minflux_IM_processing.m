%% Created by mbossi
% Mariano.Bossi@mr.mpg.de
% Last update 03.06.22
%%%%%%%%%%%%%

cfr   = itr.cfr;
efo   = itr.efo;
dcr   = itr.dcr;
loc   = itr.loc(vld,5,:);
itr2  = itr.itr(vld,:);
dataT = zeros(size(vld,2),16);
for k = 1:1:size(vld,2)
    dataT(k,1) = k;
    dataT(k,2) = single(tid(k));  
    dataT(k,3) = single(tim(k)*1e3);
    loc_k      = single(loc(k,1,:)*1e9);
    dataT(k,4) = loc_k(1);
    dataT(k,5) = loc_k(2);
    dataT(k,6) = loc_k(3);
    dataT(k,7) = single(itr2(k,end));
    dataT(k,8) = single(cfr(k,end));        % cfr
    dataT(k,9) = single(efo(k,end));        % efo
    dataT(k,10)= single(dcr(k,end));        % dcr
end
eco = itr.eco;
ecc = itr.ecc;
ecoS = sum(eco,2);
eccS = sum(ecc,2);
ecoE = double(eco(:,end));
eccE = double(ecc(:,end));
cfrE = double(cfr(:,end));
%%
tMIN = 1000*5*60;
dataTc = dataT(dataT(:,3)>=tMIN,:);
locs_pull = 4;
%%
SMid   = unique(dataTc(:,2));
dataTm = zeros(length(SMid),16);
dataTD = zeros(size(dataTc,1),2);
%cols = {'id'   'tid'  'tim'   'posX'   'posY'  'posZ'   'itr'   'cfr'   'efo' 'dcr' 'uncX' 'uncY' 'uncXY' 'locs' 'isnoise' 'IDX'}
for k = 1:length(SMid)
    SMid_k2        = SMid(k);
    dataTm_k       = dataTc(dataTc(:,2)==SMid_k2,:);
    dataTm(k,:)    = mean(dataTm_k,1);
    pos    = dataTm_k(:, 4:5);
    uncXY  = std(pos,0,1);
    uncXYM = sqrt(uncXY(1).^2 + uncXY(2).^2);
    dataTm(k,11:12)= uncXY;
    dataTm(k,13)   = uncXYM;
    dataTm(k,14)   = size(dataTm_k,1);
    if size(dataTm_k,1)>1
        epsilon = 4;
        MinPts  = 6;
        [IDX] = dbscan([dataTm_k(:,4),dataTm_k(:,5)],epsilon,MinPts);
        isnoise = IDX==-1;
        ixx = find(dataTc(:,2)==SMid_k2);
        dataTD(ixx,:) = [isnoise, IDX];
    end
end
dataTn          = dataTc;
dataTn(:,15:16) = dataTD;
dataTdbsc       = dataTn(dataTn(:,15)==0,:); 
xErr = dataTm(dataTm(:,14)<=locs_pull,11);
xErr = nonzeros(dataTm(:,11));
figure(1);subplot(2,1,1); hist(xErr,80); title('xErr'); axis tight; grid on
yErr = dataTm(dataTm(:,14)<=locs_pull,12);
yErr = nonzeros(dataTm(:,12));
subplot(2,1,2); hist(yErr,80); title('yErr'); axis tight; grid on

%%
SMid   = unique(dataTdbsc(:,2));
dist_r = [];
for k = 1:length(SMid) % first loop selects localizations by 'tid'
    SMid_k   = SMid(k);
    dataT_k  = dataTdbsc(dataTdbsc(:,2)==SMid_k,:);
    SMidDB   = unique(dataT_k(:,16));
    dataT_kk = zeros(size(SMidDB,1),16);
    for kk = 1:size(SMidDB,1)   % second loop goes inside each DBSCAN cluster
        SMid_kk        = SMidDB(kk);
        dataT_kDB      = dataT_k(dataT_k(:,16)==SMid_kk,:);
        dataT_kk(kk,:) = mean(dataT_kDB,1);
        pos   = dataT_kDB(:, 4:5);
        dist_kk = [];
        if size(pos,1)>locs_pull
            averXY    = mean(pos);
            dist_kk(:,1)  = pos(:,1) - averXY(1);
            dist_kk(:,2)  = pos(:,2) - averXY(2);
            dist_r = [dist_r; dist_kk];
        end
    end
end
[countsX, centersX] = hist(dist_r(:,1),200);
[countsY, centersY] = hist(dist_r(:,2),200);
countsX = countsX./max(countsX); countsY = countsY./max(countsY);
Gfit_objX    = fit(centersX',countsX','gauss1')
Gfit_coeffsX = coeffvalues(Gfit_objX);
Gfit_objY    = fit(centersY',countsY','gauss1')
Gfit_coeffsY = coeffvalues(Gfit_objY);
figure(2); plot(Gfit_objX,centersX,countsX); hold on
plot(Gfit_objY,centersY,countsY);
sig_m = mean([Gfit_coeffsX(3)/sqrt(2) Gfit_coeffsY(3)/sqrt(2)]);
axis tight; grid on; title(['<sigma> = ' num2str(sig_m)]); hold off; drawnow
%%
NL     = 10;
delta  = 1.644854*sig_m;
SMid   = unique(dataTdbsc(:,2));
dataTf = []; dist_f =[];
for k = 1:length(SMid)
    SMid_k    = SMid(k);
    dataTf_k  = dataTdbsc(dataTdbsc(:,2)==SMid_k,:);
    SMidDB    = unique(dataTf_k(:,16));
    for kk = 1:size(SMidDB,1)
        SMid_kk   = SMidDB(kk);
        dataTf_kk = dataTf_k(dataTf_k(:,16)==SMid_kk,:);
        if size(dataTf_kk,1)> NL
            pos_kk  = dataTf_kk(:,4:5);
            dist_kk = pdist2(pos_kk,mean(pos_kk));
            [dist_kkS,ixs] = sort(dist_kk);
            dataTf_kk = dataTf_kk(ixs,:);
            flag = size(find(dist_kkS<delta),1);
            dist_kk = [];
            if flag>0
                dataTf = [dataTf; dataTf_kk(1:flag,:)];
                posA   = dataTf_kk(1:flag,4:5);
                averA = mean(posA);
                dist_kk(:,1)  = posA(:,1) - averA(1);
                dist_kk(:,2)  = posA(:,2) - averA(2);
                dist_f = [dist_f; dist_kk];
            else
                dataTf = [dataTf; dataTf_kk(1:NL,:)];
                posA  = dataTf_kk(1:NL,4:5);
                averA = mean(posA);
                dist_kk(:,1)  = posA(:,1) - averA(1);
                dist_kk(:,2)  = posA(:,2) - averA(2);
                dist_f = [dist_f; dist_kk];
            end
        end
    end
end
%%
[countsX centersX] = hist(dist_f(:,1),200);
[countsY centersY] = hist(dist_f(:,2),200);
countsX = countsX./max(countsX); countsY = countsY./max(countsY);
Gfit_objX    = fit(centersX',countsX','gauss1');
Gfit_coeffsX = coeffvalues(Gfit_objX);
Gfit_objY    = fit(centersY',countsY','gauss1');
Gfit_coeffsY = coeffvalues(Gfit_objY);
sig_F = mean([Gfit_coeffsX(3)/sqrt(2) Gfit_coeffsY(3)/sqrt(2)])
szLS = 10;
pxLS = 0.5;
Loc_spread = hist3(dist_f, {-szLS:pxLS:szLS -szLS:pxLS:szLS});
figure(3); imagesc([-szLS +szLS],[-szLS +szLS],Loc_spread);axis image; colormap hot; title('dataTf:  Localization spread');
cl = get(get(3,'Children'),'CLim');
viscircles([0 0], sig_F,'Color','b','LineStyle',':','LineWidth',0.1);
viscircles([0 0], 2*sig_F,'Color','b','LineStyle','--','LineWidth',0.1);
hold off
%%
usedLOCS = unique(dataTf(:,2));
photonsF = [];
locTf    = [];
for k = 1:size(usedLOCS,1)
    usedLOCS(k);
    TID_k  = dataTf(dataTf(:,2)==usedLOCS(k),:);
    SMidDB = unique(TID_k(:,16));
    for kk = 1:size(SMidDB,1)
        SMid_kk   = SMidDB(kk);
        TID_kk = TID_k(TID_k(:,16)==SMid_kk,:);
        repet = size(TID_kk,1);
        ixs = TID_kk(:,1);
        eco_k8 = eco(ixs,:);
        eco_k8 = sum(eco_k8(:));    % sum(eco_k8(:,5));
        ecc_k8 = ecc(ixs,:);
        ecc_k8 = sum(ecc_k8(:));    % sum(ecc_k8(:,5));
        photonsF = [photonsF; (ecc_k8 + eco_k8)];
        locTf     = [locTf;     repet];
    end
end
loc    = mean(locTf)
ph_loc = mean(photonsF)/mean(locTf)
%%
[image1, IMsize1] = FN_gaussian_rendering(dataTf,1,0,sig_F); 
figure(4);imagesc(sqrt(image1)); axis image; colormap hot;

%%
%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = D2GaussFunction(x,xdata)
%   x0 = [Amp, x0, wx, y0, wy, offset]
F = x(1)*exp(   -((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) )    )+x(6);
end

function [image, IMsize] = FN_gaussian_rendering(pos,px,IMsize,gauss_fix_size)
pos(:,4) = pos(:,4) - min(pos(:,4));
pos(:,5) = pos(:,5) - min(pos(:,5));
if IMsize == 0
    IMsize(1) = ceil(max(pos(:,4)));
    IMsize(2) = ceil(max(pos(:,5)));
end
pos(:,4) = pos(:,4)./px;
pos(:,5) = pos(:,5)./px;
IMsize   = ceil(IMsize);
if ~gauss_fix_size
    sigma = mean(mean(pos(:,11:12)));
else
    sigma = gauss_fix_size;
end
sigma    = sigma/px;
grid_sz  = round(5*sigma/px);
pos(:,4) = pos(:,4) + 4*grid_sz./px;
pos(:,5) = pos(:,5) + 4*grid_sz./px;
IMsize   = IMsize + 8*grid_sz;
SRsize   = ceil(IMsize./px); 
SRimage  = zeros(SRsize(2),SRsize(1));
si1 = size(pos,1);
top = (grid_sz+1);
pos = pos(pos(:,4)>top & pos(:,4)<SRsize(1)-top,:);
pos = pos(pos(:,5)>top & pos(:,5)<SRsize(2)-top,:);
si2 = size(pos,1);
box_size = 1+2*grid_sz(1);
[X,Y]    = meshgrid(-grid_sz:1:grid_sz);
xdata    = cat(3,X,Y);
rows_number = size(pos,1);
wbh = waitbar(0,'rendering');
for k=1:rows_number
    osx   = floor(pos(k,4));
    osy   = floor(pos(k,5));
    x1    = pos(k,4) - osx;
    y1    = pos(k,5) - osy;
    SMgauss = 1 * exp(   -((xdata(:,:,1)-x1).^2/(2*(sigma^2)) + (xdata(:,:,2)-y1).^2/(2*(sigma^2)) )    );
    SRimage(osy-grid_sz:osy+grid_sz,osx-grid_sz:osx+grid_sz) = SRimage(osy-grid_sz:osy+grid_sz,osx-grid_sz:osx+grid_sz) + SMgauss;
    waitbar(k/rows_number,wbh,'rendering');
end
close(wbh)
grid_szSR = grid_sz./px;
SRimage = SRimage(round(4*grid_szSR)+1: end-(round(4*grid_szSR)) , round(4*grid_szSR)+1: end-(round(4*grid_szSR)));
IMsize = flip(size(SRimage).*px);
image = SRimage;
end
