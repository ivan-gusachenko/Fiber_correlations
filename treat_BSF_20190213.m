path = '/Users/ivangusachenko/DATA/LKB/FIBER_Z_CORRELATION/20190213-STEP_58mm';

%% get center drift
load([path '/SP_Z_BF57801_bf.mat'])
sz = size(SPECKLE_Z_BF);
x = 1:sz(1);
y = 1:sz(2);

[Y, X] = meshgrid(x, y);
XX = repmat(X, [1 1 sz(3)]);
YY = repmat(Y, [1 1 sz(3)]);

a = -0.22;
b = 0.33;
x0 = 226;
y0 = 200;
zi0 = round(sz(3)/2);

for i = 1:sz(3)
    % i = 5;
    X0 = a*(i-zi0)+x0;
    Y0 = b*(i-zi0)+y0;
    imagesc(SPECKLE_Z_BF(:,:,i));
    hold on
    plot(X0, Y0, 'or', 'MarkerSize', 20);
    hold off
    pause(0.1)
end

X0 = a*((1:sz(3))-zi0)+x0;
Y0 = b*((1:sz(3))-zi0)+y0;

margin = ceil(max(abs(X0(end)-X0(1)), abs(Y0(end)-Y0(1))))+2;

%% transform data with the drift
load([path '/SP_Z_L55474_focus.mat']);
SPECKLE_f = SPECKLE;
load([path '/SP_Z_L57251_speckle.mat']);
SPECKLE_s = SPECKLE;
clear SPECKLE


sz = size(SPECKLE_s);


%% check that drift looks ok
for i = 1:sz(3)
    % i = 5;
    X0 = a*(i-zi0)+x0;
    Y0 = b*(i-zi0)+y0;
    imagesc(SPs_crop_aligned(:,:,i, 43));
%     hold on
%     plot(X0, Y0, 'or', 'MarkerSize', 20);
%     hold off
    pause(0.1)
end

%%
[Yq0, Xq0] = meshgrid(x(margin:end-margin), y(margin:end-margin));
SPs_crop_aligned = zeros(sz(1)-2*margin+1, sz(2)-2*margin+1, sz(3), sz(4));
SPf_crop_aligned = zeros(sz(1)-2*margin+1, sz(2)-2*margin+1, sz(3), sz(4));

for zi = 1:sz(3)
    Xq = Xq0+X0(zi)-x0;
    Yq = Yq0+Y0(zi)-y0;
    for li = 1:sz(4)
        im_s = interp2(SPECKLE_s(:,:,zi, li), Xq, Yq);
        im_f = interp2(SPECKLE_f(:,:,zi, li), Xq, Yq);
        SPs_crop_aligned(:,:,zi,li) = im_s;
        SPf_crop_aligned(:,:,zi,li) = im_f;
%         imagesc(im_s);
%         pause(0.1);
    end
    disp(zi)
end

clear SPECKLE_s 
clear SPECKLE_f

Z = meta_Z_L.Z;
Z = (Z-mean(Z))*1000;

%% calculate correlations for speckle

COR_s = zeros(sz(3), sz(4));
COR_f = zeros(sz(3), sz(4));

REF_s = SPs_crop_aligned(:,:,round(sz(3)/2), round(sz(4)/2));
REF_f = SPf_crop_aligned(:,:,round(sz(3)/2), round(sz(4)/2));


for zi = 1:sz(3) 
    for li = 1:sz(4)
        COR_s(zi, li) = sum(sum(SPs_crop_aligned(:,:,zi,li).*REF_s))/sqrt(sum(sum(REF_s.^2)))/sqrt(sum(sum(SPs_crop_aligned(:,:,zi,li).^2)));
        COR_f(zi, li) = sum(sum(SPf_crop_aligned(:,:,zi,li).*REF_f))/sqrt(sum(sum(REF_f.^2)))/sqrt(sum(sum(SPf_crop_aligned(:,:,zi,li).^2)));
    end
end
figure(1)
imagesc(Z, meta_Z_L.Lambda, COR_s)
title('Speckle correlation')
figure(2)
imagesc(Z, meta_Z_L.Lambda, COR_f)
title('Focus correlation')




%% calculate max for focus

MAX = zeros(sz(3), sz(4));

for zi = 1:sz(3) 
    for li = 1:sz(4)
        im = imresize(SPf_crop_aligned(:,:,zi,li), 1/4, 'box');
        MAX(zi, li) = max(im(:));
    end
end
figure(3)
imagesc(Z, meta_Z_L.Lambda, MAX)
title('Focus maximum')

%%
save([path, '/RESULTS20190213.mat'], 'COR_s', 'COR_f', 'MAX');

