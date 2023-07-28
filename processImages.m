function [K2_fundamental,K2_total,K2_shot,...
    K2_read,K2_quantized,K2_spatial,cc,cc_fit] ...
    = processImages(imgs,mask,gain,dark_var,mean_img,mean_img_num,use_intensity_weight,use_lsq_fit,window_size,consider_spatial_heterogeneity)
%processImages Process image to get corrected K2.
%   imgs = 3D double (y,x,time). Only provide images within ROI to be
%   processed.
%   gain = matrix of gain map. ADU/e-
%   dark_var = dark image temporal variance. Should already have quantization
%   noise variance subtracted.
%   mean_img = 2D double (y,x)
%   use_intensity_weight = Boolean. If true, weighs K2 values by their
%   intensities.

frame_num = size(imgs,3);

x_window_start = 1:window_size:size(imgs,2);
x_window_start(size(imgs,2)-x_window_start < window_size - 1) = []; % only if there is enough buffer
y_window_start = 1:window_size:size(imgs,1);
y_window_start(size(imgs,1)-y_window_start < window_size - 1) = [];

K2_total_window = nan(length(y_window_start),length(x_window_start),frame_num);
K2_shot_window = K2_total_window; K2_read_window = K2_total_window; K2_quantization_window = K2_total_window;
K2_spatial_window = K2_total_window; cc_window = K2_total_window; cc_window_fit = K2_total_window;
cc_window_norm = K2_total_window;

good_window = false(length(y_window_start),length(x_window_start));

cc_mean = squeeze(mean(mean(imgs,1),2));
cc_mean_smooth = smooth(cc_mean); % the temporal pattern that is expected

% fitting function
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','off');
lb = [];
ub = [];
x0=[0,0];
fun = @(x,Ipattern_ts)(x(1).*Ipattern_ts+x(2));

% intensity weight
intensity_weight = mean_img./mean(mean_img(:));

% for each window
for ii=1:length(y_window_start)
    for jj=1:length(x_window_start)
        x_window = x_window_start(jj):x_window_start(jj)+window_size-1;
        y_window = y_window_start(ii):y_window_start(ii)+window_size-1;

        imgs_window = imgs(y_window,x_window,:);
        mask_window = mask(y_window,x_window);
        mask_window = mask_window(:); % extend to column
        %         var_diff_window_map = mean(var_diff(:,x_window),2);
        good_window(ii,jj) = sum(mask_window)./numel(mask_window) >= 0.5;

        if good_window(ii,jj)
            imgs_window_masked = reshape(imgs_window,[],size(imgs_window,3)); % also extend to column
            imgs_window_masked = imgs_window_masked(mask_window,:);
            cc_window_masked = mean(imgs_window_masked,1);

            if use_lsq_fit
                x = lsqcurvefit(fun,x0,cc_mean_smooth(:),cc_window_masked',lb,ub,options);
                cc_window_masked_fit = x(1).*cc_mean_smooth(:)+x(2); % the rescaled pattern for this patch
            else
                cc_window_masked_fit = cc_window_masked';
            end

            % get average gain for the window
            gain_window = gain(y_window,x_window);
            gain_window = nanmean(gain_window(:));

            dark_var_window = mean(mean(dark_var(y_window,x_window)));

            % spatial heterogeneity
            mean_img_window = mean_img(y_window,x_window);
            var_spatial = var(mean_img_window(:));
            mean_spatial = mean(mean_img_window(:));
            var_spatial = var_spatial - gain_window*mean_spatial/mean_img_num; % remove other noise contributions from mean img

            if ~consider_spatial_heterogeneity
                var_spatial = 0;
            end

            K2_total_window(ii,jj,:) = var(imgs_window_masked,0,1)'./(cc_window_masked_fit.^2);
            K2_shot_window(ii,jj,:) = shotNoiseK2(cc_window_masked_fit,gain_window);
            K2_read_window(ii,jj,:) = dark_var_window./(cc_window_masked_fit.^2);
            K2_quantization_window(ii,jj,:) = 1./12./(cc_window_masked_fit.^2);
            K2_spatial_window(ii,jj,:) = var_spatial./(mean_spatial^2)*ones(size(cc_window_masked_fit));

            cc_window(ii,jj,:) = cc_window_masked;
            cc_window_fit(ii,jj,:) = cc_window_masked_fit;
            cc_window_norm(ii,jj,:) = mean(mean(intensity_weight(y_window,x_window),1),2);

            if use_intensity_weight
                K2_total_window(ii,jj,:) = ...
                    squeeze(mean(mean(K2_total_window(ii,jj,:).*cc_window_norm(ii,jj,:).^2,1),2));
                K2_shot_window(ii,jj,:) = ...
                    squeeze(mean(mean(K2_shot_window(ii,jj,:).*cc_window_norm(ii,jj,:).^2,1),2));
                K2_read_window(ii,jj,:) = ...
                    squeeze(mean(mean(K2_read_window(ii,jj,:).*cc_window_norm(ii,jj,:).^2,1),2));
                K2_quantization_window(ii,jj,:) = ...
                    squeeze(mean(mean(K2_quantization_window(ii,jj,:).*cc_window_norm(ii,jj,:).^2,1),2));
                K2_spatial_window(ii,jj,:) = ...
                    squeeze(mean(mean(K2_spatial_window(ii,jj,:).*cc_window_norm(ii,jj,:).^2,1),2));
            end
        else
            K2_total_window(ii,jj,:) = nan;
            K2_shot_window(ii,jj,:) = nan;
            K2_read_window(ii,jj,:) = nan;
            K2_quantization_window(ii,jj,:) = nan;
            K2_spatial_window(ii,jj,:) = nan;
            cc_window(ii,jj,:) = nan;
            cc_window_fit(ii,jj,:) = nan;
            cc_window_norm(ii,jj,:) = nan;
        end
    end
end

K2_fundamental_window = K2_total_window - K2_shot_window - K2_read_window - K2_quantization_window - K2_spatial_window;

% only select good roi
K2_total_window = reshape(K2_total_window,[],size(K2_total_window,3));
K2_total_window = K2_total_window(good_window(:),:);
K2_shot_window = reshape(K2_shot_window,[],size(K2_shot_window,3));
K2_shot_window = K2_shot_window(good_window(:),:);
K2_read_window = reshape(K2_read_window,[],size(K2_read_window,3));
K2_read_window = K2_read_window(good_window(:),:);
K2_quantization_window = reshape(K2_quantization_window,[],size(K2_quantization_window,3));
K2_quantization_window = K2_quantization_window(good_window(:),:);
K2_spatial_window = reshape(K2_spatial_window,[],size(K2_spatial_window,3));
K2_spatial_window = K2_spatial_window(good_window(:),:);
K2_fundamental_window = reshape(K2_fundamental_window,[],size(K2_fundamental_window,3));
K2_fundamental_window = K2_fundamental_window(good_window(:),:);

cc_window = reshape(cc_window,[],size(cc_window,3));
cc_window = cc_window(good_window(:),:);
cc_window_fit = reshape(cc_window_fit,[],size(cc_window_fit,3));
cc_window_fit = cc_window_fit(good_window(:),:);
cc_window_norm = reshape(cc_window_norm,[],size(cc_window_norm,3));
cc_window_norm = cc_window_norm(good_window(:),:);

% get mean
K2_total = nanmean(K2_total_window,1)';
K2_shot = nanmean(K2_shot_window,1)';
K2_read = nanmean(K2_read_window,1)';
K2_quantized = nanmean(K2_quantization_window,1)';
K2_spatial = nanmean(K2_spatial_window,1)';
K2_fundamental = nanmean(K2_fundamental_window,1)';

cc = nanmean(cc_window,1)';
cc_fit = nanmean(cc_window_fit,1)';
cc_ROI_norm = nanmean(cc_window_norm,1)';
end