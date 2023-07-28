clear
close all;
clc;

[~, matlab_ver] = version;
matlab_ver = str2double(matlab_ver(end-3:end));
use_pop_ups = false; % change this to True if you want to use UI file selector to select files.

%% meta tags

meta.aux_dir = [];
meta.img_dir = '.\sample data\img_file.h5';
meta.use_intensity_weight = true;
meta.use_lsq_fit = false;

meta.mean_img_dir = [];

meta.dark_dir = '.\sample data\dark.mat';
meta.gain_map_dir = '.\resources\gain.mat';
meta.mask_dir = '.\resources\mask.mat';

%% interaction to decide files and other stuff

if use_pop_ups
    title_str = 'Get aux file.';
    if ~ispc; menu(title_str,'OK'); end
    [aux_file, aux_folder] = uigetfile(fullfile(fileparts(meta.aux_dir),'*.h5'),title_str);
    if isnumeric(aux_file)
        meta.aux_dir = []; clear aux_folder aux_file
    else
        meta.aux_dir = fullfile(aux_folder,aux_file); clear aux_folder aux_file
    end
    title_str = 'Get data file.';
    if ~ispc; menu(title_str,'OK'); end
    [img_file, img_folder] = uigetfile(fullfile(fileparts(meta.img_dir),'*_506*.h5'),title_str);
    meta.img_dir = fullfile(img_folder,img_file);

    % mean image file
    title_str = 'Choose mean image file. If not chosen, the recorded data will be used.';
    if ~ispc; menu(title_str,'OK'); end
    if ~isempty(meta.mean_img_dir)
        [mean_img_file, mean_img_folder] = uigetfile(fullfile(fileparts(meta.mean_img_dir),'*.mat'),...
            title_str);
    else
        [mean_img_file, mean_img_folder] = uigetfile('.',...
            title_str);
    end
    if isnumeric(mean_img_file)
        meta.mean_img_dir = [];
    else
        meta.mean_img_dir = fullfile(mean_img_folder,mean_img_file);
    end

    % dark img
    title_str = 'Get dark image and variance. H5 files can be used as well.';
    if ~ispc; menu(title_str,'OK'); end
    [dark_file, dark_folder] = uigetfile(fullfile(fileparts(meta.dark_dir),'*.mat'),title_str);
    meta.dark_dir = fullfile(dark_folder,dark_file);

    % mask
    title_str = 'Get mask file. If not chosen, then default mask will be used.';
    if ~ispc; menu(title_str,'OK'); end
    [mask_file, mask_folder] = uigetfile(fullfile(fileparts(meta.mask_dir),'*.mat'),title_str);
    if isnumeric(mask_file)
        clear mask_folder mask_file
    else
        meta.mask_dir = fullfile(mask_folder,mask_file); clear mask_folder mask_file
    end

    % load gain map
    [gain_map_file, gain_map_folder] = uigetfile(fullfile(fileparts(meta.gain_map_dir),'*.mat'),'Get gain map');
    meta.gain_map_dir = fullfile(gain_map_folder,gain_map_file);
else
end

[~,img_file] = fileparts(meta.img_dir);
record_time = img_file(12:end-5);

c = clock();
save_file = ['analyzed_K2_' img_file '_' ...
    num2str(c(1)) '_' num2str(c(2)) '_' num2str(c(3)) '_' num2str(c(4)) '_' ...
    num2str(c(5)) '.mat'];
t = double(h5read(meta.img_dir,'/timestamp'))*1e-9;
fr = round(1000/(t(2)-t(1)))/1000;

if contains(meta.dark_dir,'.h5')
    % get dark images and save it
    disp('Obtaining the dark image from h5 file.')
    dark_meta = h5read(meta.dark_dir,'/Metadata');
    dark_length = numel(dark_meta.timestamp);
    dark_var = zeros(2304); dark_img = zeros(2304);
    for block_start = 1:100:dark_length
        disp(['Dark frame # ' num2str(block_start)]);
        block_end = min([block_start+100-1 dark_length]);
        dark_imgs = double(readHamamatsuH5(meta.dark_dir,[block_start block_end]));
        dark_var = dark_var + var(dark_imgs,0,3).*(block_end - block_start + 1);
        dark_img = dark_img + mean(dark_imgs,3).*(block_end - block_start + 1); clear dark_imgs;
    end
    dark_img = dark_img./dark_length;
    dark_var = dark_var./dark_length;
    meanIDark = dark_img; varIDark = dark_var;
    [dark_folder, dark_file] = fileparts(meta.dark_dir);
    save(fullfile(dark_folder,[dark_file '_dark.mat']),'varIDark','meanIDark');
    dark_var = dark_var - 1/12;
else
    dark_data = load(meta.dark_dir);
    dark_var = dark_data.varIDark - 1/12;
    dark_img = dark_data.meanIDark;
end

load(meta.mask_dir);
load(meta.gain_map_dir);

% disp
disp(meta)

% remove pixels with high variance
mask = mask & dark_var < 100;

%% other parameters

block_size = 100;
window_size = 7;

block_start = 1:block_size:numel(t);

%% get mean image

% get mean time course and select good time indices
disp('get mean time course');
cc_mean = [];
for block=1:10:numel(t)-1
    img = double(readHamamatsuH5(meta.img_dir,[block block])) - dark_img;
    cc_mean = [cc_mean mean(img(:))];
end
cc_mean = cc_mean';

f_cc_mean = figure('Position',[50 50 500 400]);
plot(10:10:10*numel(cc_mean),cc_mean);
good_ind = true(size(t));

% get average image
disp('get average image')

% first quickly get average image to show
good_start = find(good_ind,1,'first');
mean_img = zeros(size(img,1),size(img,2));
block_starts = good_start:10:find(good_ind,1,'last');
for block_ind = block_starts
    img = double(readHamamatsuH5(meta.img_dir,[block_ind block_ind])) - dark_img;
    mean_img = mean_img + img;
end
mean_img = mean_img./numel(block_starts);

f1 = figure('Position',[50 50 400 300]); imagesc(mean_img); axis image; colorbar;

% use loaded image if told to do so
if ~isempty(meta.mean_img_dir)
    mean_img_data = load(meta.mean_img_dir);
    mean_img_reference = mean_img_data.mean_img;
    mean_img_reference_frame_num = sum(mean_img_data.good_ind);
else
    % if not read through all good indices and create the mean image
    mean_img = zeros(size(mean_img));
    for block_ind = good_start:block_size:find(good_ind,1,'last')
        if mod(block_ind,1000) < 100
            disp(['Getting image from frame # ' num2str(block_ind)]);
        end
        block_length = min([block_size find(good_ind,1,'last')-block_ind+1]);
        imgs = double(readHamamatsuH5(meta.img_dir,[block_ind block_ind+block_length-1]));
        imgs = double(imgs) - repmat(dark_img,[1,1,block_length]);
        mean_img = mean_img + sum(imgs,3);
    end
    mean_img = mean_img./sum(good_ind);
    mean_img_reference = mean_img;
    mean_img_reference_frame_num = sum(good_ind);
    save(fullfile([meta.img_dir(1:end-3) '_mean_img.mat']),'mean_img','good_ind');
end

%% roi

roi_positions(:,1) = [1 1 200 200];

set(0,'CurrentFigure',f1);
for roi_ind = 1:size(roi_positions,2)
    rectangle('Position',roi_positions(:,roi_ind));
end
pause(0.01);

%% get K2

block_size = 100;
disp([num2str(length(block_start)-1) ' blocks']);
K2_total_ROI = nan((length(block_start)-1)*block_size,size(roi_positions,2));
K2_fundamental_ROI = K2_total_ROI; K2_shot_ROI = K2_total_ROI;
K2_read_ROI = K2_total_ROI; K2_quantized_ROI = K2_total_ROI;
K2_spatial_ROI = K2_total_ROI;
cc_ROI = K2_total_ROI; cc_ROI_fit = K2_total_ROI;
for block=1:length(block_start)
    if mod(block,5) == 1
        disp(['Block # ' num2str(block)]);
    end

    time_ind = (block-1)*block_size + 1: block*block_size;
    imgs = readHamamatsuH5(meta.img_dir,[time_ind(1) time_ind(end)]);
    imgs = double(imgs)-repmat(dark_img,[1,1,block_size]);

    for roi_ind = 1:size(roi_positions,2)
        if mod(roi_ind,1e4) == 1
            disp(['ROI # ' num2str(roi_ind) ' : ' char(datetime)]);
        end
        
        x_start_ROI = round(roi_positions(1,roi_ind));
        x_end_ROI = round(roi_positions(1,roi_ind) + roi_positions(3,roi_ind) - 1);
        y_start_ROI = round(roi_positions(2,roi_ind));
        y_end_ROI = round(roi_positions(2,roi_ind) + roi_positions(4,roi_ind) - 1);
        mask_ROI = mask(y_start_ROI:y_end_ROI,x_start_ROI:x_end_ROI);
        if sum(mask_ROI) == 0
            K2_fundamental_ROI(time_ind,roi_ind) = nan;
            K2_total_ROI(time_ind,roi_ind) = nan;
            K2_shot_ROI(time_ind,roi_ind) = nan;
            K2_read_ROI(time_ind,roi_ind) = nan;
            K2_quantized_ROI(time_ind,roi_ind) = nan;
            K2_spatial_ROI(time_ind,roi_ind) = nan;
            cc_ROI(time_ind,roi_ind) = nan;
        else
            imgs_ROI = imgs(y_start_ROI:y_end_ROI,x_start_ROI:x_end_ROI,:);
            gain_ROI = gain(y_start_ROI:y_end_ROI,x_start_ROI:x_end_ROI);
            dark_var_ROI = dark_var(y_start_ROI:y_end_ROI,x_start_ROI:x_end_ROI);
            mean_img_ROI = mean_img_reference(y_start_ROI:y_end_ROI,x_start_ROI:x_end_ROI);

            [K2_fundamental_ROI(time_ind,roi_ind),K2_total_ROI(time_ind,roi_ind),...
                K2_shot_ROI(time_ind,roi_ind),K2_read_ROI(time_ind,roi_ind),...
                K2_quantized_ROI(time_ind,roi_ind),K2_spatial_ROI(time_ind,roi_ind),...
                cc_ROI(time_ind,roi_ind)] = ...
                processImages(imgs_ROI,mask_ROI,gain_ROI,dark_var_ROI,mean_img_ROI,mean_img_reference_frame_num,...
                meta.use_intensity_weight,meta.use_lsq_fit,window_size,true);
        end
    end
end

%% get aux data

try
    stim_data=h5read(meta.aux_dir,'/DigitalIn_headers');
    stim_time=double(stim_data.timestamp);
    stim_time=stim_time(stim_time~=0)*1e-9;  % in s
catch
    stim_time = [];
end
t_stim=20;

disp('Saving...');

save(fullfile(fileparts(meta.img_dir),save_file),...
    't','fr','good_ind','K2_total_ROI','K2_shot_ROI','K2_read_ROI',...
    'K2_quantized_ROI','K2_spatial_ROI','K2_fundamental_ROI',...
    'stim_time','t_stim','cc_ROI',...
    'cc_mean','mean_img','mask',...
    'roi_positions','window_size','gain','dark_var', ...
    'mean_img_reference','meta');

%% plot
close all;

roi_ind = 1;

t_sub = t(1:numel(K2_total_ROI(:,roi_ind)));
x_lim = [min(t_sub) max(t_sub)];

[b,a] = butter(3,21/(fr/2),"low"); % 0.01 Hz to 0.5 Hz

var_I = filtfilt(b,a,K2_total_ROI(:,roi_ind)).*(filtfilt(b,a,cc_ROI(:,roi_ind)).^2);

figure('Position',[50 50 500 400]);
subplot(2,1,1); plot(t_sub,sqrt(var_I),'LineWidth',1.5);
xlim(x_lim);
ylabel('std dev'); set(gca,'FontSize',18);
subplot(2,1,2); plot(t_sub,filtfilt(b,a,cc_ROI(:,roi_ind)),'LineWidth',1.5);
xlim(x_lim);
xlabel('time (s)');
ylabel('mean intensity'); set(gca,'FontSize',18);

Iall_smooth = smooth(filtfilt(b,a,cc_ROI(:,roi_ind)),fr);

f1 = figure('Position',[50 50 700 500]);
plot(t_sub,cc_ROI(:,roi_ind)); hold on;
plot(t_sub,Iall_smooth,'LineWidth',2);

ylabel('Intensity (ADU)'); xlabel('Time (s)')
xlim(x_lim);
set(gca,'FontSize',18)

legend('Intensity','Intensity (smoothed)');

f2 = figure('Position',[50 50 700 500]);
plot(t_sub,filtfilt(b,a,K2_total_ROI(:,roi_ind)),'LineWidth',1); hold on;
plot(t_sub,filtfilt(b,a,K2_shot_ROI(:,roi_ind)),'LineWidth',1);
plot(t_sub,filtfilt(b,a,K2_read_ROI(:,roi_ind)),'LineWidth',1);
plot(t_sub,filtfilt(b,a,K2_quantized_ROI(:,roi_ind)),'LineWidth',1);
plot(t_sub,filtfilt(b,a,K2_spatial_ROI(:,roi_ind)),'LineWidth',1);
plot(t_sub,filtfilt(b,a,K2_fundamental_ROI(:,roi_ind)),'LineWidth',1);
ylabel('K^2');
xlabel('time (s)');

legend('K^2 total','K^2 shot','K^2 read','K^2 dig','K^2 spatial','K^2 fund')
xlim(x_lim);
set(gca,'FontSize',18);

f3 = figure('Position',[50 50 700 500]);
subplot(2,1,1);
plot(t_sub,filtfilt(b,a,K2_fundamental_ROI(:,roi_ind)),'LineWidth',1);
ylabel('K_f^2');
xlabel('time (s)');
xlim(x_lim);
set(gca,'FontSize',18);
subplot(2,1,2);
plot(t_sub,filtfilt(b,a,1./K2_fundamental_ROI(:,roi_ind)),'LineWidth',1);
ylabel('BFi');
xlabel('time (s)');
xlim(x_lim);
set(gca,'FontSize',18);