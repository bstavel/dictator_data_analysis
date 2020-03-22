%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script to view the recons for the subjects with the dictator task %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load fieldtrip %%
addpath('~/Projects/fieldtrip-20191213/')
ft_defaults

%% IR9 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR9/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR9_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR9_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)

%% IR10 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR10/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR10_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR10_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)

%% IR14 %%
  %% confusing behavior, and low IQ, so we are excluding

%% IR16 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR16/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR16_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR16_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)

%% IR17 %%
  %% I cannot find any post-op scans, so no recon. But there are new elec notes
  %% from a recent meeting, within last 6 weeks. No insula so I don't think it matters
  %% but a ton of acc and some ofc, ironically.


%% IR19 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR19/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR19_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR19_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)

%% IR26 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR26/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR26_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR26_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)

%% IR28 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR28/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR28_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR28_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)

%% IR35 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR35/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR35_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR35_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)

%% IR39 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR39/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR39_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR39_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)

%% IR57 %%
  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR57/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR57_fsMR_post_acpc.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR57_elec_acpc_f.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)


  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/IR57/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/IR57_fsMR_mni_v.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/IR57_elec_mni_frv.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)





  clearvars
  % path %
  path_widcard = '/Users/bstavel/Projects/knight_server/remote/ecog/DATA_FOLDER/Stanford/ST40/3D_Images/Recon_*_20*/FT_Pipeline';
  dir_search = dir(path_widcard);
  path = dir_search(1).folder;
  % load post-op mri %
  mri = ft_read_mri(sprintf('%s/Scans/ST40_fsMR_mni_v.nii', path))
  % load electrode placement %
  load(sprintf('%s/Electrodes/ST40_elec_mni_v.mat', path))
  % create config %
  cfg = [];
  cfg.elec = elec_acpc_f;
  ft_electrodeplacement(cfg, mri)  
