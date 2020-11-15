% Add toolboxes to the path
ftDir = '~/Projects/fieldtrip-20191213/';
if exist('ft_defaults.m', 'file') == 0
    addpath(ftDir); ft_defaults;
end

subList = {'IR9', 'IR10', 'IR16', 'IR26' 'IR28', 'IR35', 'IR57'}

for idx = 1:length(subList)
    clear elec_acpc_f cortex_lh cortex_rh

    % define parameters
    patientCode = subList{idx};
    anatomicalDir = dir(['~/Projects/knight_server/remote/ecog/DATA_FOLDER/Irvine/' patientCode '/3D_Images/Recon_*/']);
    anatomicalPath =  anatomicalDir(1).folder;
    faceAlpha = .3;
    elecSize = 40;
    laterality = 'b'; % change this from 'b' (bilateral/both sides) to 'l' or 'r' to only plot left or right electrodes
    viewPoint = [-90 0]; % left view
    % viewPoint = [90 0]; % right view


    % load surfaces
    load([anatomicalPath '/FT_Pipeline/Surfaces/' patientCode '_cortex_lh.mat'])
    load([anatomicalPath '/FT_Pipeline/Surfaces/' patientCode '_cortex_rh.mat'])

    % load electrodes
    load([anatomicalPath '/FT_Pipeline/Electrodes/' patientCode '_elec_acpc_f.mat'])

    % choose only left or right electrodes
    if ~strcmpi(laterality, 'b')
        idxElecs = find(cellfun(@(x) strcmpi(x(1), laterality), elec_acpc_f.label));
        elec_acpc_f.chanpos = elec_acpc_f.chanpos(idxElecs, :);
        elec_acpc_f.chantype = elec_acpc_f.chantype(idxElecs);
        elec_acpc_f.chanunit = elec_acpc_f.chanunit(idxElecs);
        elec_acpc_f.elecpos = elec_acpc_f.elecpos(idxElecs, :);
        elec_acpc_f.label = elec_acpc_f.label(idxElecs);
        elec_acpc_f.chanside = elec_acpc_f.chanside(idxElecs);
        idxDepths = cellfun(@(x) strcmpi(x{1}(1), laterality), elec_acpc_f.depths);
        elec_acpc_f.depths = elec_acpc_f.depths(idxDepths);
        elec_acpc_f.grid_dim = elec_acpc_f.grid_dim(idxDepths, :);
    end

    % color electrodes
    elecLabel = cellfun(@(x) x(1:3), elec_acpc_f.label, 'un', 0);
    elecLabelList = unique(elecLabel);
    elecLabelId = cellfun(@(x) find(ismember(elecLabelList, x)), elecLabel);
    colPalette = lines(length(elecLabelList));
    elecColors = colPalette(elecLabelId, :);


    % plot surfaces
    hf = figure('Position', [1 1 640 640], 'Color', [1 1 1], 'Name', patientCode);
    switch laterality
        case 'b'
            h1 = ft_plot_mesh(cortex_lh, 'facealpha', faceAlpha);
            ft_plot_mesh(cortex_rh, 'facealpha', faceAlpha);
        case 'l'
            h1 = ft_plot_mesh(cortex_lh, 'facealpha', faceAlpha);
        case 'r'
            h1 = ft_plot_mesh(cortex_rh, 'facealpha', faceAlpha);
    end
    ft_plot_sens(elec_acpc_f, 'elecsize', elecSize, 'facecolor', elecColors, 'label', 'label');

    % turn the light on
    view([-90 0]); material dull; lighting gouraud; camlight;
    view([90 0]); material dull; lighting gouraud; camlight;
    view(viewPoint);

    dcm = datacursormode(hf);
    set(dcm, 'update', {@cursorOutputFunction, elec_acpc_f.label, [], []});

end
