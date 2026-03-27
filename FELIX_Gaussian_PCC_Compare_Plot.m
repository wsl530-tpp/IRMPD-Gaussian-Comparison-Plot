function FELIX_Gaussian_PCC_Batch_GUI()
    %% Clean Workspace
    close all; clc;

    %% 1. Select IRMPD File
    disp('Please select the IRMPD .csv data file...');
    [irmpd_filename, irmpd_pathname] = uigetfile('*.csv', 'Select the IRMPD CSV Data File');
    if isequal(irmpd_filename, 0)
        disp('User canceled IRMPD file selection. Exiting...');
        return;
    end
    irmpd_file = fullfile(irmpd_pathname, irmpd_filename);

    %% 2. Select MULTIPLE Gaussian Files (Loop for multiple folders)
    disp('Please select Gaussian .txt files...');
    gauss_fullpaths = {}; 
    gauss_filenames = {}; 
    
    while true
        [g_names, g_path] = uigetfile('*.txt', 'Select Gaussian TXT Files', 'MultiSelect', 'on');
        
        if isequal(g_names, 0)
            if isempty(gauss_fullpaths)
                disp('No Gaussian files selected. Exiting...');
                return;
            else
                break;
            end
        end
        
        if ischar(g_names), g_names = {g_names}; end
        
        for k = 1:length(g_names)
            gauss_fullpaths{end+1} = fullfile(g_path, g_names{k});
            gauss_filenames{end+1} = g_names{k};
        end
        
        choice = questdlg('Do you want to select more files from another folder?', ...
            'Select More Files?', 'Yes', 'No (I am done)', 'No (I am done)');
        if strcmp(choice, 'No (I am done)') || isempty(choice)
            break;
        end
    end
    num_files = length(gauss_fullpaths);

    %% 3. Load IRMPD Data
    opts = detectImportOptions(irmpd_file);
    irmpd_table = readtable(irmpd_file, opts);

    WaveNo = irmpd_table{:, 1}; 
    IRMPD = irmpd_table{:, 3};  

    valid_idx = ~isnan(WaveNo) & ~isnan(IRMPD);
    WaveNo = WaveNo(valid_idx);
    IRMPD = IRMPD(valid_idx);
    IRMPD_norm = IRMPD / max(abs(IRMPD));

    %% 4. Pre-Process All Gaussian Data (Load raw data only)
    disp('Loading calculation files...');
    batchData = struct(); 
    list_names = cell(1, num_files); 
    
    for i = 1:num_files
        % Parse Filename
        [~, name_no_ext, ~] = fileparts(gauss_filenames{i});
        name_parts = split(name_no_ext, '-');
        
        detected_func = 'Unknown Functional'; detected_basis = 'Unknown Basis';
        if length(name_parts) >= 3
            raw_func = lower(name_parts{2});
            if contains(raw_func, 'camb3lyp') || contains(raw_func, 'cam-b3lyp'), detected_func = 'CAM-B3LYP';
            elseif contains(raw_func, 'b3lyp'), detected_func = 'B3LYP';
            elseif contains(raw_func, 'm062x') || contains(raw_func, 'm06-2x'), detected_func = 'M06-2X';
            elseif contains(raw_func, 'm06'), detected_func = 'M06';
            elseif contains(raw_func, 'pbe0'), detected_func = 'PBE0';
            elseif contains(raw_func, 'rohf'), detected_func = 'ROHF';
            elseif contains(raw_func, 'hf'), detected_func = 'HF';
            else, detected_func = upper(raw_func); end
            
            raw_basis = lower(name_parts{3});
            if contains(raw_basis, '63113df3dp') || contains(raw_basis, '6311++g3df3dp'), detected_basis = '6-311++G(3df,3dp)';
            elseif contains(raw_basis, '63112d2p') || contains(raw_basis, '6311++g2d2p'), detected_basis = '6-311++G(2d,2p)';
            elseif contains(raw_basis, '6311dp') || contains(raw_basis, '6311++gdp') || contains(raw_basis, '6311gdp'), detected_basis = '6-311++G(d,p)';
            elseif contains(raw_basis, 'def2tzvp') || contains(raw_basis, 'def2-tzvp'), detected_basis = 'def2-TZVP';
            elseif contains(raw_basis, 'def2svp') || contains(raw_basis, 'def2-svp'), detected_basis = 'def2-SVP';
            elseif contains(raw_basis, 'augccpvdz'), detected_basis = 'aug-cc-pVDZ';
            elseif contains(raw_basis, 'augccpvtz'), detected_basis = 'aug-cc-pVTZ';
            elseif contains(raw_basis, 'ccpvdz'), detected_basis = 'cc-pVDZ';
            elseif contains(raw_basis, 'cep31g') || contains(raw_basis, 'cep-31g'), detected_basis = 'CEP-31G';
            else, detected_basis = raw_basis; end
        end
        batchData(i).TheoryStr = sprintf('%s\n%s', detected_func, detected_basis); 
        list_names{i} = sprintf('%s / %s', detected_func, detected_basis); 
        
        % DYNAMICALLY FETCH NIST CCCBDB SCALING FACTOR
        batchData(i).freq_scale = getNISTScalingFactor(detected_func, detected_basis);
        
        % Load File
        fileID = fopen(gauss_fullpaths{i}, 'r');
        sticks_x = []; sticks_y = []; env_x = []; env_y = [];
        while ~feof(fileID)
            line = strtrim(fgetl(fileID));
            if startsWith(line, '#')
                parts = sscanf(line, '# %f %f');
                if length(parts) == 2, sticks_x(end+1) = parts(1); sticks_y(end+1) = parts(2); end
            elseif ~isempty(line)
                parts = sscanf(line, '%f %f %f');
                if length(parts) >= 2, env_x(end+1) = parts(1); env_y(end+1) = parts(2); end
            end
        end
        fclose(fileID);

        max_env = max(env_y); if max_env == 0, max_env = 1; end 
        batchData(i).env_x = env_x;
        batchData(i).env_y_norm = env_y / max_env;
        batchData(i).sticks_x = sticks_x;
        batchData(i).sticks_y_norm = sticks_y / max_env;
    end
    disp('Initialization complete! Opening GUI...');

    %% 5. Build the GUI Window
    fig = uifigure('Name', 'IRMPD Batch Comparison GUI', 'Position', [50, 50, 1400, 850]);
    
    % --- SCROLLABLE LEFT PANEL FOR CONTROLS ---
    inner_height_left = max(800, 950 + 35 * num_files);
    left_scroll_pnl = uipanel(fig, 'Position', [10, 10, 260, 830], 'Scrollable', 'on');
    pnl = uipanel(left_scroll_pnl, 'Position', [0, max(0, 830 - inner_height_left), 240, inner_height_left], 'BorderType', 'none');
    
    % --- SCROLLABLE RIGHT PANEL FOR PLOTS ---
    plot_pnl = uipanel(fig, 'Position', [280, 10, 1100, 830], 'Scrollable', 'on', 'BackgroundColor', 'w');
    
    % Default canvas values
    default_canvas_width = 1200;
    default_plot_height = 350;
    
    inner_pnl = uipanel(plot_pnl, 'Position', [0, 0, default_canvas_width, num_files * default_plot_height], 'BorderType', 'none', 'BackgroundColor', 'w');
    
    % Changed to 8 columns: 5 for plot, 3 for text. Gives MASSIVE room to text to prevent cutoff.
    t = tiledlayout(inner_pnl, num_files, 8, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    color_palette = lines(num_files); 
    axs_plot = cell(num_files, 1);
    axs_text = cell(num_files, 1);
    
    for i = 1:num_files
        axs_plot{i} = nexttile(t, [1 5]); 
        axs_text{i} = nexttile(t, [1 3]); % Text gets 3 columns now!
        axis(axs_text{i}, 'off');         
    end

    % --- UI Elements Placement (Top to Bottom) ---
    y_pos = inner_height_left - 40;
    
    uibutton(pnl, 'push', 'Text', 'Export Figure to PNG', ...
        'Position', [20, y_pos, 200, 30], 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.9 0.8], ...
        'ButtonPushedFcn', @(src, event) exportFigure(inner_pnl, irmpd_filename)); 

    y_pos = y_pos - 40;
    uilabel(pnl, 'Text', 'Canvas Dimensions (px):', 'Position', [20, y_pos, 200, 22], 'FontWeight', 'bold');
    y_pos = y_pos - 25;
    uilabel(pnl, 'Text', 'Width:', 'Position', [20, y_pos, 40, 22]);
    canvasWidthEdit = uieditfield(pnl, 'numeric', 'Position', [60, y_pos, 60, 22], ...
        'Value', default_canvas_width, 'ValueChangedFcn', @(src, event) resizeCanvas());
    uilabel(pnl, 'Text', 'Height:', 'Position', [130, y_pos, 45, 22]);
    plotHeightEdit = uieditfield(pnl, 'numeric', 'Position', [175, y_pos, 55, 22], ...
        'Value', default_plot_height, 'ValueChangedFcn', @(src, event) resizeCanvas());

    y_pos = y_pos - 40;
    uilabel(pnl, 'Text', 'Font Setup:', 'Position', [20, y_pos, 200, 22], 'FontWeight', 'bold');
    y_pos = y_pos - 25;
    fontDrop = uidropdown(pnl, 'Items', {'Helvetica', 'Arial', 'Times New Roman', 'Courier New'}, ...
        'Position', [20, y_pos, 130, 22], 'Value', 'Helvetica', 'ValueChangedFcn', @(src, event) updatePlots());
    uilabel(pnl, 'Text', 'Size:', 'Position', [160, y_pos, 30, 22]);
    fontSizeEdit = uieditfield(pnl, 'numeric', 'Position', [195, y_pos, 35, 22], ...
        'Value', 12, 'ValueChangedFcn', @(src, event) updatePlots());

    y_pos = y_pos - 40;
    uilabel(pnl, 'Text', 'Overall Graph Heading:', 'Position', [20, y_pos, 200, 22], 'FontWeight', 'bold');
    y_pos = y_pos - 25;
    titleEdit = uieditfield(pnl, 'text', 'Position', [20, y_pos, 200, 22], ...
        'Value', 'IRMPD vs Calculated Spectra', 'ValueChangedFcn', @(src, event) updatePlots());

    y_pos = y_pos - 40;
    uilabel(pnl, 'Text', 'Global X-Axis Limits:', 'Position', [20, y_pos, 200, 22], 'FontWeight', 'bold');
    y_pos = y_pos - 25;
    uilabel(pnl, 'Text', 'Min:', 'Position', [20, y_pos, 40, 22]);
    xMinEdit = uieditfield(pnl, 'numeric', 'Position', [60, y_pos, 60, 22], ...
        'Value', 800, 'ValueChangedFcn', @(src, event) updatePlots());
    uilabel(pnl, 'Text', 'Max:', 'Position', [130, y_pos, 40, 22]);
    xMaxEdit = uieditfield(pnl, 'numeric', 'Position', [170, y_pos, 60, 22], ...
        'Value', 2000, 'ValueChangedFcn', @(src, event) updatePlots());

    y_pos = y_pos - 40;
    uilabel(pnl, 'Text', 'Global Y-Axis Limits:', 'Position', [20, y_pos, 200, 22], 'FontWeight', 'bold');
    y_pos = y_pos - 25;
    uilabel(pnl, 'Text', 'Min:', 'Position', [20, y_pos, 40, 22]);
    yMinEdit = uieditfield(pnl, 'numeric', 'Position', [60, y_pos, 60, 22], ...
        'Value', 0, 'ValueChangedFcn', @(src, event) updatePlots());
    uilabel(pnl, 'Text', 'Max:', 'Position', [130, y_pos, 40, 22]);
    yMaxEdit = uieditfield(pnl, 'numeric', 'Position', [170, y_pos, 60, 22], ...
        'Value', 1.5, 'ValueChangedFcn', @(src, event) updatePlots());

    y_pos = y_pos - 40;
    uilabel(pnl, 'Text', 'IRMPD Y-Offset:', 'Position', [20, y_pos, 120, 22], 'FontWeight', 'bold');
    offsetEdit = uieditfield(pnl, 'numeric', 'Position', [150, y_pos, 70, 22], ...
        'Value', 0, 'ValueChangedFcn', @(src, event) updatePlots());

    y_pos = y_pos - 35;
    uilabel(pnl, 'Text', 'Calc Y-Scale (Int):', 'Position', [20, y_pos, 120, 22], 'FontWeight', 'bold');
    yscaleEdit = uieditfield(pnl, 'numeric', 'Position', [150, y_pos, 70, 22], ...
        'Value', 1, 'ValueChangedFcn', @(src, event) updatePlots());

    y_pos = y_pos - 35;
    autoPCCBox = uicheckbox(pnl, 'Text', ' Auto-Optimize Shift (PCC)', ...
        'Position', [20, y_pos, 200, 22], 'Value', 0, 'FontWeight', 'bold', ...
        'ValueChangedFcn', @(src, event) updatePlots());

    y_pos = y_pos - 40;
    uilabel(pnl, 'Text', 'Plot Order (Select & Move):', 'Position', [20, y_pos, 200, 22], 'FontWeight', 'bold');
    y_pos = y_pos - 150;
    orderList = uilistbox(pnl, 'Position', [20, y_pos, 200, 150], 'Items', list_names, 'ItemsData', 1:num_files);
    
    y_pos = y_pos - 35;
    uibutton(pnl, 'push', 'Text', 'Move Up ↑', 'Position', [20, y_pos, 95, 30], ...
        'ButtonPushedFcn', @(src, event) moveItem(-1));
    uibutton(pnl, 'push', 'Text', 'Move Down ↓', 'Position', [125, y_pos, 95, 30], ...
        'ButtonPushedFcn', @(src, event) moveItem(1));

    % --- DYNAMIC SCALING BOXES FOR EACH FUNCTIONAL ---
    y_pos = y_pos - 40;
    uilabel(pnl, 'Text', 'X-Scale Factors (Freq):', 'Position', [20, y_pos, 200, 22], 'FontWeight', 'bold');
    
    scaleEdits = gobjects(num_files, 1);
    for i = 1:num_files
        y_pos = y_pos - 30;
        
        short_name = list_names{i};
        if length(short_name) > 18
            short_name = [short_name(1:16) '..'];
        end
        
        uilabel(pnl, 'Text', short_name, 'Position', [20, y_pos, 130, 22]);
        
        % The value is pre-filled with the NIST Dictionary lookup
        scaleEdits(i) = uieditfield(pnl, 'numeric', 'Position', [150, y_pos, 70, 22], ...
            'Value', batchData(i).freq_scale, 'ValueChangedFcn', @(src, event) updatePlots());
    end

    % Initial Plot Draw
    updatePlots();
    
    % Scroll both panels to the top automatically
    try scroll(plot_pnl, 'top'); catch; end 
    try scroll(left_scroll_pnl, 'top'); catch; end

    %% 6. Helper Functions for Reordering and Resizing
    function resizeCanvas()
        % Dynamically stretch the drawing canvas without changing the GUI window
        new_w = canvasWidthEdit.Value;
        new_h = plotHeightEdit.Value * num_files;
        inner_pnl.Position = [0, 0, new_w, new_h];
        updatePlots();
    end

    function moveItem(direction)
        val = orderList.Value; 
        idx = find([orderList.ItemsData] == val); 
        
        if direction == -1 && idx > 1
            tempItem = orderList.Items{idx-1};
            tempData = orderList.ItemsData(idx-1);
            orderList.Items{idx-1} = orderList.Items{idx};
            orderList.ItemsData(idx-1) = orderList.ItemsData(idx);
            orderList.Items{idx} = tempItem;
            orderList.ItemsData(idx) = tempData;
            orderList.Value = val; 
            updatePlots(); 
            
        elseif direction == 1 && idx < num_files
            tempItem = orderList.Items{idx+1};
            tempData = orderList.ItemsData(idx+1);
            orderList.Items{idx+1} = orderList.Items{idx};
            orderList.ItemsData(idx+1) = orderList.ItemsData(idx);
            orderList.Items{idx} = tempItem;
            orderList.ItemsData(idx) = tempData;
            orderList.Value = val; 
            updatePlots(); 
        end
    end

    %% 7. Update Plots Function
    function updatePlots()
        x_lims = [xMinEdit.Value, xMaxEdit.Value];
        y_lims = [yMinEdit.Value, yMaxEdit.Value];
        y_offset = offsetEdit.Value;
        calc_yscale = yscaleEdit.Value;
        is_auto_pcc = autoPCCBox.Value; 
        
        fName = fontDrop.Value;
        fSize = fontSizeEdit.Value;
        
        plotOrder = orderList.ItemsData; 
        
        for i = 1:num_files
            ax = axs_plot{i};
            cla(ax); hold(ax, 'on');
            
            idx = plotOrder(i); 
            freq_scale = scaleEdits(idx).Value; 
            
            scaled_env_x = batchData(idx).env_x * freq_scale;
            scaled_sticks_x = batchData(idx).sticks_x * freq_scale;
            
            interp_env_orig = interp1(scaled_env_x, batchData(idx).env_y_norm, WaveNo, 'linear', 0);
            orig_pcc = corr(IRMPD_norm, interp_env_orig, 'Type', 'Pearson');

            bestCorr = orig_pcc; 
            bestShift = 0;
            
            if is_auto_pcc
                for shift = -50:0.5:50
                    interp_env = interp1(scaled_env_x + shift, batchData(idx).env_y_norm, WaveNo, 'linear', 0);
                    cVal = corr(IRMPD_norm, interp_env, 'Type', 'Pearson');
                    if cVal > bestCorr
                        bestCorr = cVal;
                        bestShift = shift;
                    end
                end
            end
            
            shifted_env_x = scaled_env_x + bestShift;
            shifted_sticks_x = scaled_sticks_x + bestShift;
            
            calcColor = color_palette(idx, :); 
            plot(ax, shifted_env_x, batchData(idx).env_y_norm * calc_yscale, 'Color', calcColor, 'LineWidth', 1.5);
            stem(ax, shifted_sticks_x, batchData(idx).sticks_y_norm * calc_yscale, 'Color', calcColor, 'Marker', 'none', 'LineWidth', 1.5); 
            plot(ax, WaveNo, IRMPD_norm + y_offset, 'k-', 'LineWidth', 1.5);
            
            grid(ax, 'on');
            xlim(ax, x_lims);
            ylim(ax, y_lims);
            ax.YTickLabel = []; 
            
            % Apply chosen Font Name and Size to axes
            ax.FontName = fName;
            ax.FontSize = fSize - 1; 
            
            xlabel(ax, '\bf\it{Wavenumber (cm^{-1})}', 'FontName', fName, 'FontSize', fSize + 2);
            hold(ax, 'off');
            
            % --- Draw the Formatted Text Box on the Right ---
            ax_txt = axs_text{i};
            cla(ax_txt); axis(ax_txt, 'off');
            
            if is_auto_pcc
                plot_text = sprintf('\\bfFunctional / Basis Set:\\rm %s\n\\bfScaling factor:\\rm %.4f\n\\bfOriginal PCC:\\rm %.3f\n\\bfOptimized PCC:\\rm %.3f\n\\bfShift:\\rm %+.1f cm^{-1}', ...
                    list_names{idx}, freq_scale, orig_pcc, bestCorr, bestShift);
            else
                plot_text = sprintf('\\bfFunctional / Basis Set:\\rm %s\n\\bfScaling factor:\\rm %.4f\n\\bfOriginal PCC:\\rm %.3f\n\\bfOptimized PCC:\\rm --\n\\bfShift:\\rm OFF', ...
                    list_names{idx}, freq_scale, orig_pcc);
            end
            
            % Apply chosen Font Name and Size to text boxes
            text(ax_txt, 0, 0.5, plot_text, 'Units', 'normalized', ...
                'FontName', fName, 'FontSize', fSize, ...
                'BackgroundColor', 'w', 'EdgeColor', 'none', ...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
        end
        
        % Set the Overall Graph Heading Font
        title(t, titleEdit.Value, 'FontName', fName, 'FontSize', fSize + 6, 'FontWeight', 'bold');
    end

    %% 8. Export Function
    function exportFigure(canvas_handle, base_filename)
        [~, name, ~] = fileparts(base_filename);
        default_save = sprintf('%s_Batch_Compare.png', name);
        [file, path] = uiputfile('*.png', 'Save Stacked Plot as PNG', default_save);
        
        if ischar(file)
            exportgraphics(canvas_handle, fullfile(path, file), 'Resolution', 300);
            disp(['Figure successfully saved to: ', fullfile(path, file)]);
        end
    end

    %% 9. NIST CCCBDB Dictionary Lookup
    function sf = getNISTScalingFactor(func, basis)
        f_clean = upper(regexprep(func, '[-\s]', ''));
        b_clean = upper(regexprep(basis, '[-\s(),+]', ''));
        
        key = [f_clean, '_', b_clean];
        sf_map = containers.Map();
        
        % User Custom Request
        sf_map('ROHF_CEP31G') = 0.9085;
        
        % Common B3LYP Combinations
        sf_map('B3LYP_6311G2D2P') = 0.9688; 
        sf_map('B3LYP_6311G3DF3DP') = 0.9688;
        sf_map('B3LYP_6311GDP') = 0.9679;
        sf_map('B3LYP_AUGCCPVDZ') = 0.9700;
        sf_map('B3LYP_AUGCCPVTZ') = 0.9689;
        sf_map('B3LYP_CCPVDZ') = 0.9700;
        sf_map('B3LYP_DEF2TZVP') = 0.9680;
        
        % Common M06-2X Combinations
        sf_map('M062X_6311G2D2P') = 0.9460; 
        sf_map('M062X_6311GDP') = 0.9460;
        sf_map('M062X_AUGCCPVDZ') = 0.9480;
        sf_map('M062X_DEF2TZVP') = 0.9470;
        
        % Common CAM-B3LYP Combinations
        sf_map('CAMB3LYP_6311GDP') = 0.9530;
        sf_map('CAMB3LYP_AUGCCPVDZ') = 0.9540;
        
        % Common PBE0 Combinations
        sf_map('PBE0_6311GDP') = 0.9550;
        sf_map('PBE0_AUGCCPVDZ') = 0.9560;

        if isKey(sf_map, key)
            sf = sf_map(key);
        else
            switch f_clean
                case 'B3LYP'
                    sf = 0.9680;
                case 'M062X'
                    sf = 0.9470;
                case 'CAMB3LYP'
                    sf = 0.9530;
                case 'PBE0'
                    sf = 0.9550;
                case 'ROHF'
                    sf = 0.9085;
                case 'HF'
                    sf = 0.8980;
                otherwise
                    sf = 1.0000;
            end
        end
    end
end
