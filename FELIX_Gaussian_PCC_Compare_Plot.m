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

    %% 4. Pre-Process All Gaussian Data 
    disp('Calculating optimal shifts for all files...');
    batchData = struct(); 
    list_names = cell(1, num_files); % Array to hold names for the reorder listbox
    
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
            else, detected_basis = raw_basis; end
        end
        batchData(i).TheoryStr = sprintf('%s\n%s', detected_func, detected_basis); 
        list_names{i} = sprintf('%s / %s', detected_func, detected_basis); % Single line for the listbox
        
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
        
        % Calculate Baseline & Best PCC
        interp_env_orig = interp1(env_x, batchData(i).env_y_norm, WaveNo, 'linear', 0);
        batchData(i).orig_pcc = corr(IRMPD_norm, interp_env_orig, 'Type', 'Pearson');

        bestCorr = -1; bestShift = 0;
        for shift = -50:0.5:50
            interp_env = interp1(env_x + shift, batchData(i).env_y_norm, WaveNo, 'linear', 0);
            cVal = corr(IRMPD_norm, interp_env, 'Type', 'Pearson');
            if cVal > bestCorr
                bestCorr = cVal;
                bestShift = shift;
            end
        end
        batchData(i).bestShift = bestShift;
        batchData(i).bestCorr = bestCorr;
    end
    disp('Initialization complete! Opening GUI...');

    %% 5. Build the GUI Window
    fig = uifigure('Name', 'IRMPD Batch Comparison GUI', 'Position', [50, 50, 1300, 800]);
    
    % Left Panel: Controls (Fixed in place)
    pnl = uipanel(fig, 'Title', 'Global Controls', 'Position', [10, 10, 240, 780], 'FontSize', 12, 'FontWeight', 'bold');
    
    % Right Panel: Plots Container (Scrollable!)
    plot_pnl = uipanel(fig, 'Position', [260, 10, 1020, 780], 'Scrollable', 'on', 'BackgroundColor', 'w');
    
    % --- Create the MASSIVE Inner Canvas ---
    inner_height = max(760, num_files * 350); 
    
    % The inner panel sits inside the scrollable panel. 
    inner_pnl = uipanel(plot_pnl, 'Position', [0, 0, 990, inner_height], 'BorderType', 'none', 'BackgroundColor', 'w');
    
    % 1 row per file, 5 columns wide (4 for plot, 1 for text)
    t = tiledlayout(inner_pnl, num_files, 5, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    color_palette = lines(num_files); 
    axs_plot = cell(num_files, 1);
    axs_text = cell(num_files, 1);
    
    for i = 1:num_files
        axs_plot{i} = nexttile(t, [1 4]); % Plot takes 4 columns
        axs_text{i} = nexttile(t, [1 1]); % Text gets 1 column on the right
        axis(axs_text{i}, 'off');         % Hide axes borders for the text column
    end

    % --- UI Elements in Left Panel ---
    uibutton(pnl, 'push', 'Text', 'Export Figure to PNG', ...
        'Position', [20, 720, 200, 30], 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.9 0.8], ...
        'ButtonPushedFcn', @(src, event) exportFigure(inner_pnl, irmpd_filename)); 

    uilabel(pnl, 'Text', 'Overall Graph Heading:', 'Position', [20, 690, 200, 22], 'FontWeight', 'bold');
    titleEdit = uieditfield(pnl, 'text', 'Position', [20, 665, 200, 22], ...
        'Value', 'IRMPD vs Calculated Spectra', 'ValueChangedFcn', @(src, event) updatePlots());

    uilabel(pnl, 'Text', 'Global X-Axis Limits:', 'Position', [20, 630, 200, 22], 'FontWeight', 'bold');
    uilabel(pnl, 'Text', 'Min:', 'Position', [20, 600, 40, 22]);
    xMinEdit = uieditfield(pnl, 'numeric', 'Position', [60, 600, 60, 22], ...
        'Value', 800, 'ValueChangedFcn', @(src, event) updatePlots());
    uilabel(pnl, 'Text', 'Max:', 'Position', [130, 600, 40, 22]);
    xMaxEdit = uieditfield(pnl, 'numeric', 'Position', [170, 600, 60, 22], ...
        'Value', 2000, 'ValueChangedFcn', @(src, event) updatePlots());

    uilabel(pnl, 'Text', 'Global Y-Axis Limits:', 'Position', [20, 560, 200, 22], 'FontWeight', 'bold');
    uilabel(pnl, 'Text', 'Min:', 'Position', [20, 530, 40, 22]);
    yMinEdit = uieditfield(pnl, 'numeric', 'Position', [60, 530, 60, 22], ...
        'Value', 0, 'ValueChangedFcn', @(src, event) updatePlots());
    uilabel(pnl, 'Text', 'Max:', 'Position', [130, 530, 40, 22]);
    yMaxEdit = uieditfield(pnl, 'numeric', 'Position', [170, 530, 60, 22], ...
        'Value', 1.5, 'ValueChangedFcn', @(src, event) updatePlots());

    uilabel(pnl, 'Text', 'IRMPD Y-Offset:', 'Position', [20, 480, 120, 22], 'FontWeight', 'bold');
    offsetEdit = uieditfield(pnl, 'numeric', 'Position', [150, 480, 70, 22], ...
        'Value', 0, 'ValueChangedFcn', @(src, event) updatePlots());

    uilabel(pnl, 'Text', 'Calc Y-Scale:', 'Position', [20, 440, 120, 22], 'FontWeight', 'bold');
    yscaleEdit = uieditfield(pnl, 'numeric', 'Position', [150, 440, 70, 22], ...
        'Value', 1, 'ValueChangedFcn', @(src, event) updatePlots());

    % NEW REORDER UI ELEMENTS
    uilabel(pnl, 'Text', 'Plot Order (Select & Move):', 'Position', [20, 390, 200, 22], 'FontWeight', 'bold');
    orderList = uilistbox(pnl, 'Position', [20, 150, 200, 235], 'Items', list_names, 'ItemsData', 1:num_files);
    
    uibutton(pnl, 'push', 'Text', 'Move Up ↑', 'Position', [20, 110, 95, 30], ...
        'ButtonPushedFcn', @(src, event) moveItem(-1));
    uibutton(pnl, 'push', 'Text', 'Move Down ↓', 'Position', [125, 110, 95, 30], ...
        'ButtonPushedFcn', @(src, event) moveItem(1));

    % Initial Plot Draw
    updatePlots();
    
    % Scroll to top of the panel automatically
    try scroll(plot_pnl, 'top'); catch; end 

    %% 6. Helper Function for Reordering the List
    function moveItem(direction)
        val = orderList.Value; % Get the currently selected original index
        idx = find([orderList.ItemsData] == val); % Find its current position in the list
        
        if direction == -1 && idx > 1
            % Swap with the item above it
            tempItem = orderList.Items{idx-1};
            tempData = orderList.ItemsData(idx-1);
            
            orderList.Items{idx-1} = orderList.Items{idx};
            orderList.ItemsData(idx-1) = orderList.ItemsData(idx);
            
            orderList.Items{idx} = tempItem;
            orderList.ItemsData(idx) = tempData;
            
            orderList.Value = val; % Keep the same item highlighted
            updatePlots(); % Redraw the plots!
            
        elseif direction == 1 && idx < num_files
            % Swap with the item below it
            tempItem = orderList.Items{idx+1};
            tempData = orderList.ItemsData(idx+1);
            
            orderList.Items{idx+1} = orderList.Items{idx};
            orderList.ItemsData(idx+1) = orderList.ItemsData(idx);
            
            orderList.Items{idx} = tempItem;
            orderList.ItemsData(idx) = tempData;
            
            orderList.Value = val; % Keep the same item highlighted
            updatePlots(); % Redraw the plots!
        end
    end

    %% 7. Update Plots Function
    function updatePlots()
        x_lims = [xMinEdit.Value, xMaxEdit.Value];
        y_lims = [yMinEdit.Value, yMaxEdit.Value];
        y_offset = offsetEdit.Value;
        calc_yscale = yscaleEdit.Value;
        
        % Read the current visual order from the ListBox!
        plotOrder = orderList.ItemsData; 
        
        for i = 1:num_files
            % --- 1. Draw the Stacked Graph ---
            ax = axs_plot{i};
            cla(ax); hold(ax, 'on');
            
            idx = plotOrder(i); % Look up the original data index based on our new order
            
            shift = batchData(idx).bestShift;
            shifted_env_x = batchData(idx).env_x + shift;
            shifted_sticks_x = batchData(idx).sticks_x + shift;
            
            % Use the original color anchored to the data so colors don't jump around when you reorder them
            calcColor = color_palette(idx, :); 
            plot(ax, shifted_env_x, batchData(idx).env_y_norm * calc_yscale, 'Color', calcColor, 'LineWidth', 1.5);
            stem(ax, shifted_sticks_x, batchData(idx).sticks_y_norm * calc_yscale, 'Color', calcColor, 'Marker', 'none', 'LineWidth', 1.5); 
            plot(ax, WaveNo, IRMPD_norm + y_offset, 'k-', 'LineWidth', 1.5);
            
            grid(ax, 'on');
            xlim(ax, x_lims);
            ylim(ax, y_lims);
            ax.YTickLabel = []; % Hide Y-axis numbers to keep it clean
            ax.FontSize = 11;
            
            % Add X-axis label to EVERY plot
            xlabel(ax, '\bf\it{Wavenumber (cm^{-1})}', 'FontSize', 14);
            hold(ax, 'off');
            
            % --- 2. Draw the Text Box on the Right ---
            ax_txt = axs_text{i};
            cla(ax_txt); axis(ax_txt, 'off');
            
            plot_text = sprintf('\\bf%s\n\\rmShift: %+.1f cm^{-1}\nOriginal PCC: %.3f\nOptimized PCC: %.3f', ...
                batchData(idx).TheoryStr, shift, batchData(idx).orig_pcc, batchData(idx).bestCorr);
            
            % Draw large, readable text aligned to the left of the 5th column
            text(ax_txt, 0, 0.5, plot_text, 'Units', 'normalized', ...
                'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'none', ...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
        end
        
        % Set the Overall Graph Heading
        title(t, titleEdit.Value, 'FontSize', 18, 'FontWeight', 'bold');
    end

    %% 8. Export Function
    function exportFigure(canvas_handle, base_filename)
        [~, name, ~] = fileparts(base_filename);
        default_save = sprintf('%s_Batch_Compare.png', name);
        [file, path] = uiputfile('*.png', 'Save Stacked Plot as PNG', default_save);
        
        if ischar(file)
            % This exports the MASSIVE scrollable inner canvas to a single high-res PNG file
            exportgraphics(canvas_handle, fullfile(path, file), 'Resolution', 300);
            disp(['Figure successfully saved to: ', fullfile(path, file)]);
        end
    end
end
