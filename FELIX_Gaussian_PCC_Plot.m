function IRMPD_GUI_Viewer_SinglePlot()
    %% Clean Workspace
    close all; clc;

    %% 1. Select Input Files (Initial Load)
    disp('Please select the IRMPD .csv data file...');
    last_csv_path = getpref('FELIX_GUI_Single', 'last_csv_path', pwd);
    [irmpd_filename, irmpd_pathname] = uigetfile(fullfile(last_csv_path, '*.csv'), 'Select the IRMPD CSV Data File');
    
    if isequal(irmpd_filename, 0)
        disp('User canceled IRMPD file selection. Exiting...');
        return;
    end
    setpref('FELIX_GUI_Single', 'last_csv_path', irmpd_pathname);
    irmpd_file = fullfile(irmpd_pathname, irmpd_filename);

    disp('Please select the Gaussian .txt calculated spectrum file...');
    last_txt_path = getpref('FELIX_GUI_Single', 'last_txt_path', irmpd_pathname);
    [gauss_filename, gauss_pathname] = uigetfile(fullfile(last_txt_path, '*.txt'), 'Select the Gaussian TXT File');
    
    if isequal(gauss_filename, 0)
        disp('User canceled Gaussian file selection. Exiting...');
        return;
    end
    setpref('FELIX_GUI_Single', 'last_txt_path', gauss_pathname);
    gauss_file = fullfile(gauss_pathname, gauss_filename);

    %% 2. Auto-Detect Functional and Basis Set from Filename
    [~, name_no_ext, ~] = fileparts(gauss_filename);
    name_parts = split(name_no_ext, '-');
    
    detected_func = 'Unknown Functional';
    detected_basis = 'Unknown Basis';
    
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

    %% 3. Load IRMPD Data
    opts = detectImportOptions(irmpd_file);
    irmpd_table = readtable(irmpd_file, opts);

    WaveNo = irmpd_table{:, 1}; 
    IRMPD = irmpd_table{:, 3};  

    valid_idx = ~isnan(WaveNo) & ~isnan(IRMPD);
    WaveNo = WaveNo(valid_idx);
    IRMPD = IRMPD(valid_idx);
    
    IRMPD_norm = IRMPD / max(abs(IRMPD));

    %% 4. Load Gaussian Data
    fileID = fopen(gauss_file, 'r');
    sticks_x = []; sticks_y = [];
    env_x = []; env_y = [];

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
    env_y_norm = env_y / max_env;
    sticks_y_norm = sticks_y / max_env;
    
    %% 4.5 Calculate Original Baseline PCC (Shift = 0)
    interp_env_orig = interp1(env_x, env_y_norm, WaveNo, 'linear', 0);
    orig_pcc = corr(IRMPD_norm, interp_env_orig, 'Type', 'Pearson');

    %% 5. Build the GUI
    % Taller Figure to accommodate new Load buttons
    fig = uifigure('Name', 'IRMPD Single Plot GUI', 'Position', [100, 50, 1000, 820]);
    
    ax = uiaxes(fig, 'Position', [20, 20, 720, 780]);
    pnl = uipanel(fig, 'Title', 'Plot Controls', 'Position', [760, 20, 220, 780], 'FontSize', 12, 'FontWeight', 'bold');

    % --- NEW FILE LOADING BUTTONS ---
    uibutton(pnl, 'push', 'Text', 'Load New Calc (.txt)', ...
        'Position', [10, 725, 200, 25], 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.9 1], ...
        'ButtonPushedFcn', @(src, event) loadNewCalc());
        
    uibutton(pnl, 'push', 'Text', 'Load New IRMPD (.csv)', ...
        'Position', [10, 695, 200, 25], 'FontWeight', 'bold', 'BackgroundColor', [1 0.9 0.8], ...
        'ButtonPushedFcn', @(src, event) loadNewIRMPD());

    % --- Existing UI Elements ---
    uilabel(pnl, 'Text', 'Plot Summary/Title:', 'Position', [10, 660, 200, 22]);
    summaryEdit = uieditfield(pnl, 'text', 'Position', [10, 635, 200, 25], ...
        'Value', name_no_ext, 'ValueChangedFcn', @(src, event) updatePlot());

    uilabel(pnl, 'Text', 'Functional:', 'Position', [10, 600, 200, 22]);
    funcEdit = uieditfield(pnl, 'text', 'Position', [10, 575, 200, 25], ...
        'Value', detected_func, 'ValueChangedFcn', @(src, event) updatePlot());

    uilabel(pnl, 'Text', 'Basis Set:', 'Position', [10, 540, 200, 22]);
    basisEdit = uieditfield(pnl, 'text', 'Position', [10, 515, 200, 25], ...
        'Value', detected_basis, 'ValueChangedFcn', @(src, event) updatePlot());

    uilabel(pnl, 'Text', 'Calc Color:', 'Position', [10, 475, 80, 22]);
    colorDrop = uidropdown(pnl, 'Items', {'Blue', 'Red', 'Green', 'Magenta', 'Cyan', 'Orange', 'Purple'}, ...
        'Position', [90, 475, 120, 22], 'ValueChangedFcn', @(src, event) updatePlot());

    uilabel(pnl, 'Text', 'IRMPD Y-Offset:', 'Position', [10, 435, 120, 22]);
    offsetEdit = uieditfield(pnl, 'numeric', 'Position', [130, 435, 80, 22], ...
        'Value', 0, 'ValueChangedFcn', @(src, event) updatePlot());

    uilabel(pnl, 'Text', 'Calc Y-Scale:', 'Position', [10, 395, 120, 22]);
    yscaleEdit = uieditfield(pnl, 'numeric', 'Position', [130, 395, 80, 22], ...
        'Value', 0.3, 'ValueChangedFcn', @(src, event) updatePlot());

    uilabel(pnl, 'Text', 'Linear Shift (cm-1):', 'Position', [10, 355, 120, 22]);
    shiftEdit = uieditfield(pnl, 'numeric', 'Position', [130, 355, 80, 22], ...
        'Value', 0, 'ValueChangedFcn', @(src, event) updatePlot());

    uibutton(pnl, 'push', 'Text', 'Auto-Optimize (PCC)', ...
        'Position', [10, 315, 200, 30], 'ButtonPushedFcn', @(src, event) optimizeShift());

    pccLabel = uilabel(pnl, 'Text', 'Current PCC: --', 'Position', [10, 275, 200, 22], ...
        'FontWeight', 'bold', 'FontColor', [0 0.5 0]);

    uilabel(pnl, 'Text', 'X-Axis Limits:', 'Position', [10, 230, 200, 22], 'FontWeight', 'bold');
    uilabel(pnl, 'Text', 'Min:', 'Position', [10, 200, 40, 22]);
    xMinEdit = uieditfield(pnl, 'numeric', 'Position', [45, 200, 60, 22], ...
        'Value', 800, 'ValueChangedFcn', @(src, event) updatePlot());
    uilabel(pnl, 'Text', 'Max:', 'Position', [115, 200, 40, 22]);
    xMaxEdit = uieditfield(pnl, 'numeric', 'Position', [150, 200, 60, 22], ...
        'Value', 2000, 'ValueChangedFcn', @(src, event) updatePlot());

    uilabel(pnl, 'Text', 'Y-Axis Limits:', 'Position', [10, 150, 200, 22], 'FontWeight', 'bold');
    uilabel(pnl, 'Text', 'Min:', 'Position', [10, 120, 40, 22]);
    yMinEdit = uieditfield(pnl, 'numeric', 'Position', [45, 120, 60, 22], ...
        'Value', 0, 'ValueChangedFcn', @(src, event) updatePlot());
    uilabel(pnl, 'Text', 'Max:', 'Position', [115, 120, 40, 22]);
    yMaxEdit = uieditfield(pnl, 'numeric', 'Position', [150, 120, 60, 22], ...
        'Value', 1.5, 'ValueChangedFcn', @(src, event) updatePlot());

    % Initialize the plot
    updatePlot();

    %% 6. Nested Functions for GUI Actions
    
    % --- Function to Load a NEW Calculation on the fly ---
    function loadNewCalc()
        last_path = getpref('FELIX_GUI_Single', 'last_txt_path', pwd);
        [new_name, new_path] = uigetfile(fullfile(last_path, '*.txt'), 'Select NEW Gaussian TXT File');
        if isequal(new_name, 0), return; end % User canceled
        
        setpref('FELIX_GUI_Single', 'last_txt_path', new_path);
        new_file = fullfile(new_path, new_name);
        
        % Auto-detect new theory
        [~, n_no_ext, ~] = fileparts(new_name);
        n_parts = split(n_no_ext, '-');
        
        d_func = 'Unknown'; d_basis = 'Unknown';
        if length(n_parts) >= 3
            d_func = upper(n_parts{2}); % Simplified fallback
            d_basis = n_parts{3};
            
            % Common cleanups
            if contains(lower(d_func), 'camb3lyp'), d_func = 'CAM-B3LYP'; end
            if contains(lower(d_basis), '63113df3dp'), d_basis = '6-311++G(3df,3dp)'; end
            if contains(lower(d_basis), '63112d2p'), d_basis = '6-311++G(2d,2p)'; end
            if contains(lower(d_basis), '6311dp'), d_basis = '6-311++G(d,p)'; end
            if contains(lower(d_basis), 'def2tzvp'), d_basis = 'def2-TZVP'; end
            if contains(lower(d_basis), 'augccpvdz'), d_basis = 'aug-cc-pVDZ'; end
        end
        
        % Update GUI text boxes
        funcEdit.Value = d_func;
        basisEdit.Value = d_basis;
        summaryEdit.Value = n_no_ext;
        
        % Read new file
        fileID_new = fopen(new_file, 'r');
        sticks_x = []; sticks_y = []; env_x = []; env_y = [];
        while ~feof(fileID_new)
            line = strtrim(fgetl(fileID_new));
            if startsWith(line, '#')
                parts = sscanf(line, '# %f %f');
                if length(parts) == 2, sticks_x(end+1) = parts(1); sticks_y(end+1) = parts(2); end
            elseif ~isempty(line)
                parts = sscanf(line, '%f %f %f');
                if length(parts) >= 2, env_x(end+1) = parts(1); env_y(end+1) = parts(2); end
            end
        end
        fclose(fileID_new);

        % Normalize new data
        m_env = max(env_y); if m_env == 0, m_env = 1; end 
        env_y_norm = env_y / m_env;
        sticks_y_norm = sticks_y / m_env;
        
        % Recalculate original baseline PCC
        interp_env_orig = interp1(env_x, env_y_norm, WaveNo, 'linear', 0);
        orig_pcc = corr(IRMPD_norm, interp_env_orig, 'Type', 'Pearson');
        
        % Reset Shift to 0 and Redraw
        shiftEdit.Value = 0;
        updatePlot();
    end

    % --- Function to Load a NEW IRMPD on the fly ---
    function loadNewIRMPD()
        last_path = getpref('FELIX_GUI_Single', 'last_csv_path', pwd);
        [new_name, new_path] = uigetfile(fullfile(last_path, '*.csv'), 'Select NEW IRMPD CSV Data File');
        if isequal(new_name, 0), return; end % User canceled
        
        setpref('FELIX_GUI_Single', 'last_csv_path', new_path);
        new_file = fullfile(new_path, new_name);
        
        % Read new CSV
        opts_new = detectImportOptions(new_file);
        irmpd_table_new = readtable(new_file, opts_new);

        WaveNo = irmpd_table_new{:, 1}; 
        IRMPD = irmpd_table_new{:, 3};  

        v_idx = ~isnan(WaveNo) & ~isnan(IRMPD);
        WaveNo = WaveNo(v_idx);
        IRMPD_norm = IRMPD(v_idx) / max(abs(IRMPD(v_idx)));
        
        % Recalculate original baseline PCC against the new experimental data
        interp_env_orig = interp1(env_x, env_y_norm, WaveNo, 'linear', 0);
        orig_pcc = corr(IRMPD_norm, interp_env_orig, 'Type', 'Pearson');
        
        % Redraw
        updatePlot();
    end

    % Function to update the plot when values change
    function updatePlot()
        current_shift = shiftEdit.Value;
        y_offset = offsetEdit.Value;
        calc_yscale = yscaleEdit.Value;
        
        summary_text = summaryEdit.Value;
        func_text = funcEdit.Value;
        basis_text = basisEdit.Value;
        
        selectedColorStr = colorDrop.Value;
        colorMap = struct('Blue', 'b', 'Red', 'r', 'Green', 'g', 'Magenta', 'm', ...
                          'Cyan', 'c', 'Orange', [0.8500 0.3250 0.0980], 'Purple', [0.4940 0.1840 0.5560]);
        calcColor = colorMap.(selectedColorStr);
        
        shifted_env_x = env_x + current_shift;
        shifted_sticks_x = sticks_x + current_shift;
        
        scaled_env_y = env_y_norm * calc_yscale;
        scaled_sticks_y = sticks_y_norm * calc_yscale;
        
        interp_env = interp1(shifted_env_x, env_y_norm, WaveNo, 'linear', 0);
        corrVal = corr(IRMPD_norm, interp_env, 'Type', 'Pearson');
        pccLabel.Text = sprintf('Current PCC: %.3f', corrVal);
        
        cla(ax); 
        hold(ax, 'on');
        
        hEnv = plot(ax, shifted_env_x, scaled_env_y, 'Color', calcColor, 'LineStyle', '-', 'LineWidth', 1.5);
        stem(ax, shifted_sticks_x, scaled_sticks_y, 'Color', calcColor, 'Marker', 'none', 'LineWidth', 1.5); 
        hIRMPD = plot(ax, WaveNo, IRMPD_norm + y_offset, 'k-', 'LineWidth', 1.5);
        
        xlabel(ax, '\bf\it{Wavenumber (cm^{-1})}');
        ylabel(ax, '\bf\it{Normalized Intensity (Arb. Units)}');
              
        legend(ax, [hIRMPD, hEnv], {'IRMPD Data', 'Gaussian Calculation'}, 'Location', 'northeast');
        grid(ax, 'on');
        
        xlim(ax, [xMinEdit.Value, xMaxEdit.Value]);
        ylim(ax, [yMinEdit.Value, yMaxEdit.Value]);
        ax.YTickLabel = []; 
        
        if ~isempty(summary_text)
            text(ax, 0.5, 0.96, summary_text, 'Units', 'normalized', ...
                'FontSize', 14, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'none', 'EdgeColor', 'none');
        end

        plot_text = sprintf('Theory: %s / %s\nShift: %+.1f cm^{-1}\nOrig. PCC: %.3f\nNew PCC: %.3f', ...
                            func_text, basis_text, current_shift, orig_pcc, corrVal);
        text(ax, 0.02, 0.96, plot_text, 'Units', 'normalized', ...
            'FontSize', 11, 'FontWeight', 'bold', ...
            'BackgroundColor', 'w', 'EdgeColor', 'k', ...
            'VerticalAlignment', 'top');
            
        hold(ax, 'off');
    end

    function optimizeShift()
        bestCorr = -1;
        bestShift = 0;
        shiftRange = -200:0.1:200;
        
        pccLabel.Text = 'Calculating...';
        drawnow;
        
        for shift = shiftRange
            shifted_env_x = env_x + shift;
            interp_env = interp1(shifted_env_x, env_y_norm, WaveNo, 'linear', 0);
            cVal = corr(IRMPD_norm, interp_env, 'Type', 'Pearson');
            
            if cVal > bestCorr
                bestCorr = cVal;
                bestShift = shift;
            end
        end
        
        shiftEdit.Value = bestShift;
        updatePlot();
    end
end
