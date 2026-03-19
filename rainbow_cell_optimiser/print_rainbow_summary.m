function print_rainbow_summary(all_configs)
    % PRINT_RAINBOW_SUMMARY Detailed report including Voc, Jsc, and Generation data.
    
    total_N = length(all_configs);
    P_in = 100; % Reference input power (mW/cm^2)
    
    fprintf('\n============================================================================================\n');
    fprintf('                          RAINBOW CELL OPTIMIZATION: MULTI-SECTION SUMMARY\n');
    fprintf('============================================================================================\n');

    for n = 1:total_N
        config = all_configs{n};
        total_PCE = (config.total_power / P_in) * 100;
        
        fprintf('\n[ Configuration: %d Sections ]\n', n);
        fprintf('Total Power: %6.2f mW/cm^2 | Total PCE: %5.2f %%\n', ...
                config.total_power, total_PCE);
        fprintf('--------------------------------------------------------------------------------------------\n');
        
        % Expanded Header
        fprintf('%-4s %-14s %-8s %-8s %-8s %-6s %-10s %-8s\n', ...
                '#', 'Range (nm)', 'Eg(eV)', 'Voc(V)', 'Jsc(mA)', 'FF(%)', 'PCE(%)', 'Gen(abs)');
        
        for s = 1:n
            w1  = config.split_wavelengths(s);
            w2  = config.split_wavelengths(s+1);
            gap = config.gaps(s);
            voc = config.section_Voc(s);
            jsc = config.section_Jsc(s);
            pwr = config.section_powers(s);
            pce = (pwr / P_in) * 100;
            gen = config.section_gen_abs(s);
                      
            fprintf('%-4d %4.0f-%4.0f nm  %6.2f   %6.3f   %7.2f  %5.1f  %6.2f%%   %8.2e\n', ...
                    s, w1, w2, gap, voc, jsc, 0, pce, gen);
        end
        fprintf('--------------------------------------------------------------------------------------------\n');
    end
end