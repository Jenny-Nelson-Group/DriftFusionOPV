function [all_cell_configs, DP_table, Split_choice] = optimise_slices(n_cells, sliceDatabase)
    % OPTIMISE_SLICES Find the optimal spectral splitting for multi-junction cells.
    % Uses Dynamic Programming to find the maximum power for 1 to N sections.
    %
    % INPUTS:
    %   n_cells       - (Integer) Maximum number of sections/junctions to solve for.
    %   sliceDatabase - (3D Struct) Dimensions [gaps x start_lam x end_lam] or
    %                   (4D Struct) Dimensions [gaps x start_lam x end_lam x materials]
    %                   Must contain field .Pmax, .gap, .l_min, .l_max, .Voc, .Jsc
    %
    % OUTPUTS:
    %   all_configs   - (1xN Cell) Optimal configuration for each total section count.
    %                   Each cell contains .total_power, .gaps, .split_wavelengths, 
    %                   and .section_powers.
    %   DP_table      - (NxL Matrix) Max power for 'n' sections starting at index 'i'.
    %   Split_choice  - (NxL Matrix) Optimal split index 'k' for the current subproblem.
    %
    % Example:
    %   results = optimiseRainbow_fast_allN(5, slice_database);
    %   disp(results{3}.total_power); % Power for the 3-junction solution
    
    %% 0. Validate input format
    [num_gaps, num_lambdas, ~, num_materials] = size(sliceDatabase);

    if n_cells >= num_lambdas
        error('n_cells (%d) must be less than num_lambdas (%d).', n_cells, num_lambdas);
    end

    %% 1. Pre-Extract Power Matrix
    % P_slice(i, j) holds the max power from lambda index i to j across all gaps
    % Gap_slice(i, j) holds the specific bandgap index that produced that power

    P_slice = zeros(num_lambdas, num_lambdas);
    Gap_slice = zeros(num_lambdas, num_lambdas);
    Material_slice = zeros(num_lambdas, num_lambdas);
    
    % Find the optimal gap for each slice
    for i = 1:num_lambdas-1
        for j = i+1:num_lambdas
            slice_powers_block = [sliceDatabase(:, i, j, :).Pmax];
            slice_powers_matrix = reshape(slice_powers_block, num_gaps, num_materials);
            
            % Find global maximum in the Gap vs Material grid
            [max_p, linear_idx] = max(slice_powers_matrix(:));
            [best_g, best_m] = ind2sub(size(slice_powers_matrix), linear_idx);

            P_slice(i, j) = max_p;
            Gap_slice(i, j) = best_g;
            Material_slice(i, j) = best_m;
        end
    end

    %% 2. Initialize DP Tables
    % DP_table(n, i) = max power using 'n' sections starting at lambda index 'i'
    % Split_choice(n, i) = the lambda index 'k' we chose to end the first section
   
    DP_table = zeros(n_cells, num_lambdas);
    Split_choice = zeros(n_cells, num_lambdas);
    
    % Base Case: 1 section remaining.
    for i = 1:num_lambdas-1
        DP_table(1, i) = P_slice(i, num_lambdas);
        Split_choice(1, i) = num_lambdas;
    end
    
    %% 3. Execute DP Transitions
    % We can only start up to index (num_lambdas - n) to leave room for remaining sections
    for n = 2:n_cells
        for i = 1:(num_lambdas - n)
            max_val = -1;
            best_split = -1;

            % Check all calid split points "k" for the current section
            for k = (i+1):(num_lambdas - n + 1)
                current_val = P_slice(i, k) + DP_table(n-1, k);
                if current_val > max_val
                    max_val = current_val;
                    best_split = k;
                end
            end

            % Store the optimal choice for this subproblem
            DP_table(n, i) = max_val;
            Split_choice(n, i) = best_split;
        end
    end
    
    %% 4. Reconstruct ALL Configurations
    % reconstruct the optimal configuration for each number of cells
    all_cell_configs = cell(1, n_cells);
    
    for n_total = 1:n_cells
        % Pre-allocate the top-level total power
        cell_config = struct('total_power', DP_table(n_total, 1));
        current_idx = 1;
        
        for step = 1:n_total
            n_remaining = n_total - step + 1;
            next_idx = Split_choice(n_remaining, current_idx);
            
            % Identify and pull the exact pre-calculated entry
            g_idx = Gap_slice(current_idx, next_idx);
            m_idx = Material_slice(current_idx, next_idx);

            best_entry = sliceDatabase(g_idx, current_idx, next_idx, m_idx);
            
            % Assign values directly (No pre-allocation of zeros needed)
            cell_config.material_idx(step)        = m_idx;
            cell_config.gaps(step)                = best_entry.gap;
            cell_config.split_wavelengths(step)   = best_entry.l_min;
            cell_config.split_wavelengths(step+1) = best_entry.l_max; % Overwrites l_min of next cell with l_max of current
            cell_config.section_powers(step)      = best_entry.Pmax;
            cell_config.section_Voc(step)         = best_entry.Voc;
            cell_config.section_Jsc(step)         = best_entry.Jsc;
            cell_config.section_gen_abs(step)     = best_entry.gen_abs;
            cell_config.section_gen_frac(step)    = best_entry.gen_frac;
           
            current_idx = next_idx;
        end
        all_cell_configs{n_total} = cell_config;
    end
end
