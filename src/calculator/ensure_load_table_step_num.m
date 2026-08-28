function loadTable = ensure_load_table_step_num(loadTable)
    % ensure_load_table_step_num Keep load-step numbering aligned with table rows.
    if istable(loadTable) && (~ismember('step_num', loadTable.Properties.VariableNames) || ...
            numel(loadTable.step_num) ~= height(loadTable) || ...
            ~isequal(loadTable.step_num(:), (1:height(loadTable))'))
        loadTable.step_num = (1:height(loadTable))';
        loadTable = movevars(loadTable, 'step_num', 'Before', 1);
    end
end
