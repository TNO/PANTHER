function loadTable = ensure_load_table_step_num(loadTable)
    % ensure_load_table_step_num Add sequential step_num to legacy load tables.
    if istable(loadTable) && ~ismember('step_num', loadTable.Properties.VariableNames)
        loadTable.step_num = (1:height(loadTable))';
        loadTable = movevars(loadTable, 'step_num', 'Before', 1);
    end
end
