function string = Dummy_getModelName(model)
% DUMMY_GETMODELNAME return the model string.

% Bimodal-data models
switch model(1)
    case 1; string = 'NoL';
    case 2; string = 'L';
end

end