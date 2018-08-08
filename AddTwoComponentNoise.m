%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Rotational electrical 
% impedance tomography using electrodes with limited boundary coverage
% provides window for multimodal sensing".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional noise of two components as described in the article.
function data = AddTwoComponentNoise( data_NF, common_err, component_err )
    data = data_NF;
    for i = 1:numel(data)
        data(i) = data(i) .* ( 1 + plusminusErrorTerm(data(i)*component_err) ...
            + plusminusErrorTerm(common_err) );
    end
end

function err = plusminusErrorTerm( level )
    err = level*(2*rand - 1);
end
