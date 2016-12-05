function xArgOut = ChanVese(sContext, dImg, lMask, hDrawFcn, cParams)

switch sContext    %This distinguishes the contexts in which the plugin function is called

    % ---------------------------------------------------------------------
    % The tooltip context: Simply return a tooltip string
    case 'tooltip', xArgOut = 'Chan-Vese Segmentation';
    % End of the tooltip path
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % The input context: Return input type: 'roi', 'seed', or 'start'
    case 'input', xArgOut = 'roi';
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % The parameter context: Return a struct with a list of parameters
    case 'params'
        xArgOut = struct( ...
            'Name',         {'Lambda 1', 'Lambda 2',    'Mu', 'Delta T', 'Iterations'},...
            'Type',         {  'double',   'double','double',  'double',    'integer'}, ...
            'Min',          {         0,          0,       0,         0,            1}, ...
            'Max',          {       100,        100,     inf,       inf,           50}, ...
            'Def',          {         1,          1,       1,         4,            5});
    % End of the parameter path
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % The execution context
    case 'exe'
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Check the input arguments
        dImg = double(dImg);
        dImg = dImg./max(dImg(:));
        
        dLambda1 = cParams{1};
        dLambda2 = cParams{2};
        dNu      = cParams{3};
        dDeltaT  = cParams{4};
        iNIter   = cParams{5};
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        dPhi = double(lMask);
        dPhi = fCVReinit(dPhi - 0.5);
        
        lMask = dPhi > 0;
        if isa(hDrawFcn, 'function_handle'), hDrawFcn(lMask); end

        for iI = 1:iNIter
            dPhi = fCVIterate(dImg, dPhi, dNu, dLambda1, dLambda2, dDeltaT);
            lMask = dPhi > 0;
            if isa(hDrawFcn, 'function_handle'), hDrawFcn(lMask); end
        end
        xArgOut = lMask;
    % end of the execution path
    % ---------------------------------------------------------------------
    
    otherwise, xArgOut = [];
        
end
% of the switch statement
% -------------------------------------------------------------------------



    function dPhi = fCVIterate(dImg, dPhi, dNu, dLambda1, dLambda2, dDeltaT)
        dMuIn  = mean(dImg(dPhi >  0));
        dMuOut = mean(dImg(dPhi <= 0));
        dPhi = fCVReinit(dPhi);
        dPhi = dPhi + dDeltaT.*100.*(dLambda1.*(dImg - dMuOut).^2 - dLambda2.*(dImg - dMuIn).^2);
        dPhi = fCVReinit(dPhi);
        dPhi = ac_div_AOS_3D_dll(dPhi, ones(size(dPhi)), dDeltaT.*dNu);
    end



    function dPhi = fCVReinit(dPhi)
        dPhi0 = zy_binary_boundary_detection(uint8(dPhi > 0));
        dPhi0 = padarray(dPhi0, [1 1 1]);
        dPhi0 = ac_distance_transform_3d(dPhi0);
        dPhi0 = dPhi0(2:end - 1, 2:end - 1, 2:end - 1);
        dPhi  = dPhi0.*sign(dPhi);
    end


end