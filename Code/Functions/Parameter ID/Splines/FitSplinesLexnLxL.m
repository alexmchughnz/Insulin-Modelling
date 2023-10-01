function [P, A, b, basisSplines] = FitSplinesLexnLxL(P, splineOptions)

    CONST = Constants();
    GC = P.parameters.GC;
    SC = P.parameters.SC;
    minnL = 0;
    numFixedParameters = 2; % For xL and Lex
    
    
    %% Splines    
    defaultSplineOptions.constrain = true;
    defaultSplineOptions.knotType = "location";
    defaultSplineOptions.knots = P.data.I.time;
    defaultSplineOptions.order = 3;
    defaultSplineOptions.maxRate = 0.001;

    fields = fieldnames(defaultSplineOptions);
    for ii = 1:numel(fields)
        field = fields{ii};
        if ~isfield(splineOptions, field)
            splineOptions.(field) = defaultSplineOptions.(field);
        end
    end
    P.results.splineOptions = splineOptions;
    
    % Collect basis functions for splines.
    [P, basisSplines, tExtended, allKnots, dataKnots] = MakeSplineBasisFunctions(P, splineOptions);
    numTotalSplines = size(basisSplines, CONST.COLUMNDIM);
    numTotalParameters = numFixedParameters + numTotalSplines;
    
    
    %% Setup
    tArray = P.results.tArray;
    tInts = dataKnots';
    iiInts = SearchArray(tInts, tArray);
    
    % Plasma Insulin
    [tI, vI] = GetData(P.data.I); % [mU/L]
    ppI = griddedInterpolant(tI, vI);  % [mU/L]
    I = ppI(tArray);
    
    % Interstital Insulin
    Q = GetAnalyticalInterstitialInsulin(I, P);
    
    % Endogenous Secretion
    Uen = P.results.Uen;
    
    % Exogenous Insulin
    Uex = P.results.Uex(P);


    %% Get Coefficients
    % Analytical solutions for ISC/Lex and Q_Local/Lex.
    % Assuming ISC(0) and QLocal(0) = 0.
    cISC = exp(-SC.ks2*tArray) .* cumtrapz(tArray, exp(SC.ks2*tArray).*Uex); % factor of Lex to be IDed; not included
    cQLocal = SC.ks2*exp(-(SC.ks3+SC.kdi)*tArray) .* cumtrapz(tArray, exp((SC.ks3+SC.kdi)*tArray).*cISC);
        
    % Consider:
    % dI/dt = cL*Lex - cn*nL - cx*xL + kU*Uen - kI*I - kIQ*(I-Q)
    %   with nL = sum(nLWeight_i * shape_i).
    
    % Let cWeight_i = shape_i * cn.
    % We can express equation as:
    % dI/dt = cL*Lex - sum(cWeight_i * nLWeight_i) - cx*xL + kU*Uen - kI*I - kIQ*(I-Q)
    cn = I./(1 + GC.alphaI*I);
    cWeights =  basisSplines .* cn;
    cx = Uen/GC.VI;
    cL = SC.ks3/GC.VI * cQLocal;

    kU = 1/GC.VI;
    kI = GC.nK;
    kIQ = GC.nI./GC.VI;
    
    %% Integrate I Equation
    % I(t) - I(t0) = Lex*int{cL} - int{sum(cWeight_i * nLWeight_i)} - xL*int{cx} + kU*int{Uen} - kI*int{I} - kIQ*int{I-Q}
    % -int{cWeights}.*nLWeights - xL*int{cx} + Lex*int{cL} = I(t) - I(t0) - kU*int{Uen} + kI*int{I} + kIQ*int{I-Q} := C
    CWeights = -cumtrapz(tArray, cWeights);
    CX = -cumtrapz(tArray, cx);
    CL = cumtrapz(tArray, cL);
    
    intUTerm = kU*cumtrapz(tArray, Uen);
    intITerm = kI*cumtrapz(tArray, I);
    intIQTerm = kIQ*cumtrapz(tArray, I-Q);
    
    I0 = I(1) * ones(size(I));
    RHS = [I -I0 -intUTerm intITerm intIQTerm];
    C = sum(RHS, CONST.COLUMNDIM);
    
    %% Assemble MLR System
    % Extract values at measurement points.
    vCWeights = CWeights(iiInts, :);
    vCX = CX(iiInts, :);
    vCL = CL(iiInts, :);
    vC = C(iiInts);
    
    % Find time width of each integral wedge.
    iiFirst = 1:numel(vC) - 1;
    iiSecond = 2:numel(vC);
    
    t1 = tInts(iiFirst);
    t2 = tInts(iiSecond);
    dt = t2 - t1;
    
    % Assemble MLR system by evaluating integral between
    % sample points, and normalising by integral width (dt):
    % [CW_1(t) CW_2(t) ... CW_M(t)] * (nLW1; nLW2; ... nLWM) = [C(t)] @ measurement times only
    CWValues = (vCWeights(iiSecond,:) - vCWeights(iiFirst,:)) ./ dt;
    CXValues = (vCX(iiSecond,:) - vCX(iiFirst,:)) ./ dt;
    CLValues = (vCL(iiSecond,:) - vCL(iiFirst,:)) ./ dt;
    CValues = (vC(iiSecond) - vC(iiFirst)) ./ dt;
    
    A(:,1) = CLValues;
    A(:,2) = CXValues;
    A(:,3:numTotalParameters) = CWValues;
    b = CValues;
    
    %% Set up Splines
    % Enforce constraints on nL spline weights.
    numConstraints = numTotalSplines - 1;
    [tG, vG] = GetData(P.data.G);
    ppG = griddedInterpolant(tG, vG);
    
    nLDiffMatrix = zeros(numConstraints, numConstraints+1);
    for ii = (1+numFixedParameters):numConstraints
        % Set up difference matrix.
        nLDiffMatrix(ii, ii) = -1;
        nLDiffMatrix(ii, ii+1) = 1;
    end
    
    % Directional constraint - delta(nL) should always be opposite to delta(G).
    % Expand delta(G) pattern by at edges to constrain "extra" splines.
    dataSplineIndex = max([1 splineOptions.order]);
    dataKnots = allKnots(dataSplineIndex + (0:numConstraints))'; % Take off extra knots, keep those that define data range.
    GDiffPattern = diff(ppG(dataKnots));
    
    % nL should decrease if G is increasing, thus
    % sign(G(ii+1) - G(ii)) * (nL(ii+1) - nL(ii)) < 0
    nLDirectionA = sign(GDiffPattern) .* nLDiffMatrix;
    nLDirectionb = zeros(numConstraints, 1);
    
    % Change over time constraint - delta(nL) should be less than
    % 0.001 min^-2 (Caumo, 2007).
    maxnLRate = splineOptions.maxRate;
    tDeltas = diff(dataKnots);
    maxnLDeltas = maxnLRate * tDeltas;
    % Absolute change in nL must be less than max change, thus
    % abs(nL(ii+1) - nL(ii))  < maxnLDeltas(ii)
    % sign(nL(ii+1) - nL(ii)) * (nL(ii+1) - nL(ii)) < maxnLDeltas(ii)
    % -sign(G(ii+1) - G(ii)) * (nL(ii+1) - nL(ii)) < maxnLDeltas(ii)
    nLChangeA = -sign(GDiffPattern) .* nLDiffMatrix ;
    nLChangeb = maxnLDeltas;
    
    % Place constraints on splines.
    % "A" structure is [fixedParams, extraSplines, dataSplines, extraSplines].
    iiSplines = numFixedParameters + [1:numTotalSplines];
    
    AConstraint = zeros(2*numConstraints, numTotalParameters);
    AConstraint(:, iiSplines) = [nLDirectionA; nLChangeA];
    
    bConstraint = [nLDirectionb; nLChangeb];
    
    
    % Solve using linear solver.
    lb = minnL * ones(1, numTotalParameters);
    if splineOptions.constrain
        x = lsqlin(A, b, AConstraint, bConstraint, [], [], lb);
    else
        x = lsqlin(A, b, [], [], [], [], lb);
    end
    
    P.results.JLK = x(1);
    P.results.xL = x(2);
    nLWeights = x(3:end);
    P.results.nL = basisSplines * nLWeights;
    
    
    %% Plotting
    plotvars.basisSplines = basisSplines;
    plotvars.nLWeights = nLWeights;
    
    P = MakePlots(P, plotvars);
    
    end
    
    
    function P = MakePlots(P, plotvars)
    tag = "FitSplines";
    
    %% Splines
        P = AddFigure(P, tag, "Splines");
        
        % Plot nL.
        plt = plot(P.results.tArray, P.results.nL, 'b');
        plt.DisplayName = "nL";
        
        % Plot fitted splines.
        plot(P.results.tArray, plotvars.basisSplines .* plotvars.nLWeights', '--', ...
            'LineWidth', 1, 'HandleVisibility', 'off');
        
        ylim([0 Inf])
        
        xlabel("Time [min]")
        ylabel("nL [1/min]")
    end    

