function p = projectionparameters(approxtype, splineorder)

% PROJECTIONPARAMETERS creates a structure containing the approximation type and the spline order (if using splines) for the collocation method
% It can be used with just the first input (approxtype). In that case, splineorder is assumed to be [] (default value in CompEcon)

p.approxtype = approxtype;
if nargin==1
    p.splineorder = [];
else
    p.splineorder = splineorder;
end

end

