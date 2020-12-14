function stiffness = calcStiffness(centresVector, varargin)
%% Calculate trap stiffness from centres vector
% Either supply centresVec in units m, or conversion factor in m/[unit]

% centresVector needs to be a real nonempty vector
validateattributes(centresVector,{'numeric'},...
    {'nonempty','vector','real'},'calcStiffness','centresVector',1)

% CentresVectorM needs to be in units of m. Default to assuming it is
% already in m, and if a conversion factor is supplied, use it.
conversionFactor = 1;
if nargin == 2
    conversionFactor = varargin{1};
    validateattributes(conversionFactor,{'numeric'},...
        {'nonempty','positive','scalar','real'},'calcStiffness','conversionFactor',2)
elseif nargin > 2
    error('TOO MANY NARGINS!!!!')
end
centresVectorM = conversionFactor .* centresVector;

Kb = 1.38064852e-23; % Boltzmann constant
T = 273 + 20; % Assume trap at room temperature - maybe I should model this?
% From Sarshar et al 2014 eq 4. This gives stiffness in N/m.
stiffness = Kb * T ./ var(centresVectorM);