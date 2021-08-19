function [chi1, chi2] = rheoFDFT_Nishi(tau, msd, omega)
%Y = rheoFDFT(tau, msd, omega)
% Compute finite discrete fourier transform using Nishi et al (2018), based
% on python code from Tassieri (2018).

% Input checking
k = length(tau);
try
    validateattributes(msd, {'numeric'}, {'nrows', k, 'ncols',1})
catch
    error('Received msd with size (%i %i), and tau with %i rows', size(msd), size(tau,1))
end

% If first tau and all first msd are nonzero, correct the user
if tau(1) ~= 0 && min(msd(1,:) ~= 0)
    tau = [0; tau];
    msd = [0; msd];
% If user is inconsistent, throw toys out of the pram
elseif tau(1)~= 0 || min(msd(1,:) ~= 0)
    error('Expected first (tau msd) to be (0 0) but one was not')
end
    
% Setup
m = length(omega);
dt = tau(2);

n = size(msd, 2);

% Response function in time domain
chit = zeros(k, n);
% Real and imaginary response function in frequency domain
chi1 = zeros(m, n);
chi2 = zeros(m, n);
% 
l = floor((k-4)/2);

% Not 100% sure where this comes from. Nishi 2018 has kT in several eqs,
% but Tassieri (2018) is the source for the 0.5.
kBT = 0.5;

%% Numerical derivative
for i = 1:2
    chit(i) = 1/12 * (-25 * msd(i) + 48 * msd(i + 1) - 36 * msd(i + 2) + 16 * msd(i + 3) - 3 * msd(i + 4)) / ( kBT * 2 * dt );
end
for i = 3:k-2
    chit(i) = 1/12 * (msd(i - 2) - 8 * msd(i - 1) + 8 * msd(i + 1) - msd(i + 2)) / (kBT * 2 * dt);
end
for i = k-1:k
    chit(i) = 1/2 * (msd(i-2) - 4 * msd(i - 1) + 3 * msd(i)) / (kBT * 2 * dt);
end

%% Fourier transform

% New method: vectorized
for j = 1
    % Due to operator precedence this is the same as /3*dt, but clearer
    chi1 = chi1 + chit(j) * cos(omega * tau(j)) * dt/3;
    chi2 = chi2 + chit(j) * sin(omega * tau(j)) * dt/3;
end
for j = 2 : l+1
    chi2 = chi2 + chit(2 * j - 1) * sin(omega * tau(2 * j - 1) ) * 4 * dt / 3 + ...
        + chit(2 * j) * sin(omega * tau(2 * j) ) * 2 * dt / 3;
    chi1 = chi1 + chit(2 * j - 1) * cos(omega * tau(2 * j - 1) ) * 4 * dt / 3 + ...
        + chit(2 * j) * cos(omega * tau(2 * j) ) * 2 * dt / 3;
end
for j = l+2
    chi2 = chi2 + chit(2 * j - 1) * sin(omega * tau(2 * j - 1) ) * 4 * dt / 3 + ...
        + chit(2 * j) * sin(omega * tau(2 * j) ) * dt / 3;
    chi1 = chi1 + chit(2 * j - 1) * cos(omega * tau(2 * j - 1) ) * 4 * dt / 3 + ...
        + chit(2 * j) * cos(omega * tau(2 * j) ) * dt / 3;
end
% % Original method: nested loops
% for i = 1:m
%     if rem(i, 10) == 0
%         fprintf('\r%i / %i', i, m)
%     end
%     for j = 1
%         % Due to operator precedence this is the same as /3*dt, but clearer
%         chi1(i) = chi1(i) + chit(j) * cos(omega(i) * tau(j)) * dt/3; 
%         chi2(i) = chi2(i) + chit(j) * sin(omega(i) * tau(j)) * dt/3;
%     end
%     for j = 2 : l+1
%         chi2(i) = chi2(i) + chit(2 * j - 1) * sin(omega(i) * tau(2 * j - 1) ) * 4 * dt / 3 + ...
%             + chit(2 * j) * sin(omega(i) * tau(2 * j) ) * 2 * dt / 3;
%         chi1(i) = chi1(i) + chit(2 * j - 1) * cos(omega(i) * tau(2 * j - 1) ) * 4 * dt / 3 + ...
%             + chit(2 * j) * cos(omega(i) * tau(2 * j) ) * 2 * dt / 3;
%     end
%     for j = l+2
%             chi2(i) = chi2(i) + chit(2 * j - 1) * sin(omega(i) * tau(2 * j - 1) ) * 4 * dt / 3 + ...
%                 + chit(2 * j) * sin(omega(i) * tau(2 * j) ) * dt / 3;
%             chi1(i) = chi1(i) + chit(2 * j - 1) * cos(omega(i) * tau(2 * j - 1) ) * 4 * dt / 3 + ...
%                 + chit(2 * j) * cos(omega(i) * tau(2 * j) ) * dt / 3;
%     end
% end
% fprintf('\rDONE!\n')


end
