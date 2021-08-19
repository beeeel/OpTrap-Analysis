%% The Fourier Transform function
function [GFT,gdata_I] = rheoFDFT_irheoFT_vec(t_I,Gint_I,freqpoints,g0,ginf)
% Will's update to code from iRheoFT paper - vectorized calculation of A

gdata_I = zeros(freqpoints,3);

GFT = zeros(freqpoints,3);

wrange = logspace(log10(t_I(1)),log10(t_I(end)),freqpoints); % Freqency Range

for n = 1:freqpoints
    w = wrange(n);
    
%     % Original    
%     for k = 2:length(t_I)
%         A(k) = (((Gint_I(k)-Gint_I(k-1))/(t_I(k)-t_I(k-1)))*...
%             (exp(-1i*w*t_I(k-1))-exp(-1i*w* t_I(k))));
%     end

    % New: vectorized for t_I
    
    A = ( diff(Gint_I) ./ diff(t_I) ) ...
        .* ( exp( -1i .* w .* t_I(1:end-1) ) - exp( -1i .* w .* t_I(2:end) ) );
    
    % Original - uses t_I(2)
%     GFT1 = ((1i*w*g0 + ...
%         ((1-exp(-1i*w* t_I(2)))*((Gint_I(1)-g0)/t_I(2))) + ...
%         ginf*exp(-1i*w*t_I(end))) + ...
%         sum(A)) ...
%         / (1i*w)^2;

    % I think it should use t_I(1), but it doesn't make a visible
    % difference
    GFT1 = ((1i*w*g0 + ...
        ((1-exp(-1i*w* t_I(1)))*((Gint_I(1)-g0)/t_I(1))) + ...
        ginf*exp(-1i*w*t_I(end))) + ...
        sum(A)) ...
        / (1i*w)^2;
    Gstar = GFT1*(1i*w);
    GFT(n,:) = [w real(GFT1) imag(GFT1)];
    gdata_I(n,:) = [w real(Gstar) imag(Gstar)];
end

end

%% Check that vectorized gives same numbers
% figure(1)
% clf
% hold on
% loglog(real(A),'k')
% loglog(imag(A),'k--')
% 
% loglog(real(B),'r')
% loglog(imag(B),'r--')
% 
% figure(2)
% clf
% plot(A-B.')