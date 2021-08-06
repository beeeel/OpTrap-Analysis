%% The Fourier Transform function
function [GFT,gdata_I] = rheoFDFT_irheoFT(t_I,Gint_I,freqpoints,g0,ginf)

gdata_I = zeros(freqpoints,3);

GFT = zeros(freqpoints,3);

wrange = logspace(-2,2,freqpoints); % Freqency Range

A = zeros(1,length(t_I)-1);

for n = 1:freqpoints
    w = wrange(n);
    
    for k = 2:length(t_I)
        A(k-1) = (((Gint_I(k)-Gint_I(k-1))/(t_I(k)-t_I(k-1)))*...
            (exp(-1i*w*t_I(k-1))-exp(-1i*w* t_I(k))));
    end
    
    GFT1 = ((1i*w*g0+((1-exp(-1i*w* t_I(2)))*((Gint_I(1)-g0)/t_I(2)))+...
        ginf*exp(-1i*w*t_I(end)))+sum(A))/(1i*w)^2;
    
    Gstar = GFT1*(1i*w);
    GFT(n,:) = [w real(GFT1) imag(GFT1)];
    gdata_I(n,:) = [w real(Gstar) imag(Gstar)];
end

end
