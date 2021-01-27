function data = bead_fft_scaled(data, doPlots)
%% data = bead_fft_scaled(data, doPlots)
% Calculate Fourier transform of centres data in both directions, store
% one-sided spectra in data struct along with frequency vector in Hz

n_points = diff(data.opts.cropT)+1;

for direction = 'xy'
    % Get the FFT
    X = fft(data.pro.([ direction 'CentresM']));
    P = abs(X/n_points);
    P = P(1:end/2+1);
    P(2:end-1) = 2*P(2:end-1);
    
    data.pro.([direction 'fftM']) = P;
end

% Frequency in index (k+1) is k cycles per whole dataset, so convert to Hz
data.pro.fftFreqHz = (0:ceil((n_points-1)/2))./diff(data.raw.timeVecMs([1 end])*1e-3);

if doPlots
    fh = figure;
    fh.Name = data.fName;
    
    subplot(2,1,1)
    plot(data.pro.fftFreqHz, 1e9 * data.pro.xfftM)
    ylim([0 1])
    
    title('X centres frequency spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Intensity (nm)')
    
    subplot(2,1,2)
    plot(data.pro.fftFreqHz, 1e9 * data.pro.yfftM)
    ylim([0 1])
    
    title('Y centres frequency spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Intensity (nm)')
    drawnow
end