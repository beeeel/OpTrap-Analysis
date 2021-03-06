%% Playing with fourier for unwrapping
Angles = linspace(0,2*pi,360);
Eqn = @(a, b, phi, x) (a.^2 * cos(x + phi).^2 + b.^2 * sin(x + phi).^2);

%% Recover phase information
Maj = 100;
Min = 50;
Data = Eqn(Maj, Min, Angles', Angles);

Freq = fft(Data,[],2);

Phase = unwrap(angle(Freq(:,3)))/2;

figure(3)
hold on
plot(1:360,Angles,'xk')
plot(1:360,Phase,'r','LineWidth',3)
hold off
%% Recover amplitude information
Maj = 50:100;
Min = 50;
Data = Eqn(Maj', Min', 0, Angles);

Freq = fft(Data,[],2);

Sum = 2*(abs(Freq(:,1)));
Dif = abs(Freq(:,3))/90;

Ys = linspace(Amp(1),Amp(end),51)';

figure(3)
clf
hold on
plot(1:51, (Sum+Dif)/2, 'xk')
%plot(1:51, (Sum-Dif)/2, 'or')
plot(1:51, Maj.^2, 'k-')
plot(1:51, Min.^2, 'r-')
%plot(1:51, Amp,'r','LineWidth',3)



