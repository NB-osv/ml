% Filter specification 1
%Lowpass filter

%%
%Butterworth
Amax=0.1;
Amin=55;
Wc=2*pi*50000;
Ws=2*pi*60000;
[N, Wn]=buttord(Wc,Ws,Amax,Amin,'s');
[Z,P,G]=buttap(N); %Normalized poles
epsilon=sqrt(10^(0.1*Amax)-1);
w0p=Wc*epsilon^(-1/N);
P=P*w0p;%Denormalized poles
G=G*w0p^N;

[b, a] = zp2tf(Z, P, G); %% EDIT THIS: Use appropriate Z, P and G!
w = linspace(0, 2*Ws, 2000); %% EDIT THIS for bandpass or bandstop!
w_pass = linspace(0, Wc, 2000); %% EDIT THIS for bandpass or bandstop!
w_stop = linspace(Ws, 2*Ws, 2000); %% EDIT THIS for bandpass or bandstop!
Hb = freqs(b, a, w);
Hb_pass = freqs(b, a, w_pass);
Hb_stop = freqs(b, a, w_stop);
hFig = figure(1); %% EDIT THIS: Use different numbers!
set(hFig, 'Name', 'Butterworth') %% EDIT THIS: Approximation Type!
subplot(221)
zplane(Z/1e3, P/1e3); %% EDIT THIS: Use appropriate Z and P!
ylabel('j\omega [\times 10^3]')
xlabel('\sigma [\times 10^3]')
title('{\it s}-plane')
subplot(222)
plot(w/2/pi/1e3, 20*log10(abs(Hb)));
ylabel('dB')
xlabel('kHz')
title('Magnitude response')
subplot(223)
plot(w_pass/2/pi/1e3, 20*log10(abs(Hb_pass)));
ylabel('dB')
xlabel('kHz')
title('Passband mag. resp.')
subplot(224)
plot(w_stop/2/pi/1e3, 20*log10(abs(Hb_stop)));
ylabel('dB')
xlabel('kHz')
title('Stopband mag. resp.')

fstart = 0; %% EDIT THIS if required!
fend = 60e3; %% EDIT THIS for bandpass or bandstop!
w_vec = 2*pi*linspace(fstart, fend, 1000);
fsample = 4*fend;

gd_b = groupdelay(Z, P, G, w_vec, fsample); 

figure(2); %% EDIT THIS: May require different number!
clf
plot(w_vec/2/pi/1e3, gd_b/1e-3, '--r', 'MarkerSize', 4)
ylabel('Group delay [ms]')
xlabel('Frequency [kHz]')


%%
%Chebyshev-I
Amax=0.1;
Amin=55;
Wc=2*pi*50000;
Ws=2*pi*60000;
[N, Wn]=cheb1ord(Wc,Ws,Amax,Amin,'s');
[Z,P,G]=cheb1ap(N,Amax);
P=P*Wc;
G=G*Wc^N;

[b, a] = zp2tf(Z, P, G); %% EDIT THIS: Use appropriate Z, P and G!
w = linspace(0, 2*Ws, 2000); %% EDIT THIS for bandpass or bandstop!
w_pass = linspace(0, Wc, 2000); %% EDIT THIS for bandpass or bandstop!
w_stop = linspace(Ws, 2*Ws, 2000); %% EDIT THIS for bandpass or bandstop!
Hb = freqs(b, a, w);
Hb_pass = freqs(b, a, w_pass);
Hb_stop = freqs(b, a, w_stop);
hFig = figure(1); %% EDIT THIS: Use different numbers!
set(hFig, 'Name', 'Chebyshev-I') %% EDIT THIS: Approximation Type!
subplot(221)
zplane(Z/1e3, P/1e3); %% EDIT THIS: Use appropriate Z and P!
ylabel('j\omega [\times 10^3]')
xlabel('\sigma [\times 10^3]')
title('{\it s}-plane')
subplot(222)
plot(w/2/pi/1e3, 20*log10(abs(Hb)));
ylabel('dB')
xlabel('kHz')
title('Magnitude response')
subplot(223)
plot(w_pass/2/pi/1e3, 20*log10(abs(Hb_pass)));
ylabel('dB')
xlabel('kHz')
title('Passband mag. resp.')
subplot(224)
plot(w_stop/2/pi/1e3, 20*log10(abs(Hb_stop)));
ylabel('dB')
xlabel('kHz')
title('Stopband mag. resp.')


fstart = 0; %% EDIT THIS if required!
fend = 60e3; %% EDIT THIS for bandpass or bandstop!
w_vec = 2*pi*linspace(fstart, fend, 1000);
fsample = 4*fend;

gd_c1 = groupdelay(Z, P, G, w_vec, fsample); 

figure(2); %% EDIT THIS: May require different number!
clf
plot(w_vec/2/pi/1e3, gd_c1/1e-3, '-k', 'MarkerSize', 4)

ylabel('Group delay [ms]')
xlabel('Frequency [kHz]')


%%
%Chebyshev-II
Amax=0.1;
Amin=55;
Wc=2*pi*50000;
Ws=2*pi*60000;
[N, Wn]=cheb2ord(Wc,Ws,Amax,Amin,'s');
[Z,P,G]=cheb2ap(N,Amin);
P=P*Ws;
Z=Z*Ws;
G=G*Ws^(length(P)-length(Z));


[b, a] = zp2tf(Z, P, G); %% EDIT THIS: Use appropriate Z, P and G!
w = linspace(0, 2*Ws, 2000); %% EDIT THIS for bandpass or bandstop!
w_pass = linspace(0, Wc, 2000); %% EDIT THIS for bandpass or bandstop!
w_stop = linspace(Ws, 2*Ws, 2000); %% EDIT THIS for bandpass or bandstop!
Hb = freqs(b, a, w);
Hb_pass = freqs(b, a, w_pass);
Hb_stop = freqs(b, a, w_stop);
hFig = figure(1); %% EDIT THIS: Use different numbers!
set(hFig, 'Name', 'Chebyshev-II') %% EDIT THIS: Approximation Type!
subplot(221)
zplane(Z/1e3, P/1e3); %% EDIT THIS: Use appropriate Z and P!
ylabel('j\omega [\times 10^3]')
xlabel('\sigma [\times 10^3]')
title('{\it s}-plane')
subplot(222)
plot(w/2/pi/1e3, 20*log10(abs(Hb)));
ylabel('dB')
xlabel('kHz')
title('Magnitude response')
subplot(223)
plot(w_pass/2/pi/1e3, 20*log10(abs(Hb_pass)));
ylabel('dB')
xlabel('kHz')
title('Passband mag. resp.')
subplot(224)
plot(w_stop/2/pi/1e3, 20*log10(abs(Hb_stop)));
ylabel('dB')
xlabel('kHz')
title('Stopband mag. resp.')

fstart = 0; %% EDIT THIS if required!
fend = 60e3; %% EDIT THIS for bandpass or bandstop!
w_vec = 2*pi*linspace(fstart, fend, 1000);
fsample = 4*fend;

gd_c2 = groupdelay(Z, P, G, w_vec, fsample);
figure(2); %% EDIT THIS: May require different number!
clf
plot(w_vec/2/pi/1e3, gd_c2/1e-3, '-.b', 'MarkerSize', 4)
ylabel('Group delay [ms]')
xlabel('Frequency [kHz]')


%%
%Cauer
Amax=0.1;
Amin=55;
Wc=2*pi*50000;
Ws=2*pi*60000;
[N, Wn]=ellipord(Wc,Ws,Amax,Amin,'s');
[Z,P,G]=ellipap(N,Amax,Amin);
P=P*Wc;
Z=Z*Wc;
G=G*Ws^(length(P)-length(Z));

[b, a] = zp2tf(Z, P, G); 
w = linspace(0, 2*Ws, 2000); 
w_pass = linspace(0, Wc, 2000); 
w_stop = linspace(Ws, 2*Ws, 2000); 
Hb = freqs(b, a, w);
Hb_pass = freqs(b, a, w_pass);
Hb_stop = freqs(b, a, w_stop);
hFig = figure(1); %% EDIT THIS: Use different numbers!
set(hFig, 'Name', 'Cauer') %% EDIT THIS: Approximation Type!
subplot(221)
zplane(Z/1e3, P/1e3); %% EDIT THIS: Use appropriate Z and P!
ylabel('j\omega [\times 10^3]')
xlabel('\sigma [\times 10^3]')
title('{\it s}-plane')
subplot(222)
plot(w/2/pi/1e3, 20*log10(abs(Hb)));
ylabel('dB')
xlabel('kHz')
title('Magnitude response')
subplot(223)
plot(w_pass/2/pi/1e3, 20*log10(abs(Hb_pass)));
ylabel('dB')
xlabel('kHz')
title('Passband mag. resp.')
subplot(224)
plot(w_stop/2/pi/1e3, 20*log10(abs(Hb_stop)));
ylabel('dB')
xlabel('kHz')
title('Stopband mag. resp.')

fstart = 0; %% EDIT THIS if required!
fend = 60e3; %% EDIT THIS for bandpass or bandstop!

w_vec = 2*pi*linspace(fstart, fend, 1000);
fsample = 4*fend;
gd_ca = groupdelay(Z, P, G, w_vec, fsample);

figure(2); %% EDIT THIS: May require different number!
clf

plot(w_vec/2/pi/1e3, gd_ca/1e-3, '-+m', 'MarkerSize', 4)
hold off;

ylabel('Group delay [ms]')
xlabel('Frequency [kHz]')

------------------------------------------------------------------------------------------

% Bandpass filter
%%
%Butterworth
Amax=0.1;
Amin=40;
wc1=2*pi*50000;
wc2=2*pi*100000;
ws1=2*pi*20000;
ws2=2*pi*200000;

if wc1*wc2~=ws1*ws2
 ws1=wc1*wc2/ws2; %Öka ws1
end

wclp=wc2-wc1;
wslp=ws2-ws1;
wi2=wc1*wc2;

[Nlp, Wn]=buttord(wclp,wslp,Amax,Amin,'s');
[Zlp,Plp,Glp]=buttap(Nlp); %Normalized poles
Nlp;
epsilon=sqrt(10^(0.1*Amax)-1);
w0p=wclp*epsilon^(-1/N);
Plp=Plp*w0p;%Denormalized poles
Glp=Glp*w0p^N;
% Transformation of the lowpass filter into the bandpass filter
[Z,P]=zp2bp(Zlp,Plp,wi2);
N=2*Nlp;
Z=Z;
P=P;
G=Glp;

% Z: zeros, P: poles, G: gain
% Wc: passband edge of the lowpass filter
% Ws: stopband edge of the lowpass filter
[b, a] = zp2tf(Z, P, G); %% EDIT THIS: Use appropriate Z, P and G!
w = linspace(0, 2*ws1, 2000); %% EDIT THIS for bandpass or bandstop!
w_pass = linspace(0, wc1, 2000); %% EDIT THIS for bandpass or bandstop!
w_stop = linspace(ws1, 2*ws1, 2000); %% EDIT THIS for bandpass or bandstop!
Hb = freqs(b, a, w);
Hb_pass = freqs(b, a, w_pass);
Hb_stop = freqs(b, a, w_stop);
hFig = figure(1); %% EDIT THIS: Use different numbers!
set(hFig, 'Name', 'Butterworth') %% EDIT THIS: Approximation Type!
subplot(221)
zplane(Z/1e3, P/1e3); %% EDIT THIS: Use appropriate Z and P!
ylabel('j\omega [\times 10^3]')
xlabel('\sigma [\times 10^3]')
title('{\it s}-plane')
subplot(222)
plot(w/2/pi/1e3, 20*log10(abs(Hb)));
ylabel('dB')
xlabel('kHz')
title('Magnitude response')
subplot(223)
plot(w_pass/2/pi/1e3, 20*log10(abs(Hb_pass)));
ylabel('dB')
xlabel('kHz')
title('Passband mag. resp.')
subplot(224)
plot(w_stop/2/pi/1e3, 20*log10(abs(Hb_stop)));
ylabel('dB')
xlabel('kHz')
title('Stopband mag. resp.')
%%
%Cauer
Amax=0.1;
Amin=40;
wc1=2*pi*50000;
wc2=2*pi*100000;
ws1=2*pi*20000;
ws2=2*pi*200000;

if wc1*wc2~=ws1*ws2
 ws1=wc1*wc2/ws2; %Öka ws1
end
    
wclp=wc2-wc1;
wslp=ws2-ws1;
wi2=wc1*wc2;

% Synthesis of the lowpass filter (Cauer)
[Nlp,Wn]=ellipord(wclp,wslp,Amax,Amin,'s');
[Zlp,Plp,Glp]=ellipap(Nlp,Amax,Amin);
Nlp;
Zlp = Zlp*wclp;
Plp = Plp*wclp;
Glp = Glp*wclp^(length(Plp)-length(Zlp));
% Transformation of the lowpass filter into the bandpass filter
[Z,P]=zp2bp(Zlp,Plp,wi2);
N=2*Nlp;
Z=Z;
P=P;
G=Glp;








