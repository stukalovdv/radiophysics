while n<250000
    fdtd_classic
end
[d,fs]=audioread('mp3_last.mp3');
sound (d,fs);