
p = audioplayer(f1,Fs);
play(p,[1 (get(p,'SampleRate')*8)]);
pause(4);
p = audioplayer(f2,Fs);
play(p,[1 (get(p,'SampleRate')*8)]);
pause(4);
p = audioplayer(f3,Fs);
play(p,[1 (get(p,'SampleRate')*8)]);
pause(4);
p = audioplayer(f4,Fs);
play(p,[1 (get(p,'SampleRate')*8)]);

% f1 was just deep, oscillating sounds.
% f2 sounded like a muffled version of the original 
% f3 was similar to f2, just more high pitched. Sounded like it was coming
% through a radio 
% f4 was similar to f3, just softer and higher pitched. 

