H = tf([1 0.1 7.5],[1 0.12 9 0 0]);
w = logspace(-1,1,50);
bode(H,w)