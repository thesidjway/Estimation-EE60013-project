G = tf([5 60 100],[1 26 125 100]);
T50 = [0.05:0.05:5.00];
T100 = [0.1:0.1:5.00];
T200 = [0.2:0.2:5.00];

Yplot = step(G, T50);
figure(1)
plot(T50, Yplot);

sys50 = c2d(G, 0.05, 'zoh')
Y50 = step(sys50, T50);
sys100 = c2d(G, 0.1, 'zoh')
Y100 = step(sys100, T100);
sys200 = c2d(G, 0.2, 'zoh')
Y200 = step(sys200, T200);

U50 = ones(100, 1);
U100 = ones(50, 1);
U200 = ones(25, 1);
U50 = [0; 0; U50];
U100 = [0; 0; U100];
U200 = [0; 0; U200];
Y50 = [0; 0; 0; Y50];
Y100 = [0; 0; 0; Y100];
Y200 = [0; 0; 0; Y200];

%%%%%%%%%%%%%%% 0.05 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y50(4:103);
M50 = [];
for k = 1:100
    m = [-1*Y50(k+2) -1*Y50(k+1) -1* Y50(k) U50(k+2) U50(k+1) U50(k)];
    M50 = [M50; m];
end
Phat = inv(M50'*M50)*M50'*Y

%RLS
phatrls = [0 0 0 0 0 0]';
Prls = 10000000 * eye(6);
for k = 1:100
    mrls = [-1*Y50(k+2) -1*Y50(k+1) -1* Y50(k) U50(k+2) U50(k+1) U50(k)]';
    E = Y50(k+3) - mrls' * phatrls;
    q = Prls * mrls /(1 + mrls' * Prls * mrls);
    phatrls = phatrls + q * E;
    Prls = Prls - q * mrls'  * Prls;
end
phatrls

%%%%%%%%%%%%%%% 0.1 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y100(4:53);
M100 = [];
for k = 1:50
    m = [-1*Y100(k+2) -1*Y100(k+1) -1* Y100(k) U100(k+2) U100(k+1) U100(k)];
    M100 = [M100; m];
end
Phat = inv(M100'*M100)*M100'*Y

%RLS
phatrls = [0 0 0 0 0 0]';
Prls = 10000000 * eye(6);
for k = 1:50
    mrls = [-1*Y100(k+2) -1*Y100(k+1) -1* Y100(k) U100(k+2) U100(k+1) U100(k)]';
    E = Y100(k+3) - mrls' * phatrls;
    q = Prls * mrls /(1 + mrls' * Prls * mrls);
    phatrls = phatrls + q * E;
    Prls = Prls - q * mrls'  * Prls;
end
phatrls

%%%%%%%%%%%%%%% 0.2 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y200(4:28);
M200 = [];
for k = 1:25
    m = [-1*Y200(k+2) -1*Y200(k+1) -1* Y200(k) U200(k+2) U200(k+1) U200(k)];
    M200 = [M200; m];
end
Phat = inv(M200'*M200)*M200'*Y

%RLS
phatrls = [0 0 0 0 0 0]';
Prls = 10000000 * eye(6);
for k = 1:25
    mrls = [-1*Y200(k+2) -1*Y200(k+1) -1* Y200(k) U200(k+2) U200(k+1) U200(k)]';
    E = Y200(k+3) - mrls' * phatrls;
    q = Prls * mrls /(1 + mrls' * Prls * mrls);
    phatrls = phatrls + q * E;
    Prls = Prls - q * mrls'  * Prls;
end
phatrls

%%%%%%%%%%%%%%% Gaussian Numbers %%%%%%%%%%%%%%%
R = normrnd(0,1,1,103);
figure(2)
plot(R)

%%%%%%%%%%%%%%% Bias Estimation %%%%%%%%%%%%%%%
kmat = [0.01 0.05 0.1 0.2 0.5];
for gain = 1:5
    Ynoisy = Y50 + kmat(gain) * R;
    phatrls = [0 0 0 0 0 0]';
    Prls = 10000000 * eye(6);
    for k = 1:100
        mrls = [-1*Ynoisy(k+2) -1*Ynoisy(k+1) -1* Ynoisy(k) U50(k+2) U50(k+1) U50(k)]';
        E = Ynoisy(k+3) - mrls' * phatrls;
        q = Prls * mrls /(1 + mrls' * Prls * mrls);
        phatrls = phatrls + q * E;
        Prls = Prls - q * mrls'  * Prls;
    end
    phatrls
end

%%%%%%%%%%%%%%% PRBS %%%%%%%%%%%%%%%
seq = ltePRBS(162,103);
