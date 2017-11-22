clc;
clear;
close all;

G = tf([5 60 100],[1 26 125 100]); % Create continuous time TF as given in question
T50 = [0:0.05:5.00]; % Creating sampling times
T100 = [0:0.1:5.00];
T200 = [0:0.2:5.00];

Yplot = step(G, T50); % First plot, step response
figure(1)
plot(T50, Yplot);
title('Step Response of Continuous Time transfer function')

sys50 = c2d(G, 0.05, 'zoh') % Discretizing the transfer function for 3 different sampling times
Y50 = step(sys50, T50);
sys100 = c2d(G, 0.1, 'zoh')
Y100 = step(sys100, T100);
sys200 = c2d(G, 0.2, 'zoh')
Y200 = step(sys200, T200);

U50 = ones(100, 1); % Creating input variables
U100 = ones(50, 1);
U200 = ones(25, 1);
U50 = [0; 0; U50]; % Add zeros to make the solution full rank
U100 = [0; 0; U100];
U200 = [0; 0; U200];
Y50 = [0; 0; Y50]; % Add zeros to make the solution full rank
Y100 = [0; 0; Y100];
Y200 = [0; 0; Y200];

%%%%%%%%%%%%%%% 0.05 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y50(4:103);
M50 = [];
for k = 1:100
    m = [-1*Y50(k+2) -1*Y50(k+1) -1* Y50(k) U50(k+2) U50(k+1) U50(k)];
    M50 = [M50; m];
end
display('Estimated using pseudoinverse for T = 0.05')
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
display('Estimated using RLS for T = 0.05')
phatrls

%%%%%%%%%%%%%%% 0.1 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y100(4:53);
M100 = [];
for k = 1:50
    m = [-1*Y100(k+2) -1*Y100(k+1) -1* Y100(k) U100(k+2) U100(k+1) U100(k)];
    M100 = [M100; m];
end
display('Estimated using pseudoinverse for T = 0.1')
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
display('Estimated using RLS for T = 0.1')
phatrls

%%%%%%%%%%%%%%% 0.2 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y200(4:28);
M200 = [];
for k = 1:25
    m = [-1*Y200(k+2) -1*Y200(k+1) -1* Y200(k) U200(k+2) U200(k+1) U200(k)];
    M200 = [M200; m];
end
display('Estimated using pseudoinverse for T = 0.2')
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
display('Estimated using RLS for T = 0.2')
phatrls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Gaussian Numbers %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = normrnd(0,1,1,103);
figure(2)
plot(R)
title('Random numbers generated using Gaussian distribution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Bias Estimation %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kmat = [0.01 0.05 0.1 0.2 0.5]; 
for gain = 1:5
    Ynoisy = Y50 + kmat(gain) * R';
    phatrls = [0 0 0 0 0 0]';
    Prls = 10000000 * eye(6);
    for k = 1:100
        mrls = [-1*Ynoisy(k+2) -1*Ynoisy(k+1) -1* Ynoisy(k) U50(k+2) U50(k+1) U50(k)]';
        E = Ynoisy(k+3) - mrls' * phatrls;
        q = Prls * mrls /(1 + mrls' * Prls * mrls);
        phatrls = phatrls + q * E;
        Prls = Prls - q * mrls'  * Prls;
    end
    display('Estimated using RLS, with noise of factor')
    display(kmat(gain))
    display('for T = 0.05')
    phatrls
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PRBS %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seq50 = ltePRBS(1582,101); % Creating pseudorandom binary sequence
seq100 = ltePRBS(131, 51);
seq200 = ltePRBS(529, 26);

Y50prbs = lsim(sys50, seq50, T50);
Y100prbs = lsim(sys100, seq100, T100);
Y200prbs = lsim(sys200, seq200, T200);

%%%%%%%%%%%%%%% 0.05 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y50prbs(4:100);
M50 = [];
for k = 1:97
    m = [-1*Y50prbs(k+2) -1*Y50prbs(k+1) -1* Y50prbs(k) seq50(k+2) seq50(k+1) seq50(k)];
    M50 = [M50; m];
end
display('Estimated using pseudoinverse for T = 0.1 for PRBS')
Phat = inv(M50'*M50)*M50'*Y

%RLS
phatrls = [0 0 0 0 0 0]';
Prls = 10000000 * eye(6);
for k = 1:97
    mrls = [-1*Y50prbs(k+2) -1*Y50prbs(k+1) -1* Y50prbs(k) seq50(k+2) seq50(k+1) seq50(k)]';
    E = Y50prbs(k+3) - mrls' * phatrls;
    q = Prls * mrls /(1 + mrls' * Prls * mrls);
    phatrls = phatrls + q * E;
    Prls = Prls - q * mrls'  * Prls;
end
display('Estimated using RLS for T = 0.05 for PRBS')
phatrls

%%%%%%%%%%%%%%% 0.1 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y100prbs(4:50);
M100 = [];
for k = 1:47
    m = [-1*Y100prbs(k+2) -1*Y100prbs(k+1) -1* Y100prbs(k) seq100(k+2) seq100(k+1) seq100(k)];
    M100 = [M100; m];
end
display('Estimated using pseudoinverse for T = 0.1 for PRBS')
Phat = inv(M100'*M100)*M100'*Y

%RLS
phatrls = [0 0 0 0 0 0]';
Prls = 10000000 * eye(6);
for k = 1:47
    mrls = [-1*Y100prbs(k+2) -1*Y100prbs(k+1) -1* Y100prbs(k) seq100(k+2) seq100(k+1) seq100(k)]';
    E = Y100prbs(k+3) - mrls' * phatrls;
    q = Prls * mrls /(1 + mrls' * Prls * mrls);
    phatrls = phatrls + q * E;
    Prls = Prls - q * mrls'  * Prls;
end
display('Estimated using RLS for T = 0.1 for PRBS')
phatrls

%%%%%%%%%%%%%%% 0.2 %%%%%%%%%%%%%%%
%pseudoinverse
Y = Y200prbs(4:25);
M200 = [];
for k = 1:22
    m = [-1*Y200prbs(k+2) -1*Y200prbs(k+1) -1* Y200prbs(k) seq200(k+2) seq200(k+1) seq200(k)];
    M200 = [M200; m];
end
display('Estimated using pseudoinverse for T = 0.2 for PRBS')
Phat = inv(M200'*M200)*M200'*Y

%RLS
phatrls = [0 0 0 0 0 0]';
Prls = 10000000 * eye(6);
for k = 1:22
    mrls = [-1*Y200prbs(k+2) -1*Y200prbs(k+1) -1* Y200prbs(k) seq200(k+2) seq200(k+1) seq200(k)]';
    E = Y200prbs(k+3) - mrls' * phatrls;
    q = Prls * mrls /(1 + mrls' * Prls * mrls);
    phatrls = phatrls + q * E;
    Prls = Prls - q * mrls'  * Prls;
end
display('Estimated using RLS for T = 0.2 for PRBS')
phatrls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Determinant ratio test %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for order = [1 2 3 4 5]
   Morder = [];
   for k = order+1:103
       m = [];
       for l = 1:order
           m = [m; -Y50(k-l)];
       end
       for l = 1:order
           m = [m; U50(k-l)];
       end
       Morder  = [Morder; m'];
   end
   disp('Determinant Value for order')
   disp(order) 
   disp('is')
   disp(det(Morder' * Morder));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Simulation using order 2 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pseudoinverse
Y = Y50(3:102);
M50 = [];
for k = 1:100
    m = [-1*Y50(k+1) -1* Y50(k) U50(k+1) U50(k)];
    M50 = [M50; m];
end
display('Estimated using pseudoinverse for T = 0.05, order is 2')
Phat = inv(M50'*M50)*M50'*Y

%RLS
phatrls = [0 0 0 0]';
Prls = 10000000 * eye(4);
for k = 1:101
    mrls = [-1*Y50(k+1) -1* Y50(k) U50(k+1) U50(k)]';
    E = Y50(k+2) - mrls' * phatrls;
    q = Prls * mrls /(1 + mrls' * Prls * mrls);
    phatrls = phatrls + q * E;
    Prls = Prls - q * mrls'  * Prls;
end
display('Estimated using RLS for T = 0.05, order is 2')
phatrls
tforder2 = tf([phatrls(3) phatrls(4)],[1 phatrls(1) phatrls(2)], 0.05, 'Variable','z')
Y50order2 = lsim(tforder2, U50(2:102), T50);
figure(3)
plot(T50, Y50order2);
hold on 
plot(T50, Y50(2:102));
hold off
title('Step response of real signal & order 2 estimated signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Simulation using order 1 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pseudoinverse
Y = Y50(2:101);
M50 = [];
for k = 1:100
    m = [-1* Y50(k) U50(k)];
    M50 = [M50; m];
end
display('Estimated using pseudoinverse for T = 0.05, order is 2')
Phat = inv(M50'*M50)*M50'*Y

%RLS
phatrls = [0 0]';
Prls = 10000000 * eye(2);
for k = 1:101
    mrls = [-1* Y50(k) U50(k)]';
    E = Y50(k+1) - mrls' * phatrls;
    q = Prls * mrls /(1 + mrls' * Prls * mrls);
    phatrls = phatrls + q * E;
    Prls = Prls - q * mrls'  * Prls;
end
display('Estimated using RLS for T = 0.05, order is 2')
phatrls
tforder1 = tf([phatrls(2)],[1 phatrls(1)], 0.05, 'Variable','z')
Y50order1 = lsim(tforder1, U50(2:102), T50);
figure(4)
plot(T50, Y50order1);
hold on 
plot(T50, Y50(2:102));
hold off
title('Step response of real signal & order 1 estimated signal')



