%% Gruppe nm AA, Abgabe 28.09.2023
% Gruppenitglieder:
% Rêzan Berhan Işik
% Junui Lee
% Sokratis Panagiotidis
% Valentin Benjamin Jimenez Anders
clc; clear; close all;                                                      % Knock everything out.

%% Initialisierung
n = 5;                                                                      % Anzahl der Elemente in 1.
L = 1;                                                                      % Laenge des Balkens in m.
h = L/n;
l = 0:n-1;
n_roof = 5;                                                                 % Anzahl der zusaetzlichen Auswertungspunkte je Element in 1.
n_tilde = 7;                                                                % Ordnung der Quadratur in 1.
n_P = 50;                                                                   % Anzahl der Zeitschritte in 1.
x_tilde = linspace(0,1,n_tilde+1);

%% Aufgabe 14, initial conditions.
% beta = 1/4; gamma = 1/2; eta = 1;
beta = 1/4; gamma = 1/2; eta = .1;
% beta = 1/4+.5; gamma = 1/2+1; eta = .1;                                     % Konvergenz gegen Null.
% beta = 1/4+.001; gamma = 1/2+.2; eta = 1;                                   % Konvergenz gegen Unendlich.

mu = @(x) x.^0;                                                             % Laengenspeifische Masse in kg/m.  
E = @(x) x.^0;                                                              % Elastizitätsmodul in kg/m.
I = @(x) x.^0;                                                              % Flaeschentraegheitsmoment in m^4.
q = @(x) x.^0;                                                              % Streckenlast in N/m.

B = [0 1 0;                                                                 % Auslenkung linkes Ende in m.
     0 2 0;                                                                 % Anstieg linkes Ende in 1.
     n 3 0;                                                                 % Moment rechtes Ende in Nm.
     n 4 0];                                                                % Querkraft rechtes Ende in N.#

[ldreid,idreid,jdreid,zweilidreid,zweiljdreid,lzweid,izweid,zweilizweid] = get_indizes(n);

Mbar = getMbar(mu, n, n_tilde, x_tilde);
M = getM(Mbar,zweilidreid,zweiljdreid);

Sbar = getSbar(E, I, n, n_tilde, x_tilde);
S = getS(Sbar, zweilidreid, zweiljdreid);

qbar = getqbar(q, n, n_tilde, x_tilde);

vq = getvq(qbar, zweilizweid);
C = getC(B,n);
vn = getvn(B, n);
vd = getvd(B);

Se = getSe(S, C);
Me = getMe(M);
ve = getve(vn, vd, vq);

ae = Se\ve;          

Newmark(q, n, n_tilde, x_tilde, n_P, L, n_roof, ae, vn, vd, Se, beta, gamma, eta, Me, l, M, S, C)

Konvergenz(q, L, n, n_tilde, x_tilde, E, I, n_roof, l)

function [ldreid,idreid,jdreid,zweilidreid,zweiljdreid,lzweid,izweid,zweilizweid] = get_indizes(n)
    %% 3D-Arrays
    [i,j,l]=ndgrid(0:3,0:3,0:n-1);
    
    ldreid=l;
    idreid=i;
    jdreid=j;
    
    zweilidreid=2.*l+i;
    zweiljdreid=2.*l+j;
    
    %% 2D-Arrays
    i=i(:,1,:);
    l=l(:,1,:);
    
    lzweid=l;
    izweid=i;
    zweilizweid=2.*l+i;
    
    idreid = idreid+1;
    jdreid = jdreid+1;
    ldreid = ldreid+1;
    
    izweid = izweid+1;
    lzweid = lzweid + 1;
end

function A = getA(Abar,zweilidreid,zweiljdreid)

    Abar=reshape(Abar,[],1);
    i = reshape(zweilidreid,[],1)+1;
    j = reshape(zweiljdreid,[],1)+1;
    A = sparse(i,j,Abar);

end

function A_sparse = getAbar(count, ~, n_tilde, x_tilde)
           [~, Phi_4Di, Phi_4Dj, ~] = getphi(n_tilde, count, x_tilde);

            [~, ~, ~, ~, ~, s4D] = getstencil(n_tilde, count);
            A_produkt = Phi_4Di.*Phi_4Dj;
            
            A_sparse = sum(A_produkt.*s4D, 4);
end

function [A_sparse] = getAbarNewmark(n, n_roof)
    kUS_roof = (0:n_roof);
    j = 0:3;
    l = 0:n-1;
    
    [~, h3D] = geth(n);
    [~, ~, delta_j] = getexp(n);
    
    [J, K, L] = meshgrid(j, kUS_roof,l);
    k = kUS_roof/n_roof;

    phi0 = @(x) 1 - 3*x.^2 + 2*x.^3;
    phi1 = @(x) x - 2*x.^2 + x.^3;
    phi2 = @(x) 3*x.^2 - 2*x.^3;
    phi3 = @(x) - x.^2 + x.^3;

    phi_j = [phi0(k); phi1(k); phi2(k); phi3(k)];


    Faktor = h3D.^repmat(delta_j(1,:), [1 1 n]);                            
    A = Faktor(:,:,1).*phi_j';                                              

    K = n_roof*L + K;
    J = 2*L + J;
    
    A_reshape =reshape(repmat(A, [1 1 n]),[],1);

    Khut = reshape(K,[],1)+1;
    J_bla = reshape(J,[],1)+1;
    A_roof = sparse(Khut,J_bla,A_reshape);

% Korrektur: Überlappende Werte, welche addiert wurden, werden
% wiederhergestellt.
    Index1 = A_roof(1,1);
    Dividiere = A_roof == Index1*2;
    A_roof(Dividiere)./2;
    A_sparse = A_roof .* (~Dividiere) + (A_roof / 2) .* Dividiere;
end

function C = getC(B,n) % Gruppe nm AA

    % Erstellen von E1 und E2 als sparse-Matrizen
    % sparse( Wert in Zeile, Wert in Spalte, Wert des Wertes (lol), Anzahl
    % Zeilen, Anzahl Spalten).
    N = 2*n+2;
    K1 = B(B(:,2)==1,1);
    K2 = B(B(:,2)==2,1);

    E1 = sparse(2*K1 + 1, 1:length(K1), 1, N, length(K1)); % Ansonsten 4x1 Vector, wenn size(B,1) nicht dupliziert wird.
    E2 = sparse(2*K2 + 2, 1:length(K2), 1, N, length(K2));

    clc
    C=[E1, E2];

end

function [ddphi_3Di, ddphi_4Di, ddphi_4Dj] = getddphi(n, n_tilde, x_tilde)
    
    % x_tilde = linspace(0,1,n_tilde+1);

    ddphi0 = @(x) -6 + 12*x;
    ddphi1 = @(x) -4 + 6*x;
    ddphi2 = @(x) 6 - 12*x;
    ddphi3 = @(x) -2 + 6*x;

    ddphi = [ddphi0(x_tilde); ddphi1(x_tilde); ddphi2(x_tilde); ddphi3(x_tilde)];
    [~,idreid,jdreid,~,~,~,izweid,~] = get_indizes(n);

    %% I : 3D & 4D
    ddphi_3Di = reshape(ddphi(izweid,:), 4, n, n_tilde+1);
    ddphi_4Di = reshape(ddphi(idreid,:), 4 , 4, n, n_tilde+1);

    %% J : 4D
    ddphi_4Dj = reshape(ddphi(jdreid,:), 4 , 4, n, n_tilde+1);

end

function [getexp_2D, getexp_3D, delta_j] = getexp(n)
    [~,idreid,jdreid,~,~,~,~,~] = get_indizes(n);
    % 2D : deltai1 + deltai3    
    % Wenn Indizes der Zeilen gleich 1 oder 3
    deltai1 = idreid(idreid(:,1)==1,1) +1;
    deltai3 = idreid(idreid(:,1)==3,1) +1;

    delta_i1 = zeros(4,4);
    delta_i1(deltai1,:) = 1;

    delta_i3 = zeros(4,4);
    delta_i3(deltai3,:) = 1;

    delta_i = delta_i1 + delta_i3;

    getexp_2D = delta_i(:,1);

    deltaj1 = jdreid(jdreid(:,1)==1,1) +1;
    deltaj3 = jdreid(jdreid(:,3)==3,3) +1;

    delta_j1 = zeros(4,4);
    delta_j1(:,deltaj1) = 1;

    delta_j3 = zeros(4,4);
    delta_j3(:,deltaj3) = 1;

    delta_j = delta_j1 + delta_j3;

    getexp_3D = repmat((delta_i + delta_j), [1 1 n]);

   
end

function [hD, h3D] = geth(n)
    x = linspace(0, 1, n+1);
    hD = x(2:end)-x(1:end-1); % 2D Array 1xn_tilde
    
    h3D = reshape(hD, 1,1,[]);
end

function M = getM(Mbar,zweilidreid,zweiljdreid)

    Mbar=reshape(Mbar,[],1);
    i = reshape(zweilidreid,[],1)+1; % Funktioniert nicht ohne das +1. Wir wissen nicht wieso.
    j = reshape(zweiljdreid,[],1)+1;
    M = sparse(i,j,Mbar);

end

function [Mbar] = getMbar(mu, n, n_tilde, x_tilde)
    [~, h3D] = geth(n);
    [~, getexp_3D] = getexp(n);

    h3Dexp = h3D .^ (getexp_3D+1);

    [~, T4] = getTinv(n_tilde, n);
    mu = mu(T4);                                                            % mu = x.^0 = 1. Somit T4 kein Einfluss, und mu muss ja 1 sei.
    [~, Phi_4Di, Phi_4Dj] = getphi(n_tilde, n, x_tilde);
    Produkt = mu.*Phi_4Di.*Phi_4Dj;                                         % Aussage von Marian: Phi_4Di und Phi_4Dj sind korrekt.
    [~, ~, ~, ~, ~, s4D] = getstencil(n_tilde, n);                          % Die ganzen anderen Größen vom Stencil sind unbedeutend.
                                                                            % s4D ist lediglich s (2D) in die 4. Dimension verschoben.
    Summe = sum(Produkt.*s4D, 4);                                           % Produkt.*s4D sind alle korrekt.
    Mbar = h3Dexp.*Summe; 
end

function [Me] = getMe(M)
    Me = [M, zeros(size(M,1),2)];
    Me = [Me; zeros(2,size(Me,2))];
end

function [Phi_3Di, Phi_4Di, Phi_4Dj, phi] = getphi(n_tilde, n, k)

    [~,idreid,jdreid,~,~,~,izweid,~] = get_indizes(n);
    phi0 = @(x) 1 - 3*x.^2 + 2*x.^3;
    phi1 = @(x) x - 2*x.^2 + x.^3;
    phi2 = @(x) 3*x.^2 - 2*x.^3;
    phi3 = @(x) - x.^2 + x.^3;

    %% Attention. Hier ist nachholbedarf zur Universalität.
    phi = [phi0(k); phi1(k); phi2(k); phi3(k)];

%% I : 3D & 4D
    Phi_3Di = reshape(phi(izweid,:), 4, [], n_tilde+1); % korrekt
    Phi_4Di = reshape(phi(idreid,:), 4 , 4, [], n_tilde+1);  % korrekt
%% J : 4D
    Phi_4Dj = reshape(phi(jdreid,:), 4 , 4, [], n_tilde+1);% korrekt
end

function [w] = getplot(ae, A_sparse, n)
    N = 2*n+2;
    a = ae(1:N);
    w = A_sparse*a;
end

function [qbar] = getqbar(q, n, n_tilde, x_tilde)
    [T3, ~] = getTinv(n_tilde, n);                                         % T4 als 1x1xlxk_tilde bzw.  1x1xnxn_tilde+1
    [h, ~] = geth(n);
    [getexp_2D, ~] = getexp(n);
    [Phi_3Di, ~, ~] = getphi(n_tilde, n, x_tilde);
    [~, ~, ~, ~, s3D, ~] = getstencil(n_tilde, n);
    
    q = q(T3);
    
    Produkt_q = q.*Phi_3Di;
    
    Summe_q = sum(Produkt_q.*s3D, 3);                                       % 4xn Matrix.
    Faktor = (h.^(getexp_2D+1));
    qbar = Faktor.*Summe_q;
end

function [S] = getS(Sbar, zweilidreid, zweiljdreid)

    Sbar=reshape(Sbar,[],1);
    i = reshape(zweilidreid,[],1)+1; % Funktioniert nicht ohne das +1. Wir wissen nicht wieso.
    j = reshape(zweiljdreid,[],1)+1;
    S = sparse(i,j,Sbar);
end

function Sbar = getSbar(E, I, n, n_tilde, x_tilde)   

    [~, h3D] = geth(n);
    [~, getexp_3D] = getexp(n);

    Sbar_h3Dexp = h3D .^ (getexp_3D(:, :, 1)-3);

    [~, ddphi_4Di, ddphi_4Dj] = getddphi(n, n_tilde, x_tilde);
    [~, T_inv4D] = getTinv(n_tilde, n);                                     % T_inv4D als 1x1xlxk_tilde bzw. 1x1xnxn_tilde+1
    E = E(T_inv4D);
    I = I(T_inv4D);

    Produkt = E.*I.*ddphi_4Di.*ddphi_4Dj;

    [~, ~, ~, ~, ~, s4D] = getstencil(n_tilde, n);

    Summe = sum(Produkt.*s4D, 4);
    Sbar = Sbar_h3Dexp.*Summe;
end

function [Se] = getSe(S, C)
    Se1 = [S,C];
    Se2 = [transpose(C),zeros(2,2)];
    Se = vertcat(Se1, Se2);
end

function [V, V_plus1, s, s2D, s3D, s4D] = getstencil(n_tilde, n)
    K_tilde = (0:1:n_tilde)';
    x_tilde = linspace(0,1,n_tilde+1)';

    % V = zeros
    i_tilde = K_tilde;
    k_tilde = i_tilde;
%% Vandermondematrix
    V_tilde = repmat(x_tilde, [1 length(x_tilde)]);
    V =  V_tilde .^ repmat(i_tilde', [length(x_tilde) 1]);
    % V_vergleich = fliplr(vander(x_tilde))
    % if (isequal(V, V_vergleich) == 1)
    %     disp("ERFOLG!")
    % else
    %     disp("FAIL!")
    % end
    V_plus1 = 1./(i_tilde+1);

    s = (V'\V_plus1)';
    % reshape(phi, 1, length(i), n_tilde+1)
    s2D = reshape(s, [], 1, n_tilde+1);
    s3D = reshape(s2D, 1, 1, n_tilde+1);
    s4D = reshape(s2D, [], 1, 1, n_tilde+1);
    
end

function [T_inv3D,T_inv4D] = getTinv(n_tilde, n)

    [ldreid,idreid,jdreid,zweilidreid,zweiljdreid,lzweid,izweid,zweilizweid] = get_indizes(n);
    x_tilde = linspace(0,1,n_tilde+1);
    l = 1; 
    [hD, h3D] = geth(n);

    xl = linspace(0, l, n+1)';
    
    T_inv2D = [hD'*x_tilde + xl(1:end-1)];

    T_inv3D = reshape(T_inv2D, 1, [], n_tilde+1); % qbar
    T_inv4D = reshape(T_inv2D(ldreid,:), 4 , 4, [], n_tilde+1); %Mbar
end

function vd = getvd(B)
    a = sparse(1, 1, B(1,3), 1, 1);
    b = sparse(1, 1, B(2,3), 1, 1);
    vd = [a, b]';
end

function [ve] = getve(vn, vd, vq)
    ve = [vq+vn;vd];
end

function vn = getvn(B, n)
    % Erstellen von E1 und E2 als sparse-Matrizen
    N = 2*n+2;

    K3 = B(B(:,2)==3,1);
    K4 = B(B(:,2)==4,1);

    c3 = B(B(:,2) == 3,3);
    c4 = B(B(:,2) == 4,3);
    % 
    % E3 = sparse(N, 1, 1, N, 1); % Ansonsten 4x1 Vector, wenn size(B,1) nicht dupliziert wird.
    % E4 = sparse(N-1, 1, 1, N, 1);
    E3 = sparse(2*K3 + 2, 1:length(K3), 1, N, length(K3)); % Ansonsten 4x1 Vector, wenn size(B,1) nicht dupliziert wird.
    E4 = sparse(2*K4 + 1, 1:length(K4), 1, N, length(K4));
    % c3 = B(3,3);
    % c4 = B(4,3);
    vn=[E3*c3+E4*c4];
end

function vq = getvq(qbar,zweilizweid)

    qbar=reshape(qbar,[],1);
    i = reshape(zweilizweid,[],1)+1;
    vq=sparse(i,1,qbar);

end

function [] = Konvergenz(q, L, n, n_tilde, x_tilde, E, I, n_roof, l)
    z13 = [];

    %% Konvergenz, statisch.
        for count = 1:1000
            B = [0 1 0;                                                                                         % Auslenkung linkes Ende in m.
                 0 2 0;                                                                                         % Anstieg linkes Ende in 1.
                 count 3 0;                                                                                     % Moment rechtes Ende in Nm.
                 count 4 0];
            N = 2*count+2;
            [~,~,~,zweilidreid,zweiljdreid,~,~,zweilizweid] = get_indizes(count);
            [~, ~] = geth(count);
            qbar = getqbar(q, count, n_tilde, x_tilde);
            vn = getvn(B, count);
            vq = getvq(qbar, zweilizweid);
            vd = getvd(B);
            ve = getve(vn, vd, vq);
            
            Sbar = getSbar(E, I, count, n_tilde, x_tilde);
            S = getS(Sbar, zweilidreid, zweiljdreid);
            C = getC(B,count);
            Se = getSe(S, C);
            
            ae = Se\ve;
            alpha = ae(1:N); % Länge: N = 2*n + 2
            A_sparse = getAbar(count, n_roof, n_tilde, x_tilde);
            A = getA(A_sparse,zweilidreid,zweiljdreid);
            % w = getplot(ae, A_sparse, count);
            
            x = linspace(0, 1, count+1);                                                                      % Es stand count +1 (habe es entfernt, weil sonst die Vektorlängen nicht mehr übereinstimmen)
            % xL = 1;
            w = transpose((q(x)./(E(x).*I(x))).*((x.^4)/24 - (L*x.^3)/6 + (L^2*x.^2)/4));
            w_prime = transpose((q(x)./(E(x).*I(x))).*((x.^3)/6 - (L*x.^2)/2 + (L^2*x)/2));
            omega = zeros(N,1);
            omega(1:2:end) = w;
            omega(2:2:end) = w_prime;
            
            % Definition von A (Das gleiche wie für M bloß immer mi = 1)
            % A = getA(A_sparse, zweilidreid, zweiljdreid);                   % Probleme mit der Größe der Vektoren
            
            % w = getplot(n,  ae);
            
            %Definition der error_curve
            error_numerator = sqrt(transpose(omega-alpha)*A*(omega-alpha));
            error_denominator = sqrt(transpose(omega)*A*omega);
            error_curve = error_numerator/error_denominator;
            z13(count,1) = error_curve;
        end

%%  Fehler, arithmetisch
    figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(1:count,z13, 'r-', 'LineWidth', 2);
    ylabel('$log$ $error_{L^2}$ $in$ $1$' , 'Interpreter', 'latex')
    xlabel('$log$ $n$ $in$ $1$' , 'Interpreter', 'latex')
    legend('$Relativer$ $Fehler$ $logaritmiert$', 'Interpreter', 'latex')
    set(gca, 'FontSize', 15)

%%  Fehler, doppelt logarithmisch    
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(1:count,z13, 'g-', 'LineWidth', 2);
    ylabel('$error_{L^2}$ $in$ $1$' , 'Interpreter', 'latex')
    xlabel('$n$ $in$ $1$' , 'Interpreter', 'latex')
    legend('$Relativer$ $Fehler$', 'Interpreter', 'latex')
    set(gca, 'FontSize', 15)
end

function [] = Newmark(q, n, n_tilde, x_tilde, n_P, L, n_roof, ae, vn, vd, Se, beta, gamma, eta, Me, l, M, S, C)
    x_axis = linspace(0,L,n*n_roof+1);
    qbar = getqbar(q, n, n_tilde, x_tilde);
    [ldreid,idreid,jdreid,zweilidreid,zweiljdreid,lzweid,izweid,zweilizweid] = get_indizes(n);
    vq = getvq(qbar, zweilizweid);
    %ve = getve(vn, vd, vq);                                                % Wird nicht benutzt
    
    q = @(x) x.*0;
    qbar = getqbar(q, n, n_tilde, x_tilde);
    vq = getvq(qbar, zweilizweid);
    ve = getve(vn, vd, vq);


    t=linspace(0,n_P*1,n_P);
    z14 = [];
    ae_prime = 0; ae_double_prime = 0; %c1 = 0; c2 = 0;
    N = 2*n + 2;
    w_history = zeros(length(x_axis), n_P);
    figure; % Create a new figure before the loop

% Initialize an empty plot
h = plot(x_axis, NaN(size(x_axis)), 'k-', 'LineWidth', 2);
title({sprintf('$Bending$ $of$ $a$ $beam$ $for$ $n=%d$', n)}, 'Interpreter', 'latex');
ylabel('$\mathrm{w}$ in $m$', 'Interpreter', 'latex');
xlabel('$\mathrm{Length}$ $\mathrm{of}$ $\mathrm{beam}$ $x$ $[\mathrm{m}]$', 'Interpreter', 'latex');
legend('$Newmark$', 'Interpreter', 'latex');
axis([0 1 -.25 0.25]);

for p = 1:n_P
    % ae_star
    ae_star = ae + ae_prime * eta + (0.5 - beta) * ae_double_prime * (eta)^2;
    ae_star_prime = ae_prime + (1 - gamma) * ae_double_prime * eta;
    % ae
    ae_double_prime = (Me + Se * beta * eta^2) \ (ve - Se * ae_star);
    ae_prime = ae_star_prime + gamma * ae_double_prime * eta;
    ae = ae_star + beta * ae_double_prime * eta^2;

    A_sparse = getAbarNewmark(n, n_roof);
    w = getplot(ae, A_sparse, n);
    nu = -1 * ae(N + 1:end); % For 2 RB. Adjust as needed for the general case

    E_ges(p, 1) = (1/2) * ae_prime(1:2 * n + 2, :)' * M * ae_prime(1:2 * n + 2, :) + ((1/2) * S * ae(1:2 * n + 2, :) - vq - C * nu - vn)' * ae(1:2 * n + 2, :);

    % Store the data for animation
    w_history(:, p) = w;

    % Update the plot data
    set(h, 'YData', w_history(:, p));
    
    % Update the plot without opening a new figure
    drawnow;
    
    % Add a pause to control animation speed
    pause(0.01); % Adjust the duration as needed
end
        close all;
        
%%  Gesamtenergie
        figure('units','normalized','outerposition',[0 0 1 1]);
        plot(t,E_ges, 'b-', 'LineWidth', 2);
        ymax = max(E_ges);
        ylim([0 ymax*1.1])
        title('$Convergence$ $analysis$', 'Interpreter', 'latex')
        legend('$total$ $energy$', 'Interpreter', 'latex');
        xlabel('$\mathrm{time}\ t\ [\mathrm{s}]$', 'Interpreter', 'latex');
        ylabel('$E_{\mathrm{ges}}$ [J]', 'Interpreter', 'latex');
       % Bezeichnung in Legende ändern legend('$E_{\mathrm{ges}}$ [J]', 'Interpreter', 'latex')
        set(gca, 'FontSize', 15)
end
