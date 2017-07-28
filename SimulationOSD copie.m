%--- SIMULATION : Ordering Statistics Decoding : Soft-Decison Decoding ---%
%---------------- Application au code de Hamming (7,4,3) -----------------%

clear all;
close all;
seed = 1000; 
rand('state',seed);  % initialisation de la fonction rand()
randn('state',seed); % initialisation de la fonction randn() 

%---------------------- Les parametres du code ---------------------------

n = 7;                 % longueur du code
k = 4;                 % longueur du message
R = k/n;               % rendement du code
dmin = 3;              % distance minimale de Hamming du code
A_dmin = 7;             % Nombre de mots de poids minimale
                       % Pour le code de Hamming (7,4,3)
                       
G = [1 0 0 0 1 0 1     % Matrice generatrice
     0 1 0 0 1 1 1
     0 0 1 0 1 1 0
     0 0 0 1 0 1 1];
 
H=[1 1 0 1 1 0 0 ; 1 0 1 1 0 1 0 ; 0 1 1 1 0 0 1];

%--------------------- Les parametres de simulation -----------------------

SNRdb = [0:12];                            % SNR en dB
nb_mots = 1e4;                             % Nombre de mots a simuler

% -------------------- La boucle de simulation ----------------------------

         

Teb_osd0 = zeros(1,length(SNRdb));        % Initialisation du TEB

%while nb_mots < 10000 && Teb_osd0 < 1e-6

for kk = 1:length(SNRdb)
    
    % Affichage de l'etape de simulation
    fprintf('Simulation pour SNR = %.2f dB\n', SNRdb(kk));
   
    % Rapport signal sur bruit 
    
    SNRlin = 10.^(SNRdb(kk)./10);       % SNR en echelle lineaire(normalisé)
    sigmaw = 1/(2*R*SNRlin);            % variance du bruit
    
    % Generation des mots de codes
    
    for mots = 1:nb_mots
    
        % Generation et modulation de c
        m = rand(1,k)>0.5;
        c = mod(m*G,2);
        %c = zeros(1,n);
        x = 2*c-1; 
        
        % Passage dans le canal AWGN
        w = sqrt(sigmaw)*randn(1,n);  % generation du Bruit Blanc Gaussien
        r = x+w;                      % l'observation
     
        % Decodeur du mot recu : Decodage OSD
        c_decode0 = DecodeurOSD(r);
      
        % Comptage des erreurs osd
        nb_erreur_osd0 = sum(m ~= mess_decod_hdd);
		Teb_osd0(kk) = Teb_osd0(kk) + nb_erreur_osd0;  
        
    end  
   
end

 Teb_osd0 = Teb_osd0/(k*nb_mots);
 
% ----------- Simulations : Courbe de performance---------------

% Modulation MDP-2 

% Probabilité d'erreur binaire apres decodage ferme
mdp2_ferme = 0.5*erfc(sqrt(10.^(SNRdb/10))); 

% Probabilité d'erreur binaire apres decodage souple
mdp2_souple = 0.5*(n^(-1))*dmin*A_dmin*erfc(sqrt(dmin*R*sqrt(10.^(SNRdb/10))));

figure(1);
% Courbe et compraison :Evolution du TEB apres decodage
semilogy(SNRdb,mdp2_ferme,'--r','markersize',10,'linewidth',0.5);     % Teb transmission non-code osd,ferme
hold on

semilogy(SNRdb,mdp2_souple,'m-.','markersize',10,'linewidth',0.5);     % Teb transmission non-code osd,souple
hold on

semilogy(SNRdb,Teb_osd0,'b-*','markersize',08,'linewidth',0.5); % OSD d'ordre 0

axis([0 12 10^-7 10^0]) 

legend({'Transmission non codée ferme','Transmission non codée souple',...
         'Simulation - OSD d''ordre 0'},'Location','best','FontSize',15,...
         'FontWeight','bold','color','w');

title('Evolution du TEB en fonction de Eb/No',...
    'FontSize',15,'FontWeight','bold','color','k');
xlabel('Eb/No (dB)','FontSize',15,'FontWeight',...
    'bold','color','k');
ylabel('TEB','FontSize',15,'FontWeight','bold','color','k');
grid on 
