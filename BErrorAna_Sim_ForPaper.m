% 5 rooms example
clc
clear all


% load('Kinv_opt_2','Kinv_opt_2')
% load('Kinv_Pol_2','Kinv_Pol_2')
% R = [-10000 -10000 -10000 -10000 0 -10000
%     -10000 -10000 -10000 0 -10000 100
%     -10000 -10000 -10000 0 -10000 -10000
%     -10000 0 0 -10000 0 -10000
%     0 -10000 -10000 0 -10000 100
%     -10000 0 -10000 -10000 0 100];





R = [0 -100 -100 -100 -5 -100
    -100 0 -100 -5 -100 100
    -100 -100 0 -5 -100 -100
    -100 -5 -5 0 -100 -100
    -5 -100 -100 -5 0 100
    -100 -5 -100 -100 -5 100];


T = 1000000; % Number of iterations for learning

t_step = 100;
TVEC2 = t_step:t_step:T;
T_small = length(TVEC2);


xt = 0; % state 0;
xtdash = xt + 1; % Because matlab does not take 0 as a valid index.


beta = 0.8; % Discount factor


alpha_a = 0.6; % Polyak gain factor
alpha_b = 0.52; % Polyak gain factor
alpha_c = 0.85; % Polyak gain factor

noise = 0.2;

drho = 1;


% AlenVec = [1 2 1 3 3 3]; % Number of actions in each state
% 
% 
% NumbertoSAPair = [0 4; 1 3; 1 5; 2 3; 3 1; 3 2; 3 4; 4 0; 4 3; 4 5; 5 1; 5 4; 5 5];


AlenVec = [2 3 2 4 4 3]; % Number of actions in each state
NumbertoSAPair = [0 0; 0 4; 1 1; 1 3; 1 5; 2 2; 2 3; 3 1; 3 2; 3 3; 3 4; 4 0; 4 3; 4 4; 4 5; 5 1; 5 4; 5 5];



% Location of the state-action pair indicates the unique number which
% represents it (e.g. [0 0] = 1)


SAPairtoNumber = zeros(6,6); % Basically doing the same thing what I just mentioned; S-A -> unique identity
for ii = 1:length(NumbertoSAPair)
        SAPairtoNumber(NumbertoSAPair(ii,1)+1,NumbertoSAPair(ii,2)+1) = ii;
end


StaActLen = sum(AlenVec); % Total number of state action pairs

ones_StaActLen = ones(StaActLen,1);


psixt = zeros(StaActLen,1);
for jj = 1:6        
    if(SAPairtoNumber(xtdash,jj)>0)
        psixt(SAPairtoNumber(xtdash,jj)) = 1;        
    end
end

for count = 1:1
    
    
    % MATRIX GAINS:::::

Kinv_SNR2a = 10*eye(StaActLen,StaActLen); % Polyak averaging gain 1
Kinv_SNR2b = 10*eye(StaActLen,StaActLen); % Polyak averaging gain 2
Kinv_SNR2c = 10*eye(StaActLen,StaActLen); % Polyak averaging gain 3


Kinv_BVR = 10*eye(StaActLen,StaActLen); % Gain used by BVR and DC


Kinv_DQ_a = 10*eye(StaActLen,StaActLen);
Kinv_DQ_b = 10*eye(StaActLen,StaActLen);


    % Theta Initializations:::::

Qttheta_DQ_a = rand(StaActLen,1);  % Q value at time t; DQ
Qttheta_DQ_b = rand(StaActLen,1);  % Q value at time t; DQ


Qttheta_SNR2a = rand(StaActLen,1);  % Q value at time t; SNR2
Qttheta_SNR2b = rand(StaActLen,1);  % Q value at time t; SNR2
Qttheta_SNR2c = rand(StaActLen,1);  % Q value at time t; SNR2


Qttheta_Pol_a = rand(StaActLen,1);  % Q value at time t; SNR2
Qttheta_Pol = rand(StaActLen,1);  % Q value at time t; SNR2


Qttheta_BVRK = rand(StaActLen,1);  % Q value at time t; BVR gain


Qttheta_Speedy = rand(StaActLen,1);  % Q value at time t; Speedy
PrevQttheta_Speedy = Qttheta_Speedy;


Qttheta_Wat = rand(StaActLen,1);

     % Q functions:::::::

    
Qfun_DQ = zeros(6,6);

Qfun_SNR2a = zeros(6,6);
Qfun_SNR2b = zeros(6,6);
Qfun_SNR2c = zeros(6,6);

Qfun_Pol = zeros(6,6);

Qfun_BVRK = zeros(6,6);
Qfun_Speedy = zeros(6,6);

Qfun_Wat = zeros(6,6);



% Evolution of BE

Qttheta_DQ_BEevol = zeros(T_small,StaActLen);

Qttheta_SNR2a_BEevol = zeros(T_small,StaActLen);
Qttheta_SNR2b_BEevol = zeros(T_small,StaActLen);
Qttheta_SNR2c_BEevol = zeros(T_small,StaActLen);

Qttheta_Pol_BEevol = zeros(T_small,StaActLen);

Qttheta_BVR_BEevol = zeros(T_small,StaActLen);
Qttheta_Speedy_BEevol = zeros(T_small,StaActLen);


Qttheta_Wat_BEevol = zeros(T_small,StaActLen);



% Evolution of Average BE


Qttheta_DQ_BEAvgevol = zeros(T_small,StaActLen);

Qttheta_SNR2a_BEAvgevol = zeros(T_small,StaActLen);
Qttheta_SNR2b_BEAvgevol = zeros(T_small,StaActLen);
Qttheta_SNR2c_BEAvgevol = zeros(T_small,StaActLen);


Qttheta_Pol_BEAvgevol = zeros(T_small,StaActLen);


Qttheta_BVR_BEAvgevol = zeros(T_small,StaActLen);
Qttheta_Speedy_BEAvgevol = zeros(T_small,StaActLen);


Qttheta_Wat_BEAvgevol = zeros(T_small,StaActLen);



%%%%%%%%% BE FOR Q %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Qttheta_DQ_BEevolQ = zeros(T_small,StaActLen);

Qttheta_SNR2a_BEevolQ = zeros(T_small,StaActLen);
Qttheta_SNR2b_BEevolQ = zeros(T_small,StaActLen);
Qttheta_SNR2c_BEevolQ = zeros(T_small,StaActLen);

Qttheta_Pol_BEevolQ = zeros(T_small,StaActLen);

Qttheta_BVR_BEevolQ = zeros(T_small,StaActLen);
Qttheta_Speedy_BEevolQ = zeros(T_small,StaActLen);


Qttheta_Wat_BEevolQ = zeros(T_small,StaActLen);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:T
    
%     t
    
    % POLICY::::::::::::::::
    
%     if(rand>0.0)
    at = ChooseAction(xt);
%     else
%         QQ = (Qttheta_DQ_a + Qttheta_DQ_b)./2;
%         QQ = QQ.*psixt;
% %         QQf = zeros(6,6);
% %         for jj = 1:StaActLen
% %             QQf(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = QQ(jj);
% %         end
%         [~,samax] = max(QQ);
%         at = NumbertoSAPair(samax,2);
%             
%     end



%     xtp1 = at + NOISE ; % Next state

    
      if(rand < noise)
          xtp1 = ChooseAction(xt);
      else
          xtp1 = at;
      end
      
    
    
    
    % STEPSIZE ::::::::::::::::::::::
    
    gammat = 1./(t+1); % stepsize for Polyak
    
    gammat_alpha_a = (t+1)^(-alpha_a); % stepsize for Polyak 1
    gammat_alpha_b = (t+1)^(-alpha_b); % stepsize for Polyak 2
    gammat_alpha_c = (t+1)^(-alpha_c); % stepsize for Polyak 3
    
    
    
    xtdash = xt + 1; % Because matlab does not take 0 as a valid index.
    xtp1dash = xtp1 + 1;
    atdash = at + 1;
    
    
    
    numofsapair_x = SAPairtoNumber(xtdash,atdash); % Mapping [S A] -> identity
    
    
    % psixt value; It is an long vector with 1 at (x,a)
    
    psixtat = zeros(StaActLen,1);
    psixtat(numofsapair_x) = 1; % Basis are indicator functions
    
    
    
    % Next to look at psixtp1 values; It is a long vector with 1 at all possible x'
    
    psixtp1 = zeros(StaActLen,1);
    for jj = 1:6        
        if(SAPairtoNumber(xtp1dash,jj)>0)
        psixtp1(SAPairtoNumber(xtp1dash,jj)) = 1;        
        end
    end    % psixtp1 is an indicator function, with 1 at all the POSSIBLE future states!
    
    psixtp1_cmp = bitxor(psixtp1,ones_StaActLen);
    psixtp1_cmp_inf = - psixtp1_cmp*10000000000;
    
    
    
    
  
    %%%%%%%%%%%%%%%%% Watkins %%%%%%%%%%%%%%%%%%%%
    
    %1
    Qxt_Wat = Qttheta_Wat'*psixtat;
    Qxtp1_Wat = Qttheta_Wat.*psixtp1 + psixtp1_cmp_inf; % Q values of xtp1, with non-zero entries   
    
%     aaaaaa = (Qxtp1_Wat<0)
    
    maxQxtp1_Wat = max(Qxtp1_Wat);
    
    
    Qttheta_Wat = Qttheta_Wat  + (1.*gammat).*(psixtat.*(R(xtdash,atdash) + beta*maxQxtp1_Wat - Qxt_Wat));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%% SNR 2 b %%%%%%%%%%%%%%%%%%%%
    % 1

    %1
    Qxt_SNR2b = Qttheta_SNR2b'*psixtat;
    Qxtp1_SNR2b = Qttheta_SNR2b.*psixtp1 + psixtp1_cmp_inf; % Q values of xtp1, with non-zero entries   
    
%     aaaaaa = (Qxtp1_Wat<0)
    
    maxQxtp1_SNR2b = max(Qxtp1_SNR2b);
    
    
    Qttheta_SNR2b = min(1000,  Qttheta_SNR2b  + (70.*gammat).*(psixtat.*(R(xtdash,atdash) + beta*maxQxtp1_SNR2b - Qxt_SNR2b))   );    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%%%%%%%%%%%%%%%% SNR 2 c %%%%%%%%%%%%%%%%%%%%
    % 1
    Qxt_SNR2c = Qttheta_SNR2c'*psixtat;
    Qxtp1_SNR2c = Qttheta_SNR2c.*psixtp1 + psixtp1_cmp_inf; % Q values of xtp1, with non-zero entries 
    [max_Qxtp1_SNR2c,optact_SNR2c] = max(Qxtp1_SNR2c);  % the one with the max value
    
    % 2
    psixtp1optutp1_SNR2c = zeros(StaActLen,1);
    psixtp1optutp1_SNR2c(optact_SNR2c) = 1; % basis corresponding to the next state and it's optimal action    
    

    % 3 - A calculation
    Kinv_SNR2c = Kinv_SNR2c + gammat_alpha_c*((psixtat*psixtat' - beta*psixtat*psixtp1optutp1_SNR2c') - Kinv_SNR2c);    
    
    
    % Q update for SNR 2 (a) gain:

    Qttheta_SNR2c  =  Qttheta_SNR2c  +   (   Kinv_SNR2c\(psixtat.*gammat.*(R(xtdash,atdash) + beta*max_Qxtp1_SNR2c - Qxt_SNR2c))  );


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%%%%%%%%%%%%%%%% Polyak Averaging %%%%%%%%%%%%%%%%%%%%
    
    %1
    Qxt_Pol = Qttheta_Pol'*psixtat;
    Qxtp1_Pol = Qttheta_Pol.*psixtp1 + psixtp1_cmp_inf; % Q values of xtp1, with non-zero entries     
    

    Qttheta_Pol_a = min(1000, Qttheta_Pol_a  + gammat_alpha_a.*(psixtat.*(R(xtdash,atdash) + beta*max(Qxtp1_Pol) - Qxt_Pol)));
    Qttheta_Pol = Qttheta_Pol  + gammat.*(Qttheta_Pol_a - Qttheta_Pol);

%     Qttheta_Pol = Qttheta_Pol  + 1.*(psixtat.*(R(xtdash,atdash) + beta*max(Qxtp1_Pol) - Qxt_Pol));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    
    
    %%%%%%%%%%%%%%%%% SPEEDY Q LEARNING %%%%%%%%%%%%%%%%%%%%
    
    gammat_speedy = 1/(t+1);
    
    % 0
    PrevQxt_Speedy = PrevQttheta_Speedy'*psixtat;
    PrevQxtp1_Speedy = PrevQttheta_Speedy.*psixtp1 + psixtp1_cmp_inf;   % Previous Q values of xtp1, with non-zero entries    
    
    TQtm1 = R(xtdash,atdash) + beta*max(PrevQxtp1_Speedy);
    
    % 1
    Qxt_Speedy = Qttheta_Speedy'*psixtat;
    Qxtp1_Speedy = Qttheta_Speedy.*psixtp1 + psixtp1_cmp_inf; % Q values of xtp1, with non-zero entries     
    
    
    TQt = R(xtdash,atdash) + beta*max(Qxtp1_Speedy);
    
    
    PrevQttheta_Speedy = Qttheta_Speedy;
    
    
    Qttheta_Speedy = Qttheta_Speedy  + psixtat.*(     gammat_speedy.*(TQtm1 - Qxt_Speedy)     +        (1-gammat_speedy).*(TQt - TQtm1)    );
%     Qttheta_Speedy = Qttheta_Speedy  + psixtat.*(     gammat_speedy.*(TQtm1 - Qxt_Speedy)     +        (1-gammat_speedy).*(TQt - PrevQxt_Speedy)    );
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
   
    
    
    
    
    
    


    
    xt = xtp1;
    
    psixt = psixtp1;
    

if(mod(t,t_step)==0)

    tt = floor(t/t_step);

for jj = 1:StaActLen
    
    Qttheta_DQ = Qttheta_DQ_a; %(Qttheta_DQ_a + Qttheta_DQ_b)./2;
    
    Qfun_DQ(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = Qttheta_DQ(jj);
    
    Qfun_SNR2a(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = Qttheta_SNR2a(jj);
    Qfun_SNR2b(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = Qttheta_SNR2b(jj);
    Qfun_SNR2c(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = Qttheta_SNR2c(jj);
    
    Qfun_Pol(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = Qttheta_Pol(jj);
    
    Qfun_BVRK(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = Qttheta_BVRK(jj);
    Qfun_Speedy(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = Qttheta_Speedy(jj);

    Qfun_Wat(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = Qttheta_Wat(jj);

end




% Qfun_SNR1
% Qfun_SNR2a
% Qfun_BVRK
% Qfun_myK

% BE ANALYSIS FROM HERE:






%%%%%%%%%%%%%%% FOR Q FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



BE_DQ = zeros(1,StaActLen);

BE_SNR2a = zeros(1,StaActLen);
BE_SNR2b = zeros(1,StaActLen);
BE_SNR2c = zeros(1,StaActLen);

BE_Pol = zeros(1,StaActLen);

BE_BVRK = zeros(1,StaActLen);
BE_SpeedyQ = zeros(1,StaActLen);
BE_Wat = zeros(1,StaActLen);


for jj = 1:StaActLen
        
    x = NumbertoSAPair(jj,1);
    a = NumbertoSAPair(jj,2);
    xdash = x + 1;
    adash = a + 1;
    
    max_Qfun_xtp1_DQ = 0;
    max_Qfun_xtp1_SNR2a = 0;
    max_Qfun_xtp1_SNR2b = 0;
    max_Qfun_xtp1_SNR2c = 0;
    max_Qfun_xtp1_Pol = 0;
    max_Qfun_xtp1_BVRK = 0;
    max_Qfun_xtp1_Speedy = 0;
    max_Qfun_xtp1_Wat = 0;


    NumOfXtp1 = 0;
    
    for ll = 1:StaActLen
        
        if(NumbertoSAPair(ll,1) == x)
            
            xplus1 = NumbertoSAPair(ll,2);
            
            
            max_Qfun_xtp1_DQ = max_Qfun_xtp1_DQ + max(Qfun_DQ(xplus1+1,:));
            max_Qfun_xtp1_SNR2a = max_Qfun_xtp1_SNR2a + max(Qfun_SNR2a(xplus1+1,:));
            max_Qfun_xtp1_SNR2b = max_Qfun_xtp1_SNR2b + max(Qfun_SNR2b(xplus1+1,:));
            max_Qfun_xtp1_SNR2c = max_Qfun_xtp1_SNR2c + max(Qfun_SNR2c(xplus1+1,:));
            max_Qfun_xtp1_Pol = max_Qfun_xtp1_Pol + max(Qfun_Pol(xplus1+1,:));
            max_Qfun_xtp1_BVRK = max_Qfun_xtp1_BVRK + max(Qfun_BVRK(xplus1+1,:));
            max_Qfun_xtp1_Speedy = max_Qfun_xtp1_Speedy + max(Qfun_Speedy(xplus1+1,:));
            max_Qfun_xtp1_Wat = max_Qfun_xtp1_Wat + max(Qfun_Wat(xplus1+1,:));

            
            NumOfXtp1 = NumOfXtp1 + 1;
            
        end
        
    end
    
    max_Qfun_xtp1_DQ = max_Qfun_xtp1_DQ./NumOfXtp1;
    max_Qfun_xtp1_SNR2a = max_Qfun_xtp1_SNR2a./NumOfXtp1;
    max_Qfun_xtp1_SNR2b = max_Qfun_xtp1_SNR2b./NumOfXtp1;
    max_Qfun_xtp1_SNR2c = max_Qfun_xtp1_SNR2c./NumOfXtp1;
    max_Qfun_xtp1_Pol = max_Qfun_xtp1_Pol./NumOfXtp1;
    max_Qfun_xtp1_BVRK = max_Qfun_xtp1_BVRK./NumOfXtp1;
    max_Qfun_xtp1_Speedy = max_Qfun_xtp1_Speedy./NumOfXtp1;
    max_Qfun_xtp1_Wat = max_Qfun_xtp1_Wat./NumOfXtp1;

    
    BE_DQ(jj) = R(xdash,adash) + beta*(  (1-noise)*max(Qfun_DQ(adash,:)) + (noise)*max_Qfun_xtp1_DQ  ) - Qfun_DQ(xdash,adash);
    
    BE_SNR2a(jj) = R(xdash,adash) + beta*((1-noise)*max(Qfun_SNR2a(adash,:)) +  (noise)*max_Qfun_xtp1_SNR2a) - Qfun_SNR2a(xdash,adash);
    BE_SNR2b(jj) = R(xdash,adash) + beta*((1-noise)*max(Qfun_SNR2b(adash,:)) +  (noise)*max_Qfun_xtp1_SNR2b) - Qfun_SNR2b(xdash,adash);
    BE_SNR2c(jj) = R(xdash,adash) + beta*((1-noise)*max(Qfun_SNR2c(adash,:)) +  (noise)*max_Qfun_xtp1_SNR2c) - Qfun_SNR2c(xdash,adash);
    
    BE_Pol(jj) = R(xdash,adash) + beta*((1-noise)*max(Qfun_Pol(adash,:)) +  (noise)*max_Qfun_xtp1_Pol) - Qfun_Pol(xdash,adash);
    
    BE_BVRK(jj) = R(xdash,adash) + beta*((1-noise)*max(Qfun_BVRK(adash,:)) +  (noise)*max_Qfun_xtp1_BVRK) - Qfun_BVRK(xdash,adash);
    BE_SpeedyQ(jj) = R(xdash,adash) + beta*((1-noise)*max(Qfun_Speedy(adash,:)) +  (noise)*max_Qfun_xtp1_Speedy) - Qfun_Speedy(xdash,adash);
    
    BE_Wat(jj) = R(xdash,adash) + beta*((1-noise)*max(Qfun_Wat(adash,:)) +  (noise)*max_Qfun_xtp1_Wat) - Qfun_Wat(xdash,adash);

end


Qttheta_DQ_BEevolQ(tt) = max(abs(BE_DQ));

Qttheta_SNR2a_BEevolQ(tt) = max(abs(BE_SNR2a));
Qttheta_SNR2b_BEevolQ(tt) = max(abs(BE_SNR2b));
Qttheta_SNR2c_BEevolQ(tt) = max(abs(BE_SNR2c));


Qttheta_Pol_BEevolQ(tt) = max(abs(BE_Pol));


Qttheta_BVR_BEevolQ(tt) = max(abs(BE_BVRK));
Qttheta_Speedy_BEevolQ(tt) = max(abs(BE_SpeedyQ));

Qttheta_Wat_BEevolQ(tt) = max(abs(BE_Wat));





%%%%%%%%%%%%%%% FOR VALUE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





hStar_SNR1 = zeros(1,6);

hStar_SNR2a = zeros(1,6);
hStar_SNR2b = zeros(1,6);
hStar_SNR2c = zeros(1,6);

hStar_Pol = zeros(1,6);

hStar_BVRK = zeros(1,6);
hStar_SpeedyQ = zeros(1,6);

hStar_Wat = zeros(1,6);

for kk = 1:6


hStar_SNR1(kk) = max(Qfun_DQ(kk,:));


hStar_SNR2a(kk) = max(Qfun_SNR2a(kk,:));
hStar_SNR2b(kk) = max(Qfun_SNR2b(kk,:));
hStar_SNR2c(kk) = max(Qfun_SNR2c(kk,:));


hStar_Pol(kk) = max(Qfun_Pol(kk,:));

hStar_BVRK(kk) = max(Qfun_BVRK(kk,:));
hStar_SpeedyQ(kk) = max(Qfun_Speedy(kk,:));

hStar_Wat(kk) = max(Qfun_Wat(kk,:));

end


StaLen = 6;

BEh_SNR1 = BECalc_Policy_minu(hStar_SNR1,R,beta,StaLen,noise);


BEh_SNR2a = BECalc_Policy_minu(hStar_SNR2a,R,beta,StaLen,noise);
BEh_SNR2b = BECalc_Policy_minu(hStar_SNR2b,R,beta,StaLen,noise);
BEh_SNR2c = BECalc_Policy_minu(hStar_SNR2c,R,beta,StaLen,noise);



BEh_Pol = BECalc_Policy_minu(hStar_Pol,R,beta,StaLen,noise);


BEh_BVRK = BECalc_Policy_minu(hStar_BVRK,R,beta,StaLen,noise);
BEh_SpeedyQ = BECalc_Policy_minu(hStar_SpeedyQ,R,beta,StaLen,noise);
BEh_Wat = BECalc_Policy_minu(hStar_Wat,R,beta,StaLen,noise);






Qttheta_DQ_BEevol(tt) = max(abs(BEh_SNR1));

Qttheta_SNR2a_BEevol(tt) = max(abs(BEh_SNR2a));
Qttheta_SNR2b_BEevol(tt) = max(abs(BEh_SNR2b));
Qttheta_SNR2c_BEevol(tt) = max(abs(BEh_SNR2c));


Qttheta_Pol_BEevol(tt) = max(abs(BEh_Pol));


Qttheta_BVR_BEevol(tt) = max(abs(BEh_BVRK));
Qttheta_Speedy_BEevol(tt) = max(abs(BEh_SpeedyQ));

Qttheta_Wat_BEevol(tt) = max(abs(BEh_Wat));



% 
% 
% Qttheta_DQ_BEAvgevol(tt) = sum(abs(BEh_SNR1));
% 
% Qttheta_SNR2a_BEAvgevol(tt) = sum(abs(BEh_SNR2a));
% Qttheta_SNR2b_BEAvgevol(tt) = sum(abs(BEh_SNR2b));
% Qttheta_SNR2c_BEAvgevol(tt) = sum(abs(BEh_SNR2c));
% 
% 
% Qttheta_Pol_BEAvgevol(tt) = sum(abs(BEh_Pol));
% 
% 
% Qttheta_BVR_BEAvgevol(tt) = sum(abs(BEh_BVRK));
% Qttheta_Speedy_BEAvgevol(tt) = sum(abs(BEh_SpeedyQ));
% 
% Qttheta_Wat_BEAvgevol(tt) = sum(abs(BEh_Wat));



t

end

end

end


save('SIMData_PAPER')


TVEC2 = 1:T_small;

 

for jj = 1:StaActLen
    
    BEmatrix_SNR1(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = BE_DQ(jj);
    
    BEmatrix_SNR2a(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = BE_SNR2a(jj);
    BEmatrix_SNR2b(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = BE_SNR2b(jj);
    BEmatrix_SNR2c(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = BE_SNR2c(jj);
    
    BEmatrix_Pol(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = BE_Pol(jj);
    
    BEmatrix_BVRK(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = BE_BVRK(jj);
    BEmatrix_SpeedyQ(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = BE_SpeedyQ(jj);
    
    BEmatrix_Wat(NumbertoSAPair(jj,1)+1,NumbertoSAPair(jj,2)+1) = BE_Wat(jj);


end


BEmatrix_SNR1

BEmatrix_SNR2a
BEmatrix_SNR2b
BEmatrix_SNR2c

BEmatrix_Pol

BEmatrix_BVRK
BEmatrix_SpeedyQ

BEmatrix_Wat


TVEC1 = t_step:t_step:T;

TVEC2 = 1:T_small;


figure
plot(TVEC1,Qttheta_DQ_BEevol(TVEC2),TVEC1,Qttheta_SNR2a_BEevol(TVEC2),TVEC1,Qttheta_SNR2b_BEevol(TVEC2),TVEC1,Qttheta_SNR2c_BEevol(TVEC2),TVEC1,Qttheta_Pol_BEevol(TVEC2),TVEC1,Qttheta_BVR_BEevol(TVEC2),TVEC1,Qttheta_Speedy_BEevol(TVEC2),TVEC1,Qttheta_Wat_BEevol(TVEC2))
legend('DQ','SNR 2 \alpha = 0.7','SNR 2 \alpha = 0.8','SNR 2 \alpha = 0.9','Polyak Averaging','BVR','Speedy Q learning','Watkins')
xlabel('t')
ylabel('MaxBEh_{\theta(t)}')
title('Max BE of the optimal value function')


% figure
% plot(TVEC,Qttheta_DQ_BEAvgevol(TVEC),TVEC,Qttheta_SNR2a_BEAvgevol(TVEC),TVEC,Qttheta_SNR2b_BEAvgevol(TVEC),TVEC,Qttheta_SNR2c_BEAvgevol(TVEC),TVEC,Qttheta_Pol_BEAvgevol(TVEC),TVEC,Qttheta_BVR_BEAvgevol(TVEC),TVEC,Qttheta_Speedy_BEAvgevol(TVEC),TVEC,Qttheta_Wat_BEAvgevol(TVEC))
% legend('DQ','SNR 2 \alpha = 0.7','SNR 2 \alpha = 0.8','SNR 2 \alpha = 0.9','Polyak Averaging','BVR','Speedy Q learning','Watkins')


figure
plot(TVEC1,Qttheta_DQ_BEevolQ(TVEC2),TVEC1,Qttheta_SNR2a_BEevolQ(TVEC2),TVEC1,Qttheta_SNR2b_BEevolQ(TVEC2),TVEC1,Qttheta_SNR2c_BEevolQ(TVEC2),TVEC1,Qttheta_Pol_BEevolQ(TVEC2),TVEC1,Qttheta_BVR_BEevolQ(TVEC2),TVEC1,Qttheta_Speedy_BEevolQ(TVEC2),TVEC1,Qttheta_Wat_BEevolQ(TVEC2))
legend('DQ','SNR 2 \alpha = 0.7','SNR 2 \alpha = 0.8','SNR 2 \alpha = 0.9','Polyak Averaging','BVR','Speedy Q learning','Watkins')
xlabel('t')
ylabel('MaxBEQ_{\theta(t)}')
title('Max BE of the optimal Q function')

