function results = System_Model(optIn)
% vers3: massive MIMO uplink
% A: is the singal
% X: is the sparse channel
% M: coherence time
% N: Nr receiver antenna
% L: Nt transmit antenna

% version 4: 
% for each channel realization, we run multiple samples of noise

% version 7:
% theortical results + numerical analysis

% version 10:
% subspace estimation before DL, i.e., Y = UDV', Z = U1D1V1'. Then D1V1' is regarded as the input of DL algorithm. 

%% =====================================================================
%
%trial_DL: This function runs several algorithms on a sample
%instance of the dictionary learning problem. The function can be run with no
%arguments. See the nargin==0 block below for the required fields to call
%the function with optional parameters. The function returns a structure of
%algorithm results.
%This code can be used to produce the results for the noise free phase
%plane plots for dictionary learning in the BiG-AMP arXiv submission. 
%(This code conducts a single trial at 1 point in the phase plane.)

%Add paths
%setup_DL

%Test case
if nargin == 0
    clc
    
    %Handle random seed
    defaultStream = RandStream.getGlobalStream;
    if 1
        savedState = defaultStream.State;
        save random_state.mat savedState;
    else
        load random_state.mat %#ok<UNRCH>
    end
    defaultStream.State = savedState;
    
    %Flags to test BiG-AMP variants
    optIn.tryBigampEM = 1;
    
    %Flags to turn on comparison algorithms
    optIn.tryKsvd = 0;
    optIn.tryErspud = 0;
    optIn.trySpams = 0;

    
    %Problem dimensions
    optIn.M = 20; %coherence time
    optIn.N = optIn.M; %transmit antenna
    optIn.L = ceil(5*optIn.N*log(optIn.N)); %received antenna why set in this way?
    
    
    %Specify dictionary usage
    optIn.K = 5;
    
    %Specify maximum allowed trials
    optIn.maxTrials = 1;
    
    %SNR
    optIn.SNR = inf;
    
    %Specify coding mode (0 for OMP, 1 for TST)
    %The specified function is used to compute a coding error for the
    %resulting dictionary for each algorithm. TST may be significantly
    %slower.
    optIn.useTST = 0;
    
    %Precondition option (enable to use a preconditioner for EM-BiG-AMP)
    optIn.precondition = 0;
    
    
end

%% Problem Setup

%Turn algs on/off
tryBigampEM = optIn.tryBigampEM;
tryKsvd = optIn.tryKsvd;
tryErspud = optIn.tryErspud;
trySpams = optIn.trySpams;

%SNR
SNR = optIn.SNR;

%Coding
useTST = optIn.useTST;

%Precondition
precondition = optIn.precondition;

%Max trials
maxTrials = optIn.maxTrials;

%Define problem dimensions
M = optIn.M;
L = optIn.L;
N = optIn.N;
K = optIn.K;
P = optIn.P;
qam = optIn.qam;


%Set options
opt = BiGAMPOpt; %initialize the options object
opt.nit = optIn.nit; %limit iterations

opt.uniformVariance = optIn.uniformVariance;

%Set sizes
problem = BiGAMPProblem();
problem.M = M; %%%%% only for subspace estimation
problem.N = N;
problem.L = L;


%% Build the dictionary



%% Compute coefficient vectors

%Compute true vectors with exactly K non-zeroes in each %? does X is a
%sparse matrix, which is different from the paper where the H is sparse.
%X = randn(N,L);
% Kcol = ceil(K*L/N);
% for ll = 1:N
%     yada = randperm(L);
%     yada2 = zeros(1,L);
%     yada2(yada(1:Kcol)) = 1;
%     X(ll,:) = X(ll,:) .* yada2;
% end
X = sqrt(0.5)*complex(randn(N,L),randn(N,L));
for ll = 1:L
    yada = randperm(N);
    yada2 = zeros(N,1);
    yada2(yada(1:K)) = 1;
    X(:,ll) = X(:,ll) .* yada2;
end

% fading_factor = -20;     % the dB range
% pathLoss = rand(N, 1) * fading_factor ;
% pathLoss = 10.^(pathLoss / 10);
% pathLoss = sqrt(pathLoss / max(pathLoss)); %normalize average Avar to 1
% X = diag(pathLoss)*X ;
%Draw randomly
%% 1-bit pilot % A is the transmitted signal 
% A = zeros(M,N);
% A1 = randn(M-1,N); % singal
% A1 = sqrt(M-1)*A1*diag(1 ./ sqrt(abs(diag(A1'*A1))));
% A(1,:) = 1;
% A(2:M,:) = A1;
% scaledFactor = sqrt(M*N/real(trace(A*A')));
% A = scaledFactor*A; % normalized the signal power to M*N, each element 1.

if qam == 4
    Sam = [1+1i -1+1i -1-1i 1-1i]./sqrt(2);    
elseif qam == 16
    Sam = [1+1i -1+1i -1-1i 1-1i 3+1i -1+3i -3-1i 1-3i 3+3i -3+3i -3-3i 3-3i 1+3i -3+1i -1-3i 3-1i]./sqrt(10);
end
B=randi([0,1], sqrt(qam)*(M-P), N);
B1=randi([0,1], sqrt(qam)*P, N);
A = zeros(M,N);

B1(1:sqrt(qam) ,:)=1;
 for n=1:N
     n1=n;
   for m=(sqrt(qam) + 1) :(sqrt(qam) + ceil(log2(N)))        %qam+
        B1(m,n)=mod(n1,2); 
        n1=floor(n1/2);
   end
 end
 D = Constell_Modulate(B1, Sam);
 A(1:P,:) = D;
 A(1+P:M,:) = Constell_Modulate(B, Sam);
 P_num = Pilot_Num(D, Sam); 
% B = randi([0,1], 2, M);
% B(:,1) = 1;
% A = Constell_Modulate(B, 4);
% A = A.' ;
%Normalize the columns
%A = sqrt(M)*A*diag(1 ./ sqrt(abs(diag(A'*A))));

%Dictionary error function ? what's mean for this function 
dictionary_error_function =...
    @(q) 20*log10(norm(A -...
    q*find_permutation(A,q),'fro')/norm(A,'fro'));


%scaledFactor = sqrt(N*L/real(trace(X*X')));
%X = scaledFactor*X;

%     X_temp = randn(N,L-N);
%     for ll = 1:L-N
%         yada = randperm(N);
%         yada2 = zeros(N,1);
%         yada2(yada(1:K)) = 1;
%         X_temp(:,ll) = X_temp(:,ll) .* yada2;
%     end
%     X = [eye(N), X_temp];

%% Form the output channel

%Compute noise free output
Z = A*X;

%Define the error function
error_function = @(qval) 20*log10(norm(qval - Z,'fro') / norm(Z,'fro'));
opt.error_function = error_function;


%Determine nuw


nuw = 10^(-SNR/10)*N;
nuw = max(nuw,1e-20);
%nuw = norm(reshape(Z,[],1))^2/M/L*10^(-SNR/10);

%Initialize results as empty
results = [];
PrioriIn.noiseVar = nuw;
PrioriIn.lambda = K/N;
PrioriIn.nuX = 1;
%PrioriIn.nuX =  1 ;
PrioriIn.spar = K;
PrioriIn.A = A;
PrioriIn.X = X;
PrioriIn.D = D;
PrioriIn.P = P;
PrioriIn.B = B;
PrioriIn.Sam = Sam;
PrioriIn.P_num = P_num;

%% EM BiG AMP

if tryBigampEM
    
    %Silence
    opt.verbose = false;
 %   disp('Starting EM-BiG-AMP')
    
    %Coding
%     bestError = inf;
%     bestSparsity = inf;
%     
%     %Compute preconditioner
    if precondition
        Q = chol(inv(Y*Y'));  %? what's mean for Q?
    else
        Q = 1;
    end
    tstart = tic; % start the programme 
    Xerror = 0;
    
    
    for Nsample = 1:maxTrials
        %Noisy output channel
         Y = Z + sqrt(nuw/2)*(randn(size(Z))+1i*randn(size(Z)));
        
        % subspace estimation
%         [Uy, Dy, Vy] = svd(Y);
%         Un = Uy(:,1:N);
%         Dn = Dy(1:N,1:N);
%         Vn = Vy(:,1:N);
%         eqvZ = Dn*Vn'; %% why write it in this way? 
       % PrioriIn.Un = Un;
        eqvZ = Y;
        %Coding error
        coding_error_function = @(q) 20*log10(coding_error(eqvZ,q,K,useTST));
        
        % Define the error function
      %  QZ = Q*eqvZ; %Using (possibly) noisy data for error function
        %QZ = Q*Y;
        error_function2 = @(qval) 20*log10(norm(qval - eqvZ,'fro') / norm(eqvZ,'fro'));
        opt.error_function = error_function2;
        
        for trial = 1:1
            %Run EM-BiGAMP
            %[estFinTemp,~,~,estHistEMtemp] = EMBiGAMP_DL(Q*Y,problem,opt);
            
            %[estFinTemp,~,~,estHistEMtemp] = BiGAMP_DL_fixEM(Q*Y,problem,opt,PrioriIn);
            [estFinTemp,~,~,estHistEMtemp] = BiGAMP_EM(eqvZ, problem,opt,PrioriIn);   %Q*eqvZ
            %Correct the dictionary
            estFinTemp.Ahat = Q \ estFinTemp.Ahat;
            
            %Reset error function
            opt.error_function = error_function;
            
            %If the sparisty is better and the error is very small or better
%             if (sum(sum(estHistEMtemp.p1)) < bestSparsity) && ...
%                     ( (estHistEMtemp.errZ(end) < bestError)  ||...
%                     (estHistEMtemp.errZ(end) < -100)  )
                
                %Update
                bestSparsity = sum(sum(estHistEMtemp.p1));
                bestError = estHistEMtemp.errZ(end);
                AhatOptEM = estFinTemp.Ahat;
                XhatOptEM = estFinTemp.xhat;
                estHistEM = estHistEMtemp;
                p1EM = estHistEMtemp.p1;
                
                %Notify user
%                 disp(['Accepting new result. Error: ' num2str(bestError)...
%                     '  Average sparsity: ' num2str(bestSparsity/L)...
%                     '  Max Sparsity: ' num2str(max(sum(p1EM)))])
%            end
            %keyboard
        end
        
        
        %Save results
%        tempAhatOptEM = Un*AhatOptEM;
        tempAhatOptEM= AhatOptEM;
        %Save results
        %tempAhatOptEM = Un*AhatOptEM;
        %tempAhatOptEM = AhatOptEM;


        %% LMMSE estimator for scaling ambiguity
          pMatrix2 = find_permutation(A,tempAhatOptEM);
        %[Ahat, pMatrix2] = find_permu(tempAhatOptEM, P,Sam);
           
      
        %%
        Ahat = tempAhatOptEM*pMatrix2;
        Xhat = pinv(pMatrix2)*XhatOptEM;
        %/////////////////////////
         for m=1:M
             for n=1:N
                 [~,I]=min(abs(Sam-Ahat(m,n)));
                 Ahat(m,n)=Sam(I);
              end
         end
        Bhat = Constell_Mapping(Ahat(1+P:M, :), Sam);          

        L=length(find((B-Bhat)==0));
        errbate= (sqrt(qam)*(M-P)*N-L)/(sqrt(qam)*(M-P)*N);
        Xerror = (norm(X - Xhat,'fro')^2/norm(X,'fro')^2);

        
    end
    tEMGAMP = toc(tstart);
    
    loc = length(results) + 1;
    results{loc}.name = 'EM-BiG-AMP'; %#ok<*AGROW>
    results{loc}.Xerror = Xerror;
    results{loc}.errbate = errbate;
end

%% SPAMS


if trySpams
    
    %Build params
    spams_param = [];
    spams_param.K = N;
    spams_param.mode = 2;
    spams_param.lambda = 0.1/sqrt(N);
    spams_param.iter = 1000;
    
    
    %Trials
    bestSPAMSerror = inf;
    
    tstart = tic;
    for trial = 1:maxTrials
        
        %Do it   
        A_spamsTemp = mexTrainDL(Y,spams_param);
        
        
        SPAMSerrorTemp = dictionary_error_function(A_spamsTemp);
        
        %Update the estimate if this result was superior
        if SPAMSerrorTemp < bestSPAMSerror
            SPAMSerror = SPAMSerrorTemp;
            bestSPAMSerror = SPAMSerror;
            A_spams = A_spamsTemp;
            disp(['Updating Solution. Error was: '...
                num2str(SPAMSerror) ' dB'])
        end
        
    end
    tspams = toc(tstart);
    
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'SPAMS'; %#ok<*AGROW>
    results{loc}.err = dictionary_error_function(A_spams);
    results{loc}.time = tspams;
    results{loc}.errHist = results{loc}.err;
    results{loc}.timeHist = zeros(size(results{loc}.errHist));
    results{loc}.dict = A_spams;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
end



%% ER-SpUD


if tryErspud
    
    
    %Call it
    tic
    %[Aerspud,Xerspud] = ER_SpUD_SC(Y);
    [Aerspud,Xerspud]=dl_spud(Y);   %#ok<NASGU>
    tErspud = toc;
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'ER-SpUD (proj)'; %#ok<*AGROW>
    results{loc}.err = dictionary_error_function(Aerspud);
    results{loc}.time = tErspud;
    results{loc}.errHist = results{loc}.err;
    results{loc}.timeHist = zeros(size(results{loc}.errHist));
    results{loc}.dict = Aerspud;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
    
end


%% K-SVD

if tryKsvd
    
    %Build params
    ksvd_params = [];
    ksvd_params.data = Y;
    ksvd_params.Tdata = K;
    ksvd_params.dictsize = N;
    ksvd_params.iternum = 100;
    ksvd_params.exacty = 1;
    
    %Trials
    bestKSVDerror = inf;
    tstart = tic;
    for trial = 1:maxTrials
        
        %Do it
        
        [A_ksvdtemp,~,err_ksvdtemp] = ksvd(ksvd_params);
        
        
        %Fix error
        err_ksvdtemp = err_ksvdtemp.^2 * numel(Y);
        
        %Update the estimate if this result was superior
        if err_ksvdtemp(end) < bestKSVDerror
            bestKSVDerror = err_ksvdtemp(end);
            err_ksvd = err_ksvdtemp;
            A_ksvd = A_ksvdtemp;
            disp(['Updating Solution. Error was: '...
                num2str(10*log10(err_ksvd(end)/norm(Z,'fro')^2))])
        end
        
    end
    tksvd = toc(tstart);
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'K-SVD'; %#ok<*AGROW>
    results{loc}.err = 10*log10(err_ksvd(end)/norm(Z,'fro')^2);
    results{loc}.time = tksvd;
    results{loc}.errHist = 10*log10(err_ksvd/norm(Z,'fro')^2);
    results{loc}.timeHist = zeros(size(err_ksvd));
    results{loc}.dict = A_ksvd;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
end




%% Store the options structures in results
results{1}.optIn = optIn;
results{1}.trueDict = A;
results{1}.trueEncoding = X;


%% Show Results

if nargin == 0
    results{:} %#ok<NOPRT>
    disp('Note that dictError is the normalized error in dB for recovering')
    disp('the dictionary for each algorithm')
end

