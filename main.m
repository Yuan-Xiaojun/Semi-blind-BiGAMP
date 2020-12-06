
% This is a simple example to apply EM-BIGAMP to wireless communications

%% 1bit LMMSE to eliminate ambiguity
% the ambiguity elimination method in this programme is not practical
% available. It is just for research test.

clear all;
tic
randn('state',1);
rand('state',1);
%results = DLCommun;
%setup_DL
%Handle random seed
defaultStream = RandStream.getGlobalStream; %%£¿ the function for this code 
if 1
    
    savedState = defaultStream.State; 
    save random_state.mat savedState;
else
    load random_state.mat %#ok<UNRCH>
end
defaultStream.State = savedState; %?


Ntrange = [20]; % transmit antenna
SNRrange = [20];
SparseK = [0.30]; % sparsity ratio per column
SimLen = 1;%?
Nrrange = [128];%ceil(5.*Nrange.*log(Nrange));% [100]; % receive antenna
P = [3];
QAM = 16;

Threshold = 1e-3; % 

%Flags to test BiG-AMP variants
optIn.tryBigampEM = 1;

%Flags to turn on comparison algorithms
optIn.tryKsvd = 0;
optIn.tryErspud = 0;
optIn.trySpams = 0;
optIn.maxTrials = 1; % noise sample per channel realization


for iK = 1:length(SparseK)
    %optIn.K = SparseK(iK);
    for iN = 1:length(Nrrange)
        optIn.N = Ntrange(1); %Ntrange(iN); %
        optIn.P = P(1);
        optIn.M = 50;%optIn.N; % coherence time
        optIn.K = round(SparseK(iK)*optIn.N);
        optIn.qam = QAM;
        Llower = 180;
        Lupper = 0;
        Nr = Nrrange(iN);
        
        %while(1)
            optIn.L = Nr;
            
            %optIn.L = Nrrange(iL);%ceil(5*optIn.N*log(optIn.N)); %From Spielman et al.
            %optIn.L = ceil(5*optIn.N*log(optIn.N)); %From Spielman et al.
            errbate = zeros(1,length(SNRrange));
            XMSE = zeros(1,length(SNRrange)); % average over channel samples
            terrbate = zeros(1,length(SNRrange));
            
            for sim = 1:SimLen %?
                tAMSE = zeros(1,length(SNRrange));
                tXMSE = zeros(1,length(SNRrange));
                
                for iS = 1:length(SNRrange)
                    %Specify maximum allowed trials
                    %optIn.maxTrials = 1;
                    %SNR
                    optIn.SNR = SNRrange(iS); % un-use in this program
                    %Specify coding mode (0 for OMP, 1 for TST)
                    %The specified function is used to compute a coding error for the
                    %resulting dictionary for each algorithm. TST may be significantly
                    %slower.
                    optIn.useTST = 0;
                    
                    %Precondition option (enable to use a preconditioner for EM-BiG-AMP)
                    optIn.precondition = 0;
                   
                    optIn.nit = 100; %limit iterations
                    
                    
                    optIn.uniformVariance = 0; % 1 for uniform variance; 0 for general case
                               results = System_Model(optIn); % subspace estimation

                         terrbate(iS) = results{1}.errbate;

                         tXMSE(iS) = results{1}.Xerror;
                        
                end

                    XMSE = XMSE.*(sim-1)/sim + tXMSE./sim;
                    errbate = errbate.*(sim-1)/sim + terrbate./sim;

                    fprintf('sim = %4d, K = %2d, N = %3d, M = %3d,  L = %4d  ,P = %2d\n', sim, optIn.K, optIn.N, optIn.M, optIn.L,optIn.P);
                   
                    fprintf('SNR =');
                    for iS=1:length(SNRrange)
                        fprintf('%3d  ', SNRrange(iS));
                    end
                    fprintf('\n');
%                     fprintf('E =');
% 
%                     for iS=1:length(SNRrange)
%                         fprintf('%3d  ', E(iS));
%                     end
%                     fprintf('\n');
                    fprintf('tBerrbate =');
                    for iS=1:length(SNRrange)
                        fprintf('%3.6e  ', terrbate(iS));
                    end
                    fprintf('\n');
                    
                    fprintf('Berrbate =');
                    for iS=1:length(SNRrange)
                        fprintf('%3.6e  ', errbate(iS));
                    end
                    fprintf('\n');

                    fprintf('tXMSE =');
                    for iS=1:length(SNRrange)
                        fprintf('%3.6e  ', tXMSE(iS));
                    end
                    fprintf('\n');
                    
                    fprintf('XMSE =');
                    for iS=1:length(SNRrange)
                        fprintf('%3.6e  ', XMSE(1,iS));
                    end
                    fprintf('\n');
                    fprintf('\n');

            end
                fid = fopen('SCSE_esults.txt','a+');
                fprintf(fid,'\n');
                fprintf(fid,'T = %3d, Nt = %3d, Nr = %4d, K = %2d, simLen = %4d\n',optIn.M, optIn.N, optIn.L,optIn.K,SimLen);

                fprintf(fid,'snr = ');
                for iS = 1:length(SNRrange)
                    fprintf(fid,'%3d  ', SNRrange(iS));
                end
                fprintf(fid,'\n');

                fprintf(fid,'Berrbate = ');
                for iS = 1:length(SNRrange)
                    fprintf(fid,'%3.6e  ', errbate(iS));
                end            
                fprintf(fid,'\n');


                fprintf(fid,'XMSE =');
                for iS=1:length(SNRrange)
                    fprintf(fid,'%3.6e  ', XMSE(1,iS));
                end
                fprintf(fid,'\n');

         
    end
end
toc

            