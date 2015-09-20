function [Coh, K] = coherence_hyp(truesp, simsp, N_sample_max)
	%***BEN: THIS IS FOR ONE SEGMENT - YOU NEED TO PUT THIS IN A "SUPER-LOOP" THAT COMPUTES
	%THE COHERENCE - A COMPLEX NUMBER - FOR EACH SEGMENT AND ADDS THE COMPLEX
	%COHERENCES TOGETHER and DIVIDES BY "N" TO COMPUTE THE AVERAGE OVER
	%SEGMENTS. THEN PLOT THE ABSOLUTE VALUE, A NUMBER BETWEEN ZERO AND ONE.
	%NOTE 1 - BY PADDING, I AM MAKING ALL SPECTRAL TRACES OF EQUAL LENGTH.
	%NOTE 2 - BY USING TAPERS, I CAN DEFINE THE AVERAGING OVER TAPERS IN THE
	%SPECTRAL ESTIMATE BASED ON A FIXED BANDWIDTH (WHICH I CHOOSE TO BE
	%f_Nyquist/25 = 2 Hz); THIS PERMITS ADDTION OF SPECTRA FOR DIFFERENT LENGTH
	%TIME SERIES***
	N_sample=length(truesp) ; %THIS MUST BE CALCUALTED FOR EACH SET OF TRACES
	F_sample=100; % Sampling rate
	Pad = 2^(1+nextpow2(N_sample_max));  
	% pad to > 2-times data length - ***BEN: YOU SHOULD FIX AT THE VALUE CALCULATED FOR THE LARGEST N_sample IN THE DATA SET***
	f=linspace( 0, 1, Pad )*F_sample; %Frequency base
	NW=N_sample*0.5*(1/25); %Time-bandwidth product with bandwidth as a fraction (chose as 1/20) of f_Nyquist (0.5 in computer units) on a sample-by-sample basis
	K=fix(2*NW-1); % Degrees of freedom or 2*p-1
	[E,V] = dpss(N_sample,NW,K); % Family of multiple taper functions
	FT1=zeros(Pad,K);
	FT2=zeros(Pad,K);
	Pwr1=zeros(Pad,K);
	Pwr2=zeros(Pad,K);
	XPwr=zeros(Pad,K);
	Pwr1_tot=zeros(Pad,1);
	Pwr2_tot=zeros(Pad,1);
	XPwr_tot=zeros(Pad,1);
	Coh=zeros(Pad,1);
	for k=1:K; % sum over individual estimates of Power, with normalization corrections
	        FT1(:,k)=fft(E(:,k).*(truesp-mean(truesp)),Pad);
	        FT2(:,k)=fft(E(:,k).*(simsp-mean(simsp)),Pad);
	        Pwr1(:,k)=FT1(:,k).*conj(FT1(:,k));
	        Pwr2(:,k)=FT2(:,k).*conj(FT2(:,k));
	        XPwr(:,k)=FT1(:,k).*conj(FT2(:,k));
	        Pwr1_tot(:,1)=Pwr1_tot(:,1)+Pwr1(:,k);
	        Pwr2_tot(:,1)=Pwr2_tot(:,1)+Pwr2(:,k);
	        XPwr_tot(:,1)=XPwr_tot(:,1)+XPwr(:,k);
	end;
	Coh=XPwr_tot./sqrt(Pwr1_tot.*Pwr2_tot);
end