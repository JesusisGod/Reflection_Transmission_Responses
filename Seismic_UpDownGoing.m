function accurate_trace=Seismic_UpDownGoing(c,NT,UD_FLAG)
%% Written by Jide Ogunbo
%% Reference: Robinson, E. A., and S. Trietel, 1980, Geophysical signal analysis. Prentice-Hall, Inc., 466pp. (Chapter 13)

% c  =reflection coefficient
% NT = number of time samples
% UD_FLAG=1 for upgoing (Reflected) wave; else for downgoing (Transmitted) wave

%% Example (Using information in Table 1 of the manuscript)
% vel=[1500 3000 1500 2000 1750 2750];
% rho=[1000 2250 1000 2000 1500 2000];
% AI1=vel(1:end-1).*rho(1:end-1); AI2=vel(2:end).*rho(2:end);
% c=(AI2-AI1)./(AI2+AI1);
% NLYR=length(c);
% NT=2*NLYR; 
% TrueReflData=Seismic_UpDownGoing(c,NT,1); % Reflected data
% TrueTranData=Seismic_UpDownGoing(c,NT,2); % Transmitted data

Lc1=length(c);
if Lc1==1
    c(2)=0; 
end
NL=length(c); % Number of layers
R=c(NL);

for i=NL-1:NL-1
    nume=[c(i) R]; % Polynomial coefficients 
    denoRHS=c(i)*R;
    deno=[1 denoRHS]; 
end



if NL>2    
    for i=NL-2:-1:1
        new_nume_LHS=conv(c(i),deno); 
        new_nume_RHS=[0 nume]; % 0 indicates time shift consequence of multiplication Z                               
        new_deno_LHS=deno;
        new_deno_RHS=c(i)*[0 nume];                                         

        deno_store(1:length(new_deno_LHS),1)=new_deno_LHS;
        deno_store(1:length(new_deno_RHS),2)=new_deno_RHS;

        nume_store(1:length(new_deno_LHS),1)=new_nume_LHS;
        nume_store(1:length(new_deno_RHS),2)=new_nume_RHS;

        nume=sum(nume_store,2);
        deno=sum(deno_store,2);

        nume=nume';
        deno=deno';
    end
end

denoRHS=deno(2:end);
trace=[0 denoRHS]; 
prev_trace=1;
trace=conv(trace,prev_trace);
denoRHS=trace;
    
STORE(1)=1;
STORE(1:length(trace),2)=-trace; % handling Binomial first term =-x^1
   
if NT>2
    for j=3:NT
        pow=(j-1);
        for k=1:j-2
            denoRHS=conv(trace,denoRHS); 
        end
            STORE(1:length(denoRHS),j)=(-1)^pow*denoRHS; % Binomial part =[-x+x^2-x^3+...]
            denoRHS=trace; 
    end
end
    
    
SUMS=sum(STORE,2);
BINOMIAL(1:length(SUMS),1)=SUMS; %STORE;
    if UD_FLAG==1
        [~,~,Rtmp]=autocorr_convo(nume,BINOMIAL.');
        accurate_trace=Rtmp(1:NT,1);
    else

    %% Transmitted
        tn=1+c;       
        tns=prod(tn);
        [~,~,RtmpTransmitted]=autocorr_convo(tns,BINOMIAL.');
        trace=RtmpTransmitted(1:NT,1);
        NL1=NL-1;
        z=zeros(NL1,1); % for z^(0.5*(NL-1)) for timeshift on transmitted signal
        timeshift=z;
        if Lc1==1
            accurate_trace=trace;
        else
            accurate_trace=[timeshift;trace];
        end
    end