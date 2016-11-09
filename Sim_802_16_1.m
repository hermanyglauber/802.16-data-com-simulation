clear;
close;
num_b = 1000000; %number of bits to be simulated
k = 1;
bits=randi(2, 1, num_b)-1; %bits is a random real vector = 0 or +1 
display(bits);


%WirelessMAN-SC?
%802.16 bits encode strategy? 
%codigos concatenados Codigo externo + codigo interno: Reed-Solomon (com t ajustável), paridade, código convolucional e Block Turbo Code

%802.16 signal modulation strategy QPSK
%----------------------------inicio modulação QPSK------------------------%
for i=1:length(bits)/2
    if(bits(2*i-1)==0)
        if(bits(2*i)==0)
            mod(k)=complex(1,1); %simbolo = 00 mapeia para complexo 1+1i
        else
            mod(k)=complex(-1,1); %simbolo = 01 mapeia para complexo -1+1i
        end
    else
        if(bits(2*i)==0)
            mod(k)=complex(-1,-1); %simbolo = 10 mapeia para complexo -1-1i
        else
            mod(k)=complex(1,-1); %simbolo = 11 mapeia para complexo 1-1i
        end
    end
    k=k+1;
end
%----------------------------fim modulação QPSK---------------------------%

Eb_N0_dB = 0:1:9; %Eb/N0 range
Eb_N0_lin = 10 .^ (Eb_N0_dB/10); %linearized Eb/N0 range
ber = zeros(size(Eb_N0_lin)); %pre-allocates BER vector
Eb = 2; % the average energy per bit, since we have 1 bit per symbol and average symbol energy = 1  (( 1^2 + (-1)^2 ) / 2)

NP = Eb ./ (2*Eb_N0_lin); %real-valued noise power = N0/2
NA = sqrt(NP); %noise amplitude is the square root of the noise power
    
for i = 1:length(Eb_N0_lin)
    n = NA(i)*complex(randn(1, num_b/2), randn(1, num_b/2))*sqrt(0.5); %noise vector with proper standard deviation
    r = mod + n; % received vector
    %-------------------------demodulação do sinal------------------------%
    for k=1:length(r) %informação está no signal da parte imaginaria e real do signal distorcido recebido
        if(real(r(k))>0) 
            if(imag(r(k))>0) %simbolo recebido real positivo imaginario positivo mapeia para bits 00
                demod(2*k-1) = 0; 
                demod(2*k) = 0;
            else            %simbolo recebido real positivo imaginario negativo mapeia para bits 11
                demod(2*k-1) = 1;
                demod(2*k) = 1;
            end
        else
            if(imag(r(k))>0) %simbolo recebido real negativo imaginario positivo mapeia para bits 01
                demod(2*k-1) = 0;
                demod(2*k) = 1;
            else             %simbolo recebido real negativo imaginario negativo mapeia para bits 10
                demod(2*k-1) = 1;
                demod(2*k) = 0;
            end
        end
    end
    %------------------fim demodulação do sinal---------------------------%
    ber(i) = sum(bits ~= demod) / num_b; % counts errors and calculates BER
end

display(ber);

ber_theoretical = 0.5*erfc(sqrt(2*Eb_N0_lin)/sqrt(2)); %theoretical BER for comparison

semilogy(Eb_N0_dB, ber, 'x', Eb_N0_dB, ber_theoretical, 'r', 'LineWidth', 2, 'MarkerSize', 10);