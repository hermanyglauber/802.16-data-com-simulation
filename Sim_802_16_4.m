clear;
close;
num_b = 100000; %number of bits to be simulated
k = 1;
bits=randi(2, 1, num_b)-1; %bits is a random real vector = 0 or +1 

%WirelessMAN-SC?

%802.16 bits encode strategy? 
%codigos concatenados Codigo externo + codigo interno: Reed-Solomon (com t ajustável), paridade, código convolucional e Block Turbo Code

N = 255;  % Codeword length
K = 223;  % Message length

txData = randi([0 1],223*16*300,1);
rsEncoder = comm.RSEncoder(N,K,'BitInput',true);
encData = rsEncoder(txData);
rsDecoder = comm.RSDecoder(N,K,'BitInput',true);
%demod_decod=0;
%display(length(txData));
%display(length(encData));
%802.16 signal modulation strategy 16-QAM
%constelaçao 16-QAM complexos
const(1)=complex(-3,3);
const(2)=complex(-3,1);
const(3)=complex(-3,-3);
const(4)=complex(-3,-1);
const(5)=complex(-1,3);
const(6)=complex(-1,1);
const(7)=complex(-1,-3);
const(8)=complex(-1,-1);
const(9)=complex(3,3);
const(10)=complex(3,1);
const(11)=complex(3,-3);
const(12)=complex(3,-1);
const(13)=complex(1,3);
const(14)=complex(1,1);
const(15)=complex(1,-3);
const(16)=complex(1,-1);

%disp(const);
ctrl=1;
indice=0;

%----------------------------inicio modulação 16-QAM------------------------%
for i=1:length(encData)/4
    %calcula o indice de 4 bits agrupados para mapear com o vetor
    %constelação 16-QAM gerando vetor de numeros complexos
    for j=1:4 
        indice=indice+encData(ctrl)*(2^(4-j));
        ctrl=ctrl+1;
    end
    %sinal modulado recebe simbol da matriz constelação 16-QAM que
    %representa 4 bits
    mod(i)=const(indice+1);
    indice=0;
end

%----------------------------fim modulação 16-QAM---------------------------%

Eb_N0_dB = 0:1:10; %Eb/N0 range
Eb_N0_lin = 10 .^ (Eb_N0_dB/10); %linearized Eb/N0 range
ber = zeros(size(Eb_N0_lin)); %pre-allocates BER vector
Eb = 4; % the average energy per bit, since we have 2 bit per symbol and average symbol energy = 2  (( 1^2 + 1^2 ) in each quadrant)

NP = Eb ./ (2*Eb_N0_lin); %real-valued noise power = N0/2
NA = sqrt(NP); %noise amplitude is the square root of the noise power

dist1=0;
dist2=0;
dist3=0;
dist4=0;
for i = 1:length(Eb_N0_lin)
    n = NA(i)*complex(randn(1, length(encData)/4), randn(1, length(encData)/4))*sqrt(0.5); %noise vector with proper standard deviation
    r = mod + n; % received vector
    %-------------------------demodulação do sinal------------------------%   
    for k=1:length(r) %informação está no signal da parte imaginaria e real do signal distorcido recebido
        if(real(r(k))>0) 
            if(imag(r(k))>0) %simbolo recebido real positivo imaginario positivo
                %calculo da distancia da amostra recebida e os 4 pontos do
                %quadrante real positivo imaginario positivo da constelação
                %16-QAM
                dist1=sqrt( (real(r(k)) - 1)^2 + (imag(r(k)) - 3) ^ 2 );
                dist2=sqrt( (real(r(k)) - 3)^2 + (imag(r(k)) - 3) ^ 2 );
                dist3=sqrt( (real(r(k)) - 1)^2 + (imag(r(k)) - 1) ^ 2 );
                dist4=sqrt( (real(r(k)) - 3)^2 + (imag(r(k)) - 1) ^ 2 );
                %caso a menor distancia seja em relação ao ponto (1,3)
                %então o símbolo recebido é '1100'
                if(dist1<dist2 && dist1<dist3 && dist1<dist4)                
                    demod(4*k-3) = 1; 
                    demod(4*k-2) = 1; 
                    demod(4*k-1) = 0; 
                    demod(4*k) = 0;
                end
                %caso a menor distancia seja em relação ao ponto (3,3)
                %então o símbolo recebido é '1000'
                if(dist2<dist1 && dist2<dist3 && dist2<dist4)                
                    demod(4*k-3) = 1; 
                    demod(4*k-2) = 0; 
                    demod(4*k-1) = 0; 
                    demod(4*k) = 0;
                end
                %caso a menor distancia seja em relação ao ponto (1,1)
                %então o símbolo recebido é '1101'
                if(dist3<dist1 && dist3<dist1 && dist3<dist4)                
                    demod(4*k-3) = 1; 
                    demod(4*k-2) = 1; 
                    demod(4*k-1) = 0; 
                    demod(4*k) = 1;
                end
                %caso a menor distancia seja em relação ao ponto (3,1)
                %então o símbolo recebido é '1001'
                if(dist4<=dist1 && dist4<=dist2 && dist4<=dist3)                
                    demod(4*k-3) = 1; 
                    demod(4*k-2) = 0; 
                    demod(4*k-1) = 0; 
                    demod(4*k) = 1;
                end
                    
            else            %simbolo recebido real positivo imaginario negativo
                %calculo da distancia da amostra recebida e os 4 pontos do
                %quadrante real positivo imaginario negativo da constelação
                %16-QAM
                dist1=sqrt( (real(r(k)) - 1)^2 + (imag(r(k)) + 1) ^ 2 );
                dist2=sqrt( (real(r(k)) - 3)^2 + (imag(r(k)) + 1) ^ 2 );
                dist3=sqrt( (real(r(k)) - 1)^2 + (imag(r(k)) + 3) ^ 2 );
                dist4=sqrt( (real(r(k)) - 3)^2 + (imag(r(k)) + 3) ^ 2 );
                %caso a menor distancia seja em relação ao ponto (1,-1)
                %então o símbolo recebido é '1111'
                if(dist1<dist2 && dist1<dist3 && dist1<dist4)                
                    demod(4*k-3) = 1; 
                    demod(4*k-2) = 1; 
                    demod(4*k-1) = 1; 
                    demod(4*k) = 1;
                end
                %caso a menor distancia seja em relação ao ponto (3,-1)
                %então o símbolo recebido é '1011'
                if(dist2<dist1 && dist2<dist3 && dist2<dist4)                
                    demod(4*k-3) = 1; 
                    demod(4*k-2) = 0; 
                    demod(4*k-1) = 1; 
                    demod(4*k) = 1;
                end
                %caso a menor distancia seja em relação ao ponto (1,-3)
                %então o símbolo recebido é '1110'
                if(dist3<dist1 && dist3<dist1 && dist3<dist4)                
                    demod(4*k-3) = 1; 
                    demod(4*k-2) = 1; 
                    demod(4*k-1) = 1; 
                    demod(4*k) = 0;
                end
                %caso a menor distancia seja em relação ao ponto (3,-3)
                %então o símbolo recebido é '1010'
                if(dist4<=dist1 && dist4<=dist2 && dist4<=dist3)                
                    demod(4*k-3) = 1; 
                    demod(4*k-2) = 0; 
                    demod(4*k-1) = 1; 
                    demod(4*k) = 0;
                end
            end
        else
            if(imag(r(k))>0) %simbolo recebido real negativo imaginario positivo
                %calculo da distancia da amostra recebida e os 4 pontos do
                %quadrante real negativo imaginario positivo da constelação
                %16-QAM
                dist1=sqrt( (real(r(k)) + 3)^2 + (imag(r(k)) - 3) ^ 2 );
                dist2=sqrt( (real(r(k)) + 1)^2 + (imag(r(k)) - 3) ^ 2 );
                dist3=sqrt( (real(r(k)) + 3)^2 + (imag(r(k)) - 1) ^ 2 );
                dist4=sqrt( (real(r(k)) + 1)^2 + (imag(r(k)) - 1) ^ 2 );
                %caso a menor distancia seja em relação ao ponto (-3,3)
                %então o símbolo recebido é '0000'
                if(dist1<dist2 && dist1<dist3 && dist1<dist4)                
                    demod(4*k-3) = 0; 
                    demod(4*k-2) = 0; 
                    demod(4*k-1) = 0; 
                    demod(4*k) = 0;
                end
                %caso a menor distancia seja em relação ao ponto (-1,3)
                %então o símbolo recebido é '0100'
                if(dist2<dist1 && dist2<dist3 && dist2<dist4)                
                    demod(4*k-3) = 0; 
                    demod(4*k-2) = 1; 
                    demod(4*k-1) = 0; 
                    demod(4*k) = 0;
                end
                %caso a menor distancia seja em relação ao ponto (-3,1)
                %então o símbolo recebido é '0001'
                if(dist3<dist1 && dist3<dist1 && dist3<dist4)                
                    demod(4*k-3) = 0; 
                    demod(4*k-2) = 0; 
                    demod(4*k-1) = 0; 
                    demod(4*k) = 1;
                end
                %caso a menor distancia seja em relação ao ponto (-1,3)
                %então o símbolo recebido é '0101'
                if(dist4<=dist1 && dist4<=dist2 && dist4<=dist3)                
                    demod(4*k-3) = 0; 
                    demod(4*k-2) = 1; 
                    demod(4*k-1) = 0; 
                    demod(4*k) = 1;
                end
            else             %simbolo recebido real negativo imaginario negativo
                  %calculo da distancia da amostra recebida e os 4 pontos do
                %quadrante real negativo imaginario negativo da constelação
                %16-QAM
                dist1=sqrt( (real(r(k)) + 3)^2 + (imag(r(k)) + 1) ^ 2 );
                dist2=sqrt( (real(r(k)) + 1)^2 + (imag(r(k)) + 1) ^ 2 );
                dist3=sqrt( (real(r(k)) + 3)^2 + (imag(r(k)) + 3) ^ 2 );
                dist4=sqrt( (real(r(k)) + 1)^2 + (imag(r(k)) + 3) ^ 2 );
                %caso a menor distancia seja em relação ao ponto (-3,-1)
                %então o símbolo recebido é '0011'
                if(dist1<dist2 && dist1<dist3 && dist1<dist4)                
                    demod(4*k-3) = 0; 
                    demod(4*k-2) = 0; 
                    demod(4*k-1) = 1; 
                    demod(4*k) = 1;
                end
                %caso a menor distancia seja em relação ao ponto (-1,-1)
                %então o símbolo recebido é '0111'
                if(dist2<dist1 && dist2<dist3 && dist2<dist4)                
                    demod(4*k-3) = 0; 
                    demod(4*k-2) = 1; 
                    demod(4*k-1) = 1; 
                    demod(4*k) = 1;
                end
                %caso a menor distancia seja em relação ao ponto (-3,-3)
                %então o símbolo recebido é '0010'
                if(dist3<dist1 && dist3<dist1 && dist3<dist4)                
                    demod(4*k-3) = 0; 
                    demod(4*k-2) = 0; 
                    demod(4*k-1) = 1; 
                    demod(4*k) = 0;
                end
                %caso a menor distancia seja em relação ao ponto (-1,-3)
                %então o símbolo recebido é '0110'
                if(dist4<=dist1 && dist4<=dist2 && dist4<=dist3)                
                    demod(4*k-3) = 0; 
                    demod(4*k-2) = 1; 
                    demod(4*k-1) = 1; 
                    demod(4*k) = 0;
                end
            end
        end
    end
    %------------------fim demodulação do sinal---------------------------%
    %display(length(demod));
    demod_decod = rsDecoder(demod.');
    ber(i) = sum(txData ~= demod_decod) / num_b; % counts errors and calculates BER
    display(ber(i));
end

ber_theoretical = 0.5*erfc(sqrt(2*Eb_N0_lin)/sqrt(2)); %theoretical BER for comparison

semilogy(Eb_N0_dB, ber, 'x', Eb_N0_dB, ber_theoretical, 'r', 'LineWidth', 2, 'MarkerSize', 10);