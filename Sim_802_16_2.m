clear;
close;
num_b = 100000; %number of bits to be simulated
k = 1;
bits=randi(2, 1, num_b)-1; %bits is a random real vector = 0 or +1 

%WirelessMAN-SC?

%802.16 bits encode strategy? 
%codigos concatenados Codigo externo + codigo interno: Reed-Solomon (com t ajustável), paridade, código convolucional e Block Turbo Code


%%%% REDSOLOMON ENCODER %%%% BEGIN
function [ encoded ] = rsEncoder( msg, m, prim_poly, n, k )
    %RSENCODER Encode message with the Reed-Solomon algorithm
    % m is the number of bits per symbol
    % prim_poly: Primitive polynomial p(x). Ie for DM is 301
    % k is the size of the message
    % n is the total size (k+redundant)
    % Example: msg = uint8('Test')
    % enc_msg = rsEncoder(msg, 8, 301, 12, numel(msg));
    
    % Get the alpha
    alpha = gf(2, m, prim_poly);
    
    % Get the reed-solomon generating polynomial g(x)
    g_x = genpoly(k, n, alpha);
    
    % Multiply the information by X^(n-k), or just pad with zeros at the end to
    % get space to add the redundant information
    msg_padded = gf([msg zeros(1, n-k)], m, prim_poly);
    
    % Get the remainder of the division of the extended message by the 
    % reed-solomon generating polynomial g(x)
    [~, reminder] = deconv(msg_padded, g_x);

    % Now return the message with the redundant information
    encoded = msg_padded - reminder;

end

% Find the reed-solomon generating polynomial g(x), by the way this is the
% same as the rsgenpoly function on matlab
function g = genpoly(k, n, alpha)
    g = 1;
    % A multiplication on the galois field is just a convolution
    for k = mod(1 : n-k, n)
        g = conv(g, [1 alpha .^ (k)]);
    end
end


%%%% REDSOLOMON ENCODER %%%% END


%%%% REDSOLOMON DECODER %%%% BEGIN

function [ decoded, error_pos, error_mag, g, S ] = rsDecoder( encoded, m, prim_poly, n, k )
    %RSDECODER Decode a reed-solomon encoded message
    %   Example:
    % [dec, ~, ~, ~, ~] = rsDecoder(enc_msg, 8, 301, 12, numel(msg))
    max_errors = floor((n-k)/2);
    orig_vals = encoded.x;
    % Initialize the error vector
    errors = zeros(1, n);
    g = [];
    S = [];
    
    % Get the alpha
    alpha = gf(2, m, prim_poly);
    
    % Find the syndromes (Check if dividing the message by the generator
    % polynomial the result is zero)
    Synd = polyval(encoded, alpha .^ (1:n-k));
    Syndromes = trim(Synd);
    
    % If all syndromes are zeros (perfectly divisible) there are no errors
    if isempty(Syndromes.x)
        decoded = orig_vals(1:k);
        error_pos = [];
        error_mag = [];
        g = [];
        S = Synd;
        return;
    end
    
    % Prepare for the euclidean algorithm (Used to find the error locating
    % polynomials)
    r0 = [1, zeros(1, 2*max_errors)]; r0 = gf(r0, m, prim_poly); r0 = trim(r0);
    size_r0 = length(r0);
    r1 = Syndromes;
    f0 = gf([zeros(1, size_r0-1) 1], m, prim_poly);
    f1 = gf(zeros(1, size_r0), m, prim_poly);
    g0 = f1; g1 = f0;
    
    % Do the euclidian algorithm on the polynomials r0(x) and Syndromes(x) in
    % order to find the error locating polynomial
    while true
        % Do a long division
        [quotient, remainder] = deconv(r0, r1);    
        % Add some zeros
        quotient = pad(quotient, length(g1));
        
        % Find quotient*g1 and pad
        c = conv(quotient, g1);
        c = trim(c);
        c = pad(c, length(g0));
        
        % Update g as g0-quotient*g1
        g = g0 - c;
        
        % Check if the degree of remainder(x) is less than max_errors
        if all(remainder(1:end - max_errors) == 0)
            break;
        end
        
        % Update r0, r1, g0, g1 and remove leading zeros
        r0 = trim(r1); r1 = trim(remainder);
        g0 = g1; g1 = g;
    end
    
    % Remove leading zeros
    g = trim(g);
    
    % Find the zeros of the error polynomial on this galois field
    evalPoly = polyval(g, alpha .^ (n-1 : -1 : 0));
    error_pos = gf(find(evalPoly == 0), m);
    
    % If no error position is found we return the received work, because
    % basically is nothing that we could do and we return the received message
    if isempty(error_pos)
        decoded = orig_vals(1:k);
        error_mag = [];
        return;
    end
    
    % Prepare a linear system to solve the error polynomial and find the error
    % magnitudes
    size_error = length(error_pos);
    Syndrome_Vals = Syndromes.x;
    b(:, 1) = Syndrome_Vals(1:size_error);
    for idx = 1 : size_error
        e = alpha .^ (idx*(n-error_pos.x));
        err = e.x;
        er(idx, :) = err;
    end
    
    % Solve the linear system
    error_mag = (gf(er, m, prim_poly) \ gf(b, m, prim_poly))';
    % Put the error magnitude on the error vector
    errors(error_pos.x) = error_mag.x;
    % Bring this vector to the galois field
    errors_gf = gf(errors, m, prim_poly);
    
    % Now to fix the errors just add with the encoded code
    decoded_gf = encoded(1:k) + errors_gf(1:k);
    decoded = decoded_gf.x;
    
end

% Remove leading zeros from galois array
function gt = trim(g)
    gx = g.x;
    gt = gf(gx(find(gx, 1) : end), g.m, g.prim_poly);
end

% Add leading zeros
function xpad = pad(x,k)
    len = length(x);
    if (len<k)
        xpad = [zeros(1, k-len) x];
    end
end
%%%% REDSOLOMON DECODER %%%% END

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
for i=1:length(bits)/4
    %calcula o indice de 4 bits agrupados para mapear com o vetor
    %constelação 16-QAM gerando vetor de numeros complexos
    for j=1:4 
        indice=indice+bits(ctrl)*(2^(4-j));
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
    n = NA(i)*complex(randn(1, num_b/4), randn(1, num_b/4))*sqrt(0.5); %noise vector with proper standard deviation
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
    ber(i) = sum(bits ~= demod) / num_b; % counts errors and calculates BER
end

ber_theoretical = 0.5*erfc(sqrt(2*Eb_N0_lin)/sqrt(2)); %theoretical BER for comparison

semilogy(Eb_N0_dB, ber, 'x', Eb_N0_dB, ber_theoretical, 'r', 'LineWidth', 2, 'MarkerSize', 10);