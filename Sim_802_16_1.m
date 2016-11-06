clear;
close;
num_b = 1000; %number of bits to be simulated
bits = randi(2, 1, num_b)-1; %bits is a random real vector = 0 or +1 

%802.16 bits encode strategy? 

%802.16 signal modulation strategy?

display(bits);

Eb_N0_dB = 0:1:9; %Eb/N0 range
Eb_N0_lin = 10 .^ (Eb_N0_dB/10); %linearized Eb/N0 range
ber = zeros(size(Eb_N0_lin)); %pre-allocates BER vector
Eb = 1; % the average energy per bit, since we have 1 bit per symbol and average symbol energy = 1  (( 1^2 + (-1)^2 ) / 2)

NP = Eb ./ (2*Eb_N0_lin); %real-valued noise power = N0/2
NA = sqrt(NP); %noise amplitude is the square root of the noise power
    
for i = 1:length(Eb_N0_lin)
    n = NA(i)*randn(1, num_b); %noise vector with proper standard deviation
    r = bits + n; % received vector
    demod = sign(r); % the information is in the sign
    ber(i) = sum(bits ~= demod) / num_b; % counts errors and calculates BER
end

ber_theoretical = 0.5*erfc(sqrt(2*Eb_N0_lin)/sqrt(2)); %theoretical BER for comparison
semilogy(Eb_N0_dB, ber, 'x', Eb_N0_dB, ber_theoretical, 'r', 'LineWidth', 2, 'MarkerSize', 10);