function [ final_uni,uni_hist,hist_time ] = unitaryEvolutionTotal_MATLAB( hammy,opt_params )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%uni = zeros(512,1);

init_uni = eye(16); % U(t0,t0) = 1
unitary = zeros(256,1);

% Code for doing conversion between 16x16 matrices and a 256x1 array. This
% kind of code pops up in several of these files. Just look for moda and
% diva as indicators that it's doing a conversion.
for a = 1:1:256
    if mod(a,16) >= 1
        moda = mod(a,16);
    else
        moda = 16;
    end
    diva = floor((a/16)-0.001)+1;
    unitary(a) = init_uni(moda,diva);
end

% init_deriv = hammy(1e-20);
% for a = 257:1:512
%     if mod(a,16) >= 1
%         moda = mod(a,16);
%     else
%         moda = 16;
%     end
%     diva = floor(((a-256)/16)-0.001)+1;
%     unitary(a) = init_deriv(moda,diva);
% end


eqs = @(t,uni)equationMatrixMaker(hammy,uni,t);

% coupled_eqs = @(t,uni) functionator_MATLAB(eqs,256,1,t,uni);

% Use ode45 for non-stiff equations, and ode15s for stiff equations
options = odeset('RelTol',1e-5,'AbsTol',1e-6);
[t,unitary_hist] = ode45(eqs,[1e-20 opt_params.tot_time-1e-18],unitary,options);
% [t,unitary_hist] = ode45(coupled_eqs,[3.99999999999e-5 opt_params.tot_time],unitary);

final_uni = zeros(16,16);

for a = 1:1:256
    if mod(a,16) >= 1
        moda = mod(a,16);
    else
        moda = 16;
    end
    diva = floor((a/16)-0.001)+1;
    final_uni(moda,diva) = unitary_hist(length(t),a);
end

uni_hist = zeros(16,16,length(t));

for a = 1:1:length(t)
    for b = 1:1:256
        if mod(b,16) >= 1
           moda = mod(b,16);
        else
           moda = 16;
        end
        diva = floor((b/16)-0.001)+1;
        uni_hist(moda,diva,a) = unitary_hist(a,b);   
    end
end

hist_time = t;

end

