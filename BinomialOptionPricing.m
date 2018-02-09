% HW2
% Arjun Singh

clear
clc

% Define all variables
mu = 0.05;
sigma = 0.2;
T = 1; % End time
r = 0.03; % Risk free interest rate
S0 = 5; % Share price at t = 0
m = 730; % Number of time steps
dt = T/m;

% Define u,v and p
u = 1 + sigma*sqrt(dt)*sqrt(1+(mu^2*dt)/sigma^2);
v = 1 - sigma*sqrt(dt)*sqrt(1+(mu^2*dt)/sigma^2);
p = 0.5 * (1 + mu*sqrt(dt)/(sigma*sqrt(1+(mu^2*dt)/sigma^2)));
type = 1; % 0 for call, put otherwise


strikePrices = 0:1:10; % Vector for strike prices
optionPrices = zeros(1,length(strikePrices)); % Vector to store option prices
meanExpiration = zeros(1,length(strikePrices));

for strike = 1:length(strikePrices)

    prices = zeros(m+1,m+1);
    prices(1,1) = S0; % Matrix to store asset prices
    
    for i = 2:m+1
        prices(1:i-1,i) = prices(1:i-1,i-1)*u;
        prices(i,i) = prices(i-1,i-1)*v;  % Multiplying asset prices by u and v
    end

    vals = zeros(size(prices)); % Matrix to store option prices
    if type == 0
        vals(:,end) = max(prices(:,end)-strikePrices(strike),0); % Find intrinsic value of option
    else
        vals(:,end) = max(strikePrices(strike)-prices(:,end),0);
        
    end
    
    prob = zeros(m+1,1);
    for i = 1:m+1
        prob(i,1) = (nchoosek(m+1,i))*p^i*(1-p)^(m+1-i); % Probability of last price
    end
    
    meanExpiration(strike) = sum(prob.*(vals(:,end))); % Expected value
    
    for i = m:-1:1
        vals(1:i,i) = exp(-r*dt)*(p*vals(1:i,i+1) + (1-p)*vals(2:i+1,i+1)); % Backtrack to find current option price
    end
    
    optionPrices(strike) = vals(1); % Store current value in vector
end

% Part a
plot(strikePrices, optionPrices)
xlabel('Strike Price')
ylabel('Value of European Option')
figure
% 
% Part b
profitPercent = (meanExpiration - optionPrices)./optionPrices;
plot(strikePrices, profitPercent)
xlabel('Strike Price')
ylabel('Profit Percentage')
