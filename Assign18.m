%  Name                   	Student id                email
% +------------------------+------------------------+-------------------------
% |                        |                        |
% +------------------------+------------------------+-------------------------
% |                        |                        |
% +------------------------+------------------------+-------------------------
% I (enlisted above) declare that:
%   1. Our assignment will be our own work.
%   2. We shall not make solutions to the assignment available to anyone else.
%   3. We shall not engage in any other activities that will dishonestly improve my results or dishonestly improve or hurt the results of others.
clear; clc; tic;       % clear memory & screen
rng(1429);             % reset seed
n=50;                  fprintf('n=%i\n',n);
REP =1000;             fprintf('REP=%i\n',REP);
BOOTREP=499;           fprintf('BOOTREP=%i\n',BOOTREP);
m = 0;
s = 1.3;
mu = exp(m+1/2*s^2);   fprintf('mu=%8.4f\n',mu);
beta = sin(mu);        fprintf('beta=%8.4f\nReplications:',beta);

xbar =zeros(REP,1);    % average original sample
bhat =zeros(REP,1);    % estimate of beta
SE   =zeros(REP,1);    % standard error bhat (asymptotic)
trat  =zeros(REP,1);   % t-ratio
LCLasym =zeros(REP,1); % Lower confidence limit (asym)
UCLasym =zeros(REP,1); % Upper confidence limit (asym)

for i=1:REP
  if mod(i,100)==0; fprintf('%i ',i); end
  X=exp(random('Normal',m,s,[n,1]));
  xbar(i)=mean(X);
  bhat(i)=sin(xbar(i));
  SE(i)=0.14;  % CODE THE SE YOURSELF!!!
  trat(i)=(bhat(i)-beta)/SE(i);
  LCLasym(i) =bhat(i)-1.96*SE(i);
  UCLasym(i) =bhat(i)+1.96*SE(i);
end; fprintf('\n');

CoverageFreqasym=mean((beta>LCLasym) & (beta<UCLasym));
fprintf('Coverage freq.(asym):    %6.3f\n',CoverageFreqasym);
toc