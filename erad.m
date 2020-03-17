function [S,R,Z,t] = erad(Npop,birthrate,alpha,beta,zeta,delta,killRate,tRaid,tmax,dt)
% [S,R,Z,t] = erad(Npop,alpha,beta,zeta,delta,killRate,tRaid,tmax,dt)
% simulates the time-evolution of the human population, dead and zombie
% population after propagation of zombiism and elimination by raids. 
% This function is inspired by the original function created by Philip Munz [1]. 
% Credit goes for the computational model goes, therefore to [1]. 
% The simulation relies on a system of coupled ordinary differential equations.
% 
% Input
% 
% Npop: scalar: initial population
% birthrate: scalar: population birth rate
% alpha : scalar: "zombie destruction" rate.
% beta: scalar:  transformation into a zombie rate.
% zeta: scalar:  zombie resurrection rate.
% delta: scalar:  "natural population death" rate.
% killRate: scalar: kill rate.
% tRaid: scalar: vector [1xNraid]: time at which raid are conducted
% tmax: scalar: max duration of the simulation
% dt: scalar: time step
% 
% Output
% 
% S: vector [1xN]: number of susceptibles ( human population)
% R: vector [1xN]:  number of dead
% Z: vector [1xN]:  number of zombies
% 
% Reference:
% Munz, P., Hudea, I., Imad, J., & Smith, R. J. (2009).
% When zombies attack!: mathematical modelling of an outbreak of zombie infection.
% Infectious Disease Modelling Research Progress, 4, 133-150.
%
%
% Author: E. Cheynet -UiB - last modified: 16-03-2020
% 
% see also zombies.m

%%

N = round(1+tmax/dt); % total number of time steps
Y = zeros(3,N);
Y(:,1) = [Npop;0;0]; % SIZR

t = (0:1:N-1)*dt; % time
modelFun = @(Y,A,F) A*Y + F;


for ii=1:N-1
    A = getA(zeta,delta);
    SZ = Y(1,ii)*Y(2,ii);
    F = [birthrate-beta;beta-alpha;alpha].*SZ;
    Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
    
    if min(abs(t(ii)-tRaid))<= 0.1
        Y(2,ii+1) =  max(0,Y(2,ii) -killRate*Y(2,ii));
        
        [ ~,ind] = (min(abs(t(ii)-tRaid)));
        tRaid(ind)=[];
    end
    
end

t = (0:1:N-1)*dt; % time
S = Y(1,1:N);
Z = Y(2,1:N);
R = Y(3,1:N);


    function [A] = getA(zeta,delta)
        S0 = [-delta,0,0];
        Z0 = [0,0,zeta];
        R0 = [delta,0,-zeta,];
        A = [S0;Z0;R0];
    end
    function [Y] = RK4(Fun,Y,A,F,dt)
        
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A,F);
        k_2 = Fun(Y+0.5*dt*k_1,A,F);
        k_3 = Fun(Y+0.5*dt*k_2,A,F);
        k_4 = Fun(Y+k_3*dt,A,F);
        % output
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    end



end


