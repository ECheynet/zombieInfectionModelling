function [S,I,Z,R,Q,t] = zombies(Npop,alpha,beta,zeta,delta,rho,kappa,sigma,gamma,c,birthrate,tmax,dt)
%  [S,I,Z,R,Q,t] = zombies(Npop,alpha,beta,zeta,delta,rho,kappa,sigma,gamma,c,birthrate,tmax,dt)
% simulates the time-evolution of the human population, dead and zombie
% population after propagation of zombiism. This function is inspired by the original function created by Philip Munz [1]. Credit goes for the
% computational model goes, therefore to [1]. The simulation relies on a
% system of coupled ordinary differential equations.
%
% Input
%
% Npop: scalar: initial population
% alpha : scalar: "zombie destruction" rate.
% beta: scalar:  transformation into a zombie rate.
% zeta: scalar:  zombie resurrection rate.
% delta: scalar:  "natural population death" rate.
% rho: scalar: proportional to incumbation time
% kappa: scalar: Rate of infected human being quarantined
% sigma: scalar:  Rate of zombie being quarantined
% gamma: scalar:  rate of quanrantined person killed while trying to escape the quarantined zone
% c: scalar:  rate of recovery (a cure exists for infected humans who did not turn into zombies!)
% birthrate: scalar:the birth rate of the population 
% tmax: scalar: max duration of the simulation
% dt: scalar: time step
%
%
% Output
%
% S: vector [1xN] Non-infected human population outisde the quarantine zone
% I: vector [1xN] Infected human population outisde the quarantine zone
% R: vector [1xN] Number of deads
% Z: vector [1xN] Zombie population outisde the quarantine zone
% Q:  vector [1xN] Quarantined population (zombies + humans)
% t: vector [1xN] non-dimensional time
% 
% Reference:
% Munz, P., Hudea, I., Imad, J., & Smith, R. J. (2009).
% When zombies attack!: mathematical modelling of an outbreak of zombie infection.
% Infectious Disease Modelling Research Progress, 4, 133-150.
%
%
% Author: E. Cheynet -UiB - last modified: 16-03-2020
% 
% see also erad.m

%% Initial conditions
N = round(1+tmax/dt);
Y = zeros(5,N);
Y(:,1) = [Npop;0;0;0;0]; % SIZR




%%

modelFun = @(Y,A,F) A*Y + F;
% ODE reYution
for ii=1:N-1
    
    
    A = getA(zeta,delta,rho,kappa,sigma,gamma,c);
    SZ = Y(1,ii)*Y(3,ii);
    F = [birthrate-beta;beta;-alpha; alpha;0].*SZ;
    Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
    
    if round(diff(Y(3,ii:ii-1)))==0
        t = (0:1:ii-1)*dt; % time
        S = Y(1,1:ii);
        I = Y(2,1:ii);
        Z = Y(3,1:ii);
        R = Y(4,1:ii);
        Q = Y(5,1:ii);
        
        fprintf('Simulation finished \n')
        return
    end
end

t = (0:1:N-1)*dt; % time
S = Y(1,1:N);
I = Y(2,1:N);
Z = Y(3,1:N);
R = Y(4,1:N);
Q = Y(5,1:N);

    function [A] = getA(zeta,delta,rho,kappa,sigma,gamma,c)
        S1 = [-delta,c,0,0,0];
        I2 = [0,-rho-delta-kappa-c,0,0,0];
        Z3 = [0,rho,-sigma,zeta,0];
        R4 = [delta,delta,0,-zeta,gamma];
        Q5 = [0,kappa,sigma,0,-gamma];
        A = [S1;I2;Z3;R4;Q5];
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


