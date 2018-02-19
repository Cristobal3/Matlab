


function Final

% solves the differential equation xdotdot + 0.05*xdot + x^3 = 7.5*sin(t)
% this is the Duffing equation
% for certain parameters it is a chaotic system
% theta=pi/2:0.01:0;

% M=[1/3 0; 0 3/2+cos(t)+(1/3)*(cos(t))^2];
% F=[0.5*cos(theta)-(0.5+(1/3)*cos(theta))*sin(theta)*(q(3))^2;
%     (1+(1/3)*cos(theta))*sin(theta)*q(3)*q(4)-4)];
calltoODE = @(t,q)Finalode(t,q); % request to solve the ODE
options = odeset('RelTol',1e-5,'AbsTol',1e-10); % solver options
timespan = [0 5]; % length of time to integrate equation
IC = [0 0 0 0 0 0 0 0]; % initial conditions for x, dotx

% request to solve ODE
% time contains column vector of solution times
% solution contains column vectors of x, dotx
[time,solution]=ode45(calltoODE,timespan,IC,options); % request to solve ODE

% plot x, xdot as functions of time
figure
plot(time,solution(:,2))
ylabel('Desired variable')
xlabel('Time (s)')
hold on
plot(time, solution(:,3))
plot(time, solution(:,4))
plot(time, solution(:,1))
legend('Phi','Psi','theta1','theta')

figure
plot(time,solution(:,6))
ylabel('Desired variable')
xlabel('Time (s)')
hold on
plot(time, solution(:,7))
plot(time, solution(:,8))
plot(time, solution(:,5))
legend('Phidot','Psidot','theta1dot','thetadot')

%% Animation
L1 = .62;
L2 = .66;

figure
for t=1:length(time)
    clf
    z = [L1*cos(solution(t,2)) 0];
    x = [L1*sin(solution(t,2)) 0];
    y = [L1*sin(solution(t,2))*sin(solution(t,1)) 0];
    %abs(L2*sin(pi+solution(t,3))+L1*cos(solution(t,2)))
    %abs(L1*sin(solution(t,2))+L2*cos(pi + solution(t,3)))
    z1 = [L1*cos(solution(t,2)) L2*sin(solution(t,3))+L1*cos(solution(t,2))];
    x1 = [L1*sin(solution(t,2)) L1*sin(solution(t,2))+L2*cos(solution(t,3))];
    y1 = [L1*sin(solution(t,2))*sin(solution(t,1)) L1*sin(solution(t,2))*...
        sin(solution(t,1))+ L2*sin(solution(t,1))];
    plot3(x,y,z)
    axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
    %view(0,90)
    xlabel('X label')
    ylabel('Y label')
    zlabel('Z Label')
    hold on
    plot3(x1,y1,z1)
    hold off
    pause(0.01)
end


%% example ODE
    function qdot = Finalode(t,q)
        
        %     % this is the first-order form of the original ODE
        %     q(1) = theta;
        %     q(2) = phi;
        %     q(3) = psi;
        %     q(4) = theta1;
        %     q(5) = thetadot;
        %     q(6) = phidot;
        %     q(7) = psidot
        %     q(8) = theta1dot;
        %     % therefore
        
        %%
        %Constants
        m1 = 1;
        m2 = 1;
        m3 = 1;
        I1 = 1;
        I2 = 1;
        I3 = 1;
        L1 = 1;
        L2 = 1;
        R1 = 1;
        R2 = 1;
        % in vector form
        M1 = [0.5*m1*R1^2 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
        F1 =[-2;0;0;0];
        M2=[I2 I2*cos(q(2)) 0  0; I2*cos(q(2)) I2+m2*(L1^2/4) 0 0; 0 0 0 0; 0 0 0 0 ];
        F2=[-I2*q(6)^2*sin(q(2));
            -0.5*m2*(9.81*(L1/2)*sin(q(2))); 0; 0];
        M3 = [ I3 0 0 0;  0 I3+m3*L1^2 I3+m3*sin(q(2)-q(3)) I3; 0 I3+sin(q(2)-q(3)) I3+m3+(L2^2/4) I3;....
            0 I3 I3 I3];
        F3 = [0;-m3*((q(6)*q(6)-q(7)^2)*cos(q(2)-q(3)))+0.5*m3*(q(6)*q(7)*cos(q(2)-q(3))-9.81*L1*sin(q(2)));...
            -m3*(q(6)*cos(q(2)-q(3))*(q(6)-q(7)))+m3*(-9.81*(L2/2)*cos(q(3))-q(6)*q(7)*cos(q(2)-q(3)))+L2/2*9.81*m3; 0];
        M = M1 + M2 + M3;
        F = F1 + F2 + F3;
        qdot = [ q(5);
            q(6);
            q(7);
            q(8)
            M\F];
        
        
        
    end

end
