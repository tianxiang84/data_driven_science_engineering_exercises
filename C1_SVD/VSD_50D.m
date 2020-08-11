%% SVD for a 50-segment mass-spring system

%% Run a full simulation and store the full state as a 100 x T matrix
numSeg = 50; %[-] number of segments
dt = 0.01;   %[sec] time step
numStep = 600*ceil(1.0/dt); %[-] number of time step to run
dl = 1.0; %[meter] segment length
k = 3.0;
c = 1.0;
m = 1.0;

% Initial state
fullState = zeros(2*numSeg,numStep); % full state matrix
x = dl*((1:1:numSeg)-1); %[meter] initial position
v = 0*x;                 %[m/s] initial velocity
fullState(:,1) = [x';v'];

% Loop the dynamics
for i=1:1:numStep-1
    currT = i*dt;
    
    intF = zeros(numSeg+1,1);
    for j=2:numSeg
        intF(j) = k * (x(j)-x(j-1)-dl);
    end
    intF(1) = -1.0;
    intF(end) = 0.0;
    
    extF = zeros(numSeg,1);
    accCorr = zeros(numSeg,1);
    for j=1:numSeg
        extF(j) = -c*v(j);
        accCorr(j) = (-c)*dt*0.5/m;
    end
    
    for j=1:numSeg
        v(j) = v(j) + dt*((intF(j+1)-intF(j)+extF(j))/m)/(1.0-accCorr(j));
    end
    %dv/dt = f/m + (df/dv dv/dt dt 0.5)/m
    %(1 - *df/dv dt 0.5 / m*) dv/dt = f/m
    
    for j=1:numSeg
        x(j) = x(j) + dt*v(j);
    end
    
%     if mod(i,100)==0
%         plot(currT,x(1),'bo',currT,x(end),'ro'); hold on;
%     end
%     pause(0.001);
    
    fullState(:,i+1) = [x';v'];
end % end of dynamics time step loop


[U,S,V] = svd(fullState,'econ');
V = V';
rank = 4;
U = U(:,1:rank); % Get the first component only
S = S(1:rank,1:rank); 
V = V(1:rank,:);
approxState4 = U*S*V;

rank = 3;
U = U(:,1:rank); % Get the first component only
S = S(1:rank,1:rank); 
V = V(1:rank,:);
approxState3 = U*S*V;

rank = 2;
U = U(:,1:rank); % Get the first component only
S = S(1:rank,1:rank); 
V = V(1:rank,:);
approxState2 = U*S*V;

rank = 1;
U = U(:,1:rank); % Get the first component only
S = S(1:rank,1:rank); 
V = V(1:rank,:);
approxState1 = U*S*V;

% figure 
% for i=1:1:numStep
%     if mod(i,500)==0
%         plot(fullState(1:50,i),fullState(1:50,i)*0,'bo-'); hold on;
%         plot(approxState1(1:50,i),approxState1(1:50,i)*0+1,'ro-'); hold on;
%         plot(approxState2(1:50,i),approxState2(1:50,i)*0+2,'ro-'); hold on;
%         plot(approxState3(1:50,i),approxState3(1:50,i)*0+3,'ro-'); hold on;
%         plot(approxState4(1:50,i),approxState4(1:50,i)*0+4,'ro-'); hold on;
%         ylim([-0.5,4.5]);
%         xlim([-10,100]);
%         xlabel('x [meter]','FontSize',16);
%         ylabel('rank [-]','FontSize',16);
%         set(gca,'FontSize',16);
%         pause(0.1);
%     end
%     hold off;
% end


for i=1:1:numStep
    e1(i) = norm(fullState(:,i) - approxState1(:,i));
    e2(i) = norm(fullState(:,i) - approxState2(:,i));
    e3(i) = norm(fullState(:,i) - approxState3(:,i));
    e4(i) = norm(fullState(:,i) - approxState4(:,i));
end 
semilogy(e1,'k'); hold on;
semilogy(e2,'b'); hold on;
semilogy(e3,'r'); hold on;
semilogy(e4,'g'); hold on;
%% Perform dimension reduction and plot the results