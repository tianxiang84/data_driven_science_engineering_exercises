%% 2D demonstration of the SVD
clear; clc;

%% Raw data
x=0:0.1:100;     %[meter] x coord vs time
y=4*(0:0.1:100); %[meter] y coord vs time
y = x + randn(1,length(x)); %[meter] put noise in the x coordinate
y = y + randn(1,length(y)); %[meter] put noise in the y coordinate
Traj = [x;y]; %[meter;meter] noisy 2D trajectory
plot(Traj(1,5:10:end),Traj(2,5:10:end),'bo-'); hold on;
xlabel('x [meter]','FontSize',16);
ylabel('y [meter]','FontSize',16);
set(gca,'FontSize',16);
grid on;

%% SVD dimension reduction
[U,S,V]=svd(Traj); % Do a single value decomposition
rank = 1;
U = U(:,1:rank); % Get the first component only
S = S(1:rank,1:rank); 
V = V';
V = V(1:rank,:);
newTraj = U*S*V; %[meter;meter] coordinates of the reduced dimension trajectory
plot(newTraj(1,1:10:end),newTraj(2,1:10:end),'r-','LineWidth',2.0); hold on;
xlabel('x [meter]','FontSize',16);
ylabel('y [meter]','FontSize',16);
set(gca,'FontSize',16);
grid on;

figure
plot(V);