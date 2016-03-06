%% Created BY: Prashant Shekhar
% Emulator for Eruption height and Partile flux
% Credits: MPerk Library

%% Getting ready
clc; 
clear;

%% Reading the data
x = csvread('lhs_samples.csv');
y1 = csvread('eruption_height.csv');
y2 = csvread('particle_flux.csv');

y = horzcat(y1,y2);


%% generating points for sparse representation
% for first output

for i = 1:16
    xpred(i,1) = min(x(:,1)) + rand(1)*(max(x(:,1))-min(x(:,1)));
    xpred(i,2) = min(x(:,2)) + rand(1)*(max(x(:,2))-min(x(:,2)));
end


F = [ones(64,1) x(:,1) x(:,2) diag(x(:,1)*(x(:,2))')];   
        
Fpred = [ones(16,1) xpred(:,1) xpred(:,2) diag(xpred(:,1)*(xpred(:,2))')];

ypred1 = mperk('X',x, ...
                'Y',y(:,1),...
                'CrossValidation','Yes',...
                'CorrelationEstimationMethod','REML',...
                'CorrelationFamily','Gaussian',...
                'Xpred',xpred,...
                'RegressionModel',F,...
                'PredRegressionModel', Fpred);

% for second output            
                

ypred2 = mperk('X',x, ...
                'Y',y(:,2),...
                'CrossValidation','Yes',...
                'CorrelationEstimationMethod','REML',...
                'CorrelationFamily','Gaussian',...
                'Xpred',xpred,...
                'RegressionModel',F,...
                'PredRegressionModel', Fpred);
            

            
Ypred = horzcat(ypred1.preds.ypreds,ypred2.preds.ypreds);


%% Visualizing the output
%Output 1: Eruption Height

figure
tri = delaunay(x(:,1),x(:,2));
trisurf(tri,x(:,1),x(:,2),y(:,1));

% Clean it up
l = light('Position',[-50 -15 29]);
lighting none
shading interp
colorbar EastOutside

hold on
scatter3(xpred(:,1),xpred(:,2),Ypred(:,1),100,'k','filled');
xlabel('Water Fraction');
ylabel('Temperature(K)');
zlabel('Eruption height(m)');
title('Emulator for eruption height');


%Output 2: Particle Flux

figure
tri = delaunay(x(:,1),x(:,2));
trisurf(tri,x(:,1),x(:,2),y(:,2));

% Clean it up
l = light('Position',[-50 -15 29]);
lighting none
shading interp
colorbar EastOutside

hold on
scatter3(xpred(:,1),xpred(:,2),Ypred(:,2),100,'k','filled');
xlabel('Water Fraction');
ylabel('Temperature(K)');
zlabel('Particle Flux');
title('Emulator for Particle Flux');









 
            
