function [Best,Bestf,conv,conve,timf]=algen(g,xmin,xmax,Npop,Ngen,pmut,Nsic,txt)
% Genetic Algorithm (GA)
% Input:
% g=function to be optimized
% xmin = Lower Bound Vector of Parameter Values
% xmax= Upper Bound Vector of Parameter Values
% Ngen= Number of generations
% Npop= Population Size
% pmut= probability of mutation of the components of each child
% Nsic= Stop Criterion: If the best element don't change after Nsic generations,stop the AG
% txt=1 or 0 to report the state of the algorithm
%
% Output
% Bestf= minimum of f
% Beste= argminimum
% conv = evolution of the best element of each generation
% timf = runtime(s)
tic
if txt==1
disp(['Genetic Algorithm']);
end
Ndim=max(size(xmin));
contcp=0; % Counter for Stopping Criterion
Besteo=NaN;
for it=1:Ngen
if it==1
% Initial population
vec=rand(Npop,Ndim);
Pop=xmin.*vec+xmax.*(1-vec);
end
% Population Assessment
for i=1:Npop
evalfP(i)=g(Pop(i,:)); %Initial Population Assessment
end
% Elitism
[a,b]=min(evalfP);
Best=Pop(b,:);
Bestf=a;
conv(it)=a;
conve(it,1:Ndim)= Best;
if Best==Besteo
contcp=contcp+1;
else
contcp=0;
end
if txt==1
disp(['Gen : ' num2str(it) '/' num2str(Ngen) ' | Best element: ' ...
num2str(Bestf) ' | Gen without changes: ' num2str(contcp)])
end
if contcp>=Nsic
disp(['Stop: Cannot improve best element']);
break % Stop the AG iterations
end
% Selection of individuals
% Creating the Selection Probability
for i=1:Npop
prob(i)=(1/(evalfP(i)-min(evalfP)+1))/(sum(1./(evalfP-min(evalfP)+1)));
end
%Selection of 2*Npop parents
for i=1:2*Npop
n=rand;
cont=0; j=0;
while(cont<n)
j=j+1;
cont=cont+prob(j);
end
parent(i,:)=Pop(j,:);
end
%Creation of descendants
for i=1:Npop
rr=rand;
child(i,:)=rr*parent(2*i-1,:)+(1-rr)*parent(2*i,:);
end
%Mutation of descendants
for i=1:Npop
for j=1:Ndim
if rand<pmut
rr=rand;
child(i,j)=xmin(j).*rr+xmax(j).*(1-rr);
end
end
end
child(end,:)=Best;
Pop=child;
Besteo=Best;
if it==Ngen
disp(['Last Generation Stop']);
end
end
% Evaluation of the last population
for i=1:Npop
evalfP(i)=g(Pop(i,:)); %Initial Population evaluation
end
[a,b]=min(evalfP);
Best=Pop(b,:);
Bestf=a;
conv(it+1)=a;
conve(it+1,1:Ndim)= Best;
timf=toc; % Runtime
end