function [Best,Bestf,conv,conve,timf]=algen(g,xmin,xmax,Npop,Ngen,pmut,Nsic,txt)
% Algoritmo genético (GA)
% Entrada:
% g = función a optimizar
% xmin = vector de cotas inferiores de los parámetros
% xmax = vector de cotas superiores de los parámetros
% Ngen = número de generaciones
% Npop = tamaño de la población
% pmut = probabilidad de mutación de los componentes de cada hijo
% Nsic = criterio de parada: si el mejor elemento no cambia tras Nsic generaciones, se detiene el AG
% txt = 1 o 0 para mostrar/ocultar el estado del algoritmo
%
% Salida:
% Bestf = mínimo de f
% Best  = argmínimo
% conv  = evolución del mejor valor de cada generación
% conve = evolución del mejor elemento de cada generación
% timf  = tiempo de ejecución (s)
tic
if txt==1
disp(['Genetic Algorithm']);
end
Ndim=max(size(xmin));
contcp=0; % Contador para el criterio de parada
Besteo=NaN;
for it=1:Ngen
if it==1
% Población inicial
vec=rand(Npop,Ndim);
Pop=xmin.*vec+xmax.*(1-vec);
end
% Evaluación de la población
for i=1:Npop
evalfP(i)=g(Pop(i,:)); % Evaluación de la población inicial
end
% Elitismo
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
break % Detener las iteraciones del AG
end
% Selección de individuos
% Creación de la probabilidad de selección
for i=1:Npop
prob(i)=(1/(evalfP(i)-min(evalfP)+1))/(sum(1./(evalfP-min(evalfP)+1)));
end
% Selección de 2*Npop padres
for i=1:2*Npop
n=rand;
cont=0; j=0;
while(cont<n)
j=j+1;
cont=cont+prob(j);
end
parent(i,:)=Pop(j,:);
end
% Creación de descendientes
for i=1:Npop
rr=rand;
child(i,:)=rr*parent(2*i-1,:)+(1-rr)*parent(2*i,:);
end
% Mutación de descendientes
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
% Evaluación de la última población
for i=1:Npop
evalfP(i)=g(Pop(i,:)); % Evaluación de la población
end
[a,b]=min(evalfP);
Best=Pop(b,:);
Bestf=a;
conv(it+1)=a;
conve(it+1,1:Ndim)= Best;
timf=toc; % Tiempo de ejecución
end
