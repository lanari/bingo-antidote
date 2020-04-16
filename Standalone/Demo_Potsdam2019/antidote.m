%% ANTIDOTE-Standalone
% Antidote Test 2.0 (Last modification 17.05.15 - eduester)
% This script loads the input file and calls BINGO over a given P-T range.
%
clc, close all, clear all
more off
%
%
%==========================================================================
%% USER INPUT
%
% Define inout file:
fid = fopen('input_xmaptools2.txt');
% Define bulk composition file:
WorkVariXMap.bcomp='THERIN';
% Define database file:
WorkVariXMap.dbs = 'JUN92.bs';
% Define range of temperature [?C]:
Tmin=500     ;
Tmax=700     ;
% Define range of pressure [bar]:
Pmin=4000    ;
Pmax=10000   ;
% -------------------------------
% T-P grid
[T,press] = meshgrid(Tmin:40:Tmax,Pmin:500:Pmax);
%
%==========================================================================
%
%
% read input from "input-file":
WorkVariXMap.NbPhases = str2num(fgetl(fid));
WorkVariXMap.NbEl = str2num(fgetl(fid));
WorkVariXMap.Els = strread(fgetl(fid),'%s')';

for i=1:WorkVariXMap.NbPhases
    WorkVariXMap.COMPoxy(i,:) = strread(fgetl(fid),'%f');
end

WorkVariXMap.Names = strread(fgetl(fid),'%s')';
WorkVariXMap.VolFrac = strread(fgetl(fid),'%f')';
WorkVariXMap.Format = strread(fgetl(fid),'%s')';

WorkVariXMap.VolFrac = WorkVariXMap.VolFrac ./ sum(WorkVariXMap.VolFrac);
fclose(fid);

% Optional: Define parameter for optimization: 
WorkVariXMap.Minim='Qtot1';

% Convert Oxides to atoms[mol]: 
WorkVariXMap=MinConverter(WorkVariXMap);
%WorkVariXMap.COMP

%%
% initialization of variables:
A1=zeros(size(T));
V1=zeros(size(T));
C1=zeros(size(T));
TOTAL=zeros(size(T));

% display output
disp(' ')
disp('...  BINGO-ANTIDOTE standalone ...')
disp(' ')
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp(['>>> New Run: ',datestr(now),'  <<<'])
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp(' ')
%%
% START Calculation:

h = waitbar(0,'Please wait...');
for i = 1:1:size(T,1) 
    for j = 1 :1:  size(T,2)
        % CALL BINGO
        [Ymin, Evaluation]=Bingo([T(i,j),press(i,j)],WorkVariXMap,0);
        % store BINGO OUTPUT
        disp([T(i,j) press(i,j)/1000])
        A1(i,j)=Evaluation.assemblage;
        V1(i,j)=Evaluation.Volume;
        C1(i,j)=Evaluation.Compositions;
        TOTAL(i,j)=Evaluation.Total;
        MIN(i,j,:)=Evaluation.MinComp;     
    end
    waitbar(i/size(T,1) ,h)
    drawnow
end
close(h)

%% 
% Minimizing Qtot

X0 = [T(find(TOTAL==max(max(TOTAL)))), press(find(TOTAL==max(max(TOTAL))))];
WorkVariXMap.Minim='Qtot1';
minphase=WorkVariXMap;
f = @(x) Bingo(x,minphase,0);
        options = optimset('fminsearch'); options=optimset(options,'TolX',0.0001,'TolFun',0.0001,'display','off','MaxFunEvals',300,'MaxIter',1000);
[XminT1,Val]=fminsearch(f, X0, options);

% optimization output:
 [C2, Evaluation]=Bingo(XminT1,WorkVariXMap);
 Evaluation

% X0 = [T(find(V1==max(max(V1)))),press(find(V1==max(max(V1))))];
% WorkVariXMap.Minim='Qvol1';
% minphase=WorkVariXMap;
% f = @(x) Bingo(x,minphase,0);
%         options = optimset('fminsearch'); options=optimset(options,'TolX',0.0001,'TolFun',0.0001,'display','off','MaxFunEvals',300,'MaxIter',1000);
% [XminV1,Val]=fminsearch(f, X0, options);
% 
% X0 = [T(find(C1==max(max(C1)))),press(find(C1==max(max(C1))))];
% WorkVariXMap.Minim='Qcmp1';
% minphase=WorkVariXMap;
% f = @(x) Bingo(x,minphase,0);
%         options = optimset('fminsearch'); options=optimset(options,'TolX',0.0001,'TolFun',0.0001,'display','off','MaxFunEvals',300,'MaxIter',1000);
% [XminC1,Val]=fminsearch(f, X0, options);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 1 - Quality factors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
colormap(flipud(hot))
caxis([0 100])
colorbar
subplot(2,2,1), h=pcolor(T,press/1000,A1);
set(h,'EdgeColor','none');
xlabel('Temperature [?C]')
ylabel('pressure [kbar]')
title('Evaluation Assemblage')
colormap(flipud(hot))
caxis([0 100]);
colorbar
hold on
%plot(T(find(A1==max(max(A1)))),press(find(A1==max(max(A1))))/1000,'bo')
hold off

subplot(2,2,2)
[c,h]=contourf(T,press/1000,V1,0:5:100);
clabel(c,h,'FontSize',12,'Color','k','LabelSpacing',200);
%h=pcolor(T,press/1000,V1);
%set(h,'EdgeColor','none')
xlabel('Temperature [?C]')
ylabel('pressure [kbar]')
title('Evaluation Volume fraction')
colormap(flipud(hot))
caxis([50 100]);
colorbar
hold on
plot(T(find(V1==max(max(V1)))),press(find(V1==max(max(V1))))/1000,'bo')
hold off

subplot(2,2,3)
%h=pcolor(T,press/1000,C1);
%set(h,'EdgeColor','none')
[c,h]=contourf(T,press/1000,C1,0:5:100);
clabel(c,h,'FontSize',12,'Color','k','LabelSpacing',200);
xlabel('Temperature [?C]')
ylabel('pressure [kbar]')
title('Evaluation Mineral composition')
colormap(flipud(hot))
caxis([50 100]);
colorbar
hold on
plot(T(find(C1==max(max(C1)))),press(find(C1==max(max(C1))))/1000,'bo')
hold off

subplot(2,2,4)
[c,h]=contourf(T,press/1000,TOTAL,0:10:100);
clabel(c,h,'FontSize',12,'Color','k','LabelSpacing',200);
%h=pcolor(T,press/1000,TOTAL);
%set(h,'EdgeColor','none')
xlabel('Temperature [?C]')
ylabel('pressure [kbar]')
title('TOTAL Likelihood')
colormap(flipud(hot))
caxis([50 100]);
colorbar
hold on
plot(T(find(TOTAL==max(max(TOTAL)))),press(find(TOTAL==max(max(TOTAL))))/1000,'bo')
plot(XminT1(1),XminT1(2)./1000,'go')
hold off

%%
% Piechart
figure(2)
subplot(1,2,1)
pie(Evaluation.VolXMap,Evaluation.PhasesNames)
title('Observed')
subplot(1,2,2)
pie(Evaluation.VolMod,Evaluation.PhasesNames)
title(['Modeled at ',num2str(round(XminT1(1))),'?C and ',num2str(round(XminT1(2))./1000),'kbar'])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 2 - Qcmp of minerals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=1:WorkVariXMap.NbPhases;
if mod(WorkVariXMap.NbPhases,3)==0
    v=floor(WorkVariXMap.NbPhases./3);
else
    v=floor(WorkVariXMap.NbPhases./3)+1;
end

for i=1:WorkVariXMap.NbPhases
figure(3)
subplot(v,3,u(i))
h=pcolor(T,press/1000,MIN(:,:,i));
set(h,'EdgeColor','none');
colormap(jet)
title(WorkVariXMap.Names(i))
caxis([0 100]);

end
colorbar
for i=1:WorkVariXMap.NbPhases
figure(4)
subplot(v,3,u(i))
[c,h]=contourf(T,press/1000,MIN(:,:,i),0:10:100);
clabel(c,h,'FontSize',12,'Color','k','LabelSpacing',200);
%colormap(flipud(hot))
colormap(jet)
title(WorkVariXMap.Names(i))
caxis([0 100]);
end

%%
% Minimizing for a given parameter
if 0
%% Simplex optimization

% define initial values for optimization:
 X0 = [610,7000];
% define value that shoould be optimized:
 WorkVariXMap.Minim='Qtot1';
% give input parameters:
 minphase=WorkVariXMap;
% run optimization: 
 f = @(x) Bingo(x,minphase,0);
        options = optimset('fminsearch'); options=optimset(options,'TolX',0.0001,'TolFun',0.0001,'display','off','MaxFunEvals',300,'MaxIter',1000);
 [Xmin,Val]=fminsearch(f, X0, options);

% optimization output:
 [C2, Evaluation]=Bingo(Xmin,WorkVariXMap);
 Evaluation
end
