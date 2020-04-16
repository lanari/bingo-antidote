%% BINGO-Standalone
%
% Bingo 2.0 (Last modification 17.05.15 - eduester)
%
%
function [OUT,Evaluation] = Bingo(X,WorkVariXMap,DoWePrint)

setenv('DYLD_LIBRARY_PATH', '/usr/local/bin')
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')

if ispc
    InvMet.path = 'BA_ProgramsWIN\theriakd.exe';
    %C:\TheriakDominoWIN_2017\Programs\theriakd
else
    InvMet.path = 'BA_ProgramsMAC/theriakd';    
end    
InvMet.In = 'in';
InvMet.bcomp =WorkVariXMap.bcomp;
InvMet.dbs= WorkVariXMap.dbs;


 if length(X) == 2
 
    D_Temp = num2str(abs(X(1)));
    D_Press = num2str(abs(X(2)));
 else
     fprintf('%s\n','--------------ERROR---------------');
     fprintf('%s\n','number of input arguments is not 2');
     fprintf('%s\n','Bingo uses default:               ');
     fprintf('%s\n','temperature = 500 ?C              ');
     fprintf('%s\n','pressure= 4 kbar)                 ');
     fprintf('%s\n','----------------------------------');
     D_Temp = '700';
     D_Press = '7000';
 end

% -------------------------------------------------------------------------
% Call Theriakd with the instruction: BINGO
Command = 'BINGO'; 
INFILE=char( InvMet.dbs, D_Temp,  D_Press , Command);
dlmwrite(InvMet.In , INFILE,'delimiter','');
address = [InvMet.path ' ' InvMet.In ' ' InvMet.bcomp];
[wum,yum]=system(address);
[WorkVariMod] = ReadResTheriakd(yum);
[WorkVariXMap] = CutMinNames(WorkVariXMap);

% -------------------------------------------------------------------------
%% READ XTT_Min_files

[DefMin,DefGf] = BingoInitialization(InvMet);
%%
% -------------------------------------------------------------------------
% Evaluation CHECK:
if exist('DoWePrint')
    DoWePrint = DoWePrint;
else    
    DoWePrint = 1;
end

if DoWePrint
    disp('===================')
    disp('...  BINGO RUN  ...')
    disp('===================')
    disp(' ')
end

% Step 1 - Assemblage
[Evaluation.assemblage,Link] = ComputeQualityAssemblage(WorkVariMod,WorkVariXMap,DoWePrint);

% Step 2 - Volume
[Evaluation] = ComputeQualityVolume(WorkVariMod,WorkVariXMap,Link,Evaluation,DoWePrint);

% Step 3 - Compositions
[Evaluation] = ComputeQualityCompositions(WorkVariMod,WorkVariXMap,Link,Evaluation,DoWePrint,DefMin);



% store Quality factors:
A1=Evaluation.assemblage;
V1=Evaluation.Volume;
C1=Evaluation.Compositions;
TOTAL2= (Evaluation.assemblage+Evaluation.assemblage/100*Evaluation.Volume + Evaluation.assemblage/100*Evaluation.Compositions)/3;
Evaluation.Total=TOTAL2;


% check which parameter should be outputed for optimization:
V1=100-V1;
C1=100-C1;
TOTAL2=100-TOTAL2;
disp(['...  BINGO  =>  ', num2str(Evaluation.Total)])


if WorkVariXMap.Minim=='Qcmp1'
    Evaluation.MinimVar='Qcmp1';
    OUT=C1;
elseif WorkVariXMap.Minim=='Qvol1'
    Evaluation.MinimVar='Qvol1';
    OUT=V1;
elseif WorkVariXMap.Minim=='Qcmp2'
    Evaluation.MinimVar='Qcmp2';
    OUT=TOTAL2.*C1;
elseif WorkVariXMap.Minim=='Qvol2'
    Evaluation.MinimVar='Qvol2';
    OUT=TOTAL2.*V1;
else
    Evaluation.MinimVar='Qtot';
    OUT=TOTAL2;
end 

if isnan(OUT)
  OUT=0;
end

% end of MAIN Function
return

% =========================================================================
%% SUB-FUNCTIONS:
% =========================================================================
function [Result,Link] = ComputeQualityAssemblage(WorkVariMod,WorkVariXMap,DoWePrint)
% ComputeQualityAssemblage is a function to estimate the Evaluation of the  
% model for the stable assemblage.
%
% To be described later once the strategy is set.
% 
% Last change by Pierre Lanari (30.11.2015)
%
PhasesTher = WorkVariMod.Names;
PhasesXMap = WorkVariXMap.Names;

% Generate the variable LINK to be used later
ComptPhase=0;
Link.PhasesNames = {''};

%eduester
for i=1:length(PhasesTher)
    if ~ismember(Link.PhasesNames,PhasesTher{i})
        ComptPhase = ComptPhase+1;
        Link.PhasesNames{ComptPhase} = PhasesTher{i};
        PhasesTher{i}=PhasesTher{i};
    end
end
for i=1:length(PhasesXMap)
    PhasesXMap{i}=PhasesXMap{i};
    if ~ismember(Link.PhasesNames,PhasesXMap{i})
        ComptPhase = ComptPhase+1;
        Link.PhasesNames{ComptPhase} = PhasesXMap{i};
    end
end
[Link.TherIsIn,Link.TherIndices] = ismember(Link.PhasesNames,PhasesTher);
[Link.XMapIsIn,Link.XMapIndices] = ismember(Link.PhasesNames,PhasesXMap);


NbPhasesInvolved = length(Link.TherIsIn);
NbMissingPhases = 2*NbPhasesInvolved-(length(find(Link.TherIsIn))+length(find(Link.XMapIsIn)));

NbPhasesMax = max([length(PhasesTher),length(PhasesXMap)]);
NbMatches = length(find(ismember(PhasesTher,PhasesXMap)));

Result =  NbMatches/NbPhasesInvolved*100;


if DoWePrint
    Code1 = '%s\t\t';
    Code2 = '%s\t\t';
    for i=1:length(Link.PhasesNames)
        Len = length(Link.PhasesNames{i});
        if Len < 8
            for j=Len:8
                Link.PhasesNames{i}=[char(Link.PhasesNames{i}),' '];
            end
        end    
        Code1=[Code1,'%s\t'];
        Code2=[Code2,'%f\t'];
    end
    Code1(end) = 'n';
    Code2(end) = 'n';
    
    fprintf('%s\n','##### Evaluation critrion (1) ASSEMBLAGE ##### ');
    fprintf(Code1,'Phases:',Link.PhasesNames{:});
    fprintf(Code2,'THER:',Link.TherIsIn);
    fprintf(Code2,'XMAP:',Link.XMapIsIn);
    fprintf('%s\n','  ');
    fprintf('%s\n','-------------');
    fprintf('%s\t%.0f\n','n =',length(PhasesTher));
    fprintf('%s\t%.0f\n','m =',length(PhasesXMap));
    fprintf('%s\t%.0f\n','l =',NbMatches);
    fprintf('%s\t%.2f\n','Q =',Result);  
    fprintf('%s\n','-------------');
end

return


function [Evaluation] = ComputeQualityVolume(WorkVariMod,WorkVariXMap,Link,Evaluation,DoWePrint)
% ComputeVolumeAssemblage is a function to estimate the Evaluation of the  
% model for the volume of stable phases.
%
% To be described later once the strategy is set.
% 
% Last change by Pierre Lanari (30.11.2015)
%

VolTherMatch = zeros(size(Link.PhasesNames));
VolXMapMatch = zeros(size(Link.PhasesNames));

VolTherMatch(find(Link.TherIsIn)) = WorkVariMod.VolFrac(Link.TherIndices(find(Link.TherIndices)));
VolXMapMatch(find(Link.XMapIsIn)) = WorkVariXMap.VolFrac(Link.XMapIndices(find(Link.XMapIndices)));

VolTherMatchNorm = VolTherMatch ./ sum(VolTherMatch);
VolXMapMatchNorm = VolXMapMatch ./ sum(VolXMapMatch);


Result_Diff = 100.*sqrt(sum((VolTherMatch+VolXMapMatch)/2.*(1-abs(VolTherMatch-VolXMapMatch)./(VolTherMatch+VolXMapMatch)).^2));


Evaluation.Volume = Result_Diff;

if DoWePrint
    Code1 = '%s\t\t';
    Code2 = '%s\t\t';
    for i=1:length(Link.PhasesNames)
        Len = length(Link.PhasesNames{i});
        if Len < 8
            for j=Len:8
                Link.PhasesNames{i}=[char(Link.PhasesNames{i}),' '];
            end
        end    
        Code1=[Code1,'%s\t'];
        Code2=[Code2,'%f\t'];
    end
    Code1(end) = 'n';
    Code2(end) = 'n'; 
    
    disp(' ')
    fprintf('%s\n','##### Evaluation critrion (2) VOLUME FRACTIONS ##### ');
    fprintf(Code1,'Phases:',Link.PhasesNames{:})
    fprintf(Code2,'THER:',VolTherMatch);
    fprintf(Code2,'XMAP:',VolXMapMatch);
    fprintf(Code2,'abs(D):',abs(VolTherMatch-VolXMapMatch));
    fprintf('%s\n','  ');
    fprintf('%s\n','-------------');
    fprintf('%s\t%.2f\n','Q =',Result_Diff);  
    fprintf('%s\n','-------------');
end

%store Vol. fractions for pie chart:
Evaluation.PhasesNames=Link.PhasesNames;
Evaluation.VolXMap=VolXMapMatch;
Evaluation.VolMod=VolTherMatch;
return


function [Evaluation] = ComputeQualityCompositions(WorkVariMod,WorkVariXMap,Link,Evaluation,DoWePrint,DefMin)
% ComputeQualityCompositions is a function to estimate the Evaluation of the  
% model for the composition of stable phases.
%
% To be described later once the strategy is set.
% 
% Last change by Erik Duesterhoeft (08.03.2019)
%
% CutOffThreshold (apfu value)
CutOffThreshold = 0.005;

NbElemsTher = WorkVariMod.NbEl;   
ElemsListTher = WorkVariMod.Els;

NbElemsXMap = WorkVariXMap.NbEl;
ElemsListXMap = WorkVariXMap.Els;

NbElems = NbElemsTher;
ElemsList = ElemsListTher;

[Ok,WhereIsMem] = ismember(ElemsListTher,ElemsListXMap);

WorkVariXMap.COMPok = zeros(WorkVariXMap.NbPhases,NbElemsTher);
WorkVariXMap.COMPok(:,find(Ok)) = WorkVariXMap.COMP(:,WhereIsMem(find(Ok)));


% Correction of O for H in the structural formula:
WhereH = find(ismember(ElemsList,'H'));
if WhereH
    TheHValues = WorkVariMod.COMP(:,WhereH);
    WorkVariMod.COMP(1:WorkVariMod.NbPhases,WhereH) = zeros(WorkVariMod.NbPhases,1);
else
    TheHValues = 0;
end

WhereO = find(ismember(ElemsList,'O'));
if WhereO
    WorkVariMod.COMP(1:WorkVariMod.NbPhases,WhereO) = WorkVariMod.COMP(1:WorkVariMod.NbPhases,WhereO)-0.5*TheHValues;
end


CompTherMatch = zeros(size(Link.PhasesNames,2),NbElems);
CompXMapMatch = zeros(size(Link.PhasesNames,2),NbElems);
YesMat = zeros(size(Link.PhasesNames,2),NbElems);


if size(WorkVariMod.COMP,1) < 2
    keyboard
end 

CompTherMatch(find(Link.TherIsIn),:) = WorkVariMod.COMP(Link.TherIndices(find(Link.TherIndices)),:);
CompXMapMatch(find(Link.XMapIsIn),:) = WorkVariXMap.COMPok(Link.XMapIndices(find(Link.XMapIndices)),:);




% 4.0 Perspectives version (Version 8) no weights (eduester)
if 1

    Qual = zeros(size(CompTherMatch,1),1);
    Missfit = zeros(size(CompTherMatch,1),1);
    weight = zeros(size(CompTherMatch,1),1);
    
    for i=1:size(CompTherMatch,1)
        NbAtXMap = sum(CompXMapMatch(i,2:end-1)); % I excluded oxygen (2:end) and E (end)
        NbAtTher = sum(CompTherMatch(i,2:end-1)); % and E (end)
         
        if NbAtXMap > 0 && NbAtTher > 0 
            

            % Extract the element list (to be compared) from BinPhaseDef
            PhasesTher = WorkVariMod.Names;
            PhasesDB=DefMin.DBNames';
            for ii=1:length(PhasesDB)
                    if length(PhasesDB{ii})<4
                    PhasesDB{ii}=[PhasesDB{ii}(1:3) ' '];
                    else
                    PhasesDB{ii}=PhasesDB{ii};
                    end
            end

            WhereInBinDef = find(ismember(strtrim(PhasesDB),PhasesTher{i}));
            ElsSel=DefMin.Els;
            [Yes,WhereElemInList] =ismember(WorkVariMod.Els,cellstr(ElsSel{WhereInBinDef}));
            YesMat(i,:)=Yes;

            
            
            DIFFabs = abs(CompTherMatch(i,2:end-1)-CompXMapMatch(i,2:end-1));
            
            Fac1=1.*ones(size(DIFFabs)); % plateau at 1 sigma
            Fac2=5.*ones(size(DIFFabs)); %
            %guessed uncertainties
            sigma=2.0;
            stds=normcdf2(1,0,1/sigma)-(1-normcdf2(1,0,1/sigma));
            Unc1s = CompXMapMatch(i,2:end-1).*(1-stds);
            Unc1s(find(Unc1s<0.01))=0.01; 
            
            

            DIFF2=DIFFabs-(Unc1s./Fac1);
            WhereNEGZ=find(DIFF2<=0); %where negative
            DIFF2(WhereNEGZ)=zeros(size(WhereNEGZ));
 
            if sum(DIFF2>Fac2.*Unc1s)>0
            WherePOSZ=find(DIFF2>Fac2.*Unc1s); %where positive
            DIFF2(WherePOSZ)=(Fac2(WherePOSZ).*Unc1s(WherePOSZ)).*ones(size(WherePOSZ)); %hier war fehler Unc1s(1)
            end
            
            
            
            if (sum(DIFFabs)==0)
            wj=1./(1-(DIFFabs./1));
            else
            wj=1./(1-(DIFFabs./sum(DIFFabs)));
            end
            % here we deactivate weights:
            wj=Yes(2:end-1);
            %wj=wj.*Yes(2:end-1);

            

            Qual(i) =100.*sum(((1-DIFF2./(Fac2.*Unc1s)).^(CompTherMatch(i,2:end-1)+1)).*(wj)./sum(wj));
            %Qual(i) =100.*mean(((1-DIFF2./(Fac2.*Unc1s)).^(CompTherMatch(i,2:end-1)+1)));  
            
            

        end
        
        weight(i) = max([CompXMapMatch(i,1),CompTherMatch(i,1)]);
        Missfit(i) = sum(abs(CompTherMatch(i,:)-CompXMapMatch(i,:)))*weight(i);

    end
    
end



     
if sum(Qual)  > 0
    
    % Qcmp with normalized VolFrac from Theriak (Qi):    
    VolTherMatch(find(Link.TherIsIn)) = WorkVariMod.VolFrac(Link.TherIndices(find(Link.TherIndices)));
    VolTherMatchNorm=VolTherMatch(find(Qual))./sum(VolTherMatch(find(Qual)));
    Evaluation.Compositions=sum(Qual(find(Qual)).*VolTherMatchNorm'); % normalized, only the ones that matches.
    
        
else
    Evaluation.Compositions = 0;
end


% store Qcmp of mineral phases:
for i=1:length(WorkVariXMap.Names)
    if ismember(WorkVariXMap.Names(i),WorkVariMod.Names)  
        Evaluation.MinComp(i)=Qual(find(ismember(strtrim(Link.PhasesNames),WorkVariXMap.Names(i))));
    else
        Evaluation.MinComp(i)=0;
    end
end

% print output to command window/screen:
if DoWePrint
    
    Code1 = '%s\t\t';
    Code2 = '%s\t\t';
    for i=1:NbElems    
        Code1=[Code1,'%s\t\t'];
        Code2=[Code2,'%f\t'];
    end
    Code1(end) = 'n';
    Code2(end) = 'n'; 
    
    
    disp(' ')
    fprintf('%s\n','##### Evaluation critrion (3) PHASE COMPOSITIONS ##### ');
    disp('-')
    for i=1:length(Link.PhasesNames)
        
        disp(char(Link.PhasesNames{i}))

        fprintf(Code1,'Elements:',ElemsList{:})
        fprintf(Code2,'THER:',CompTherMatch(i,:));
        fprintf(Code2,'XMAP:',CompXMapMatch(i,:));
        fprintf(Code2,'abs(D):',abs(CompTherMatch(i,:)-CompXMapMatch(i,:)));
        www=[0,(1./(1-(abs(CompTherMatch(i,2:end-1)-CompXMapMatch(i,2:end-1))./sum(abs(CompTherMatch(i,2:end-1)-CompXMapMatch(i,2:end-1)))))),0];
        if (sum(abs(CompTherMatch(i,2:end-1)-CompXMapMatch(i,2:end-1)))==0)
            www=[0,(1./(1-(abs(CompTherMatch(i,2:end-1)-CompXMapMatch(i,2:end-1))))),0];
        end

        fprintf(Code2,'w(abs):',www);
        fprintf(Code2,'ElsSel:',YesMat(i,:));        
        www=(www.*YesMat(i,:))./sum(www.*YesMat(i,:));
        www(isnan(www))=0;
        fprintf(Code2,'w(rel):',www);
        fprintf('%s\t\t%.2f\n','Qi',Qual(i));
        disp('-')
    end
    
    fprintf('%s\n','-------------');
    fprintf('%s\t%.2f\n','Q =',Evaluation.Compositions);  
    fprintf('%s\n','-------------');
    
end


return




function [WorkVariMod] = ReadResTheriakd(OutputTheriakd)
% ReadResTheriakd is a function that read the result of theriakd and
% generate a variable WorkVariMod.
%
% -> Theriakd result: 
%   -----------------------------------
%   | Nb of phases                    | 
%   | Nb of system components         |
%   | List of system components       |
%   | Matrix of COMPOSITION of phases |
%   | List of phase names             |
%   | Volumes persentages of solids   |
%   | Density of solids               |
%   | Nb of solid solutions           |
%   -----------------------------------
%
% -> ReadResTheriakd output variable: 
%
% WorkVariMod.NbPhases      number of stable solid phases 
% WorkVariMod.NbEl          number of components
% WorkVariMod.Els           list of component names (format: database)
% WorkVariMod.COMP          composition of stable solid phases ...
%                           row order: .Names  / Column order: .Els
% WorkVariMod.VolFrac       volume fractions of solids with normalized sum
% WorkVariMod.Dens          specific density of solids (to be use for LEB)
% WorkVariMod.NbSS          number of solids: solid solution
% WorkVariMod.NbPP          number of solids: pure phase
%
% Last change by Pierre Lanari (22.11.2015)
%


TestInput = strread(OutputTheriakd,'%s');

WorkVariMod(1).NbPhases = str2num(char(TestInput(1)));
if isempty(WorkVariMod(1).NbPhases)
    disp(['ERROR: ',OutputTheriakd])
    keyboard
end    

WorkVariMod(1).NbEl = str2num(char(TestInput(2)));

for i=1:WorkVariMod(1).NbEl
    WorkVariMod(1).Els{i} = char(TestInput(2+i));
end

WereForCompo = 3 + WorkVariMod(1).NbEl;
WereForNames = 3 + WorkVariMod(1).NbEl + WorkVariMod(1).NbPhases*WorkVariMod(1).NbEl;
WereForVol = 3 + WorkVariMod(1).NbEl + WorkVariMod(1).NbPhases*WorkVariMod(1).NbEl + WorkVariMod(1).NbPhases;
WereForDens = 3 + WorkVariMod(1).NbEl + WorkVariMod(1).NbPhases*WorkVariMod(1).NbEl + WorkVariMod(1).NbPhases + WorkVariMod(1).NbPhases;

for i=1:WorkVariMod(1).NbPhases
    WorkVariMod(1).COMP(i,:) = str2num(char(TestInput(WereForCompo+(i*WorkVariMod(1).NbEl-WorkVariMod(1).NbEl):WereForCompo+(i*WorkVariMod(1).NbEl)-1)))';
    WorkVariMod(1).Names{i} = char(TestInput(WereForNames+i-1));
    WorkVariMod(1).VolFrac(i) = str2num(char(TestInput(WereForVol+i-1)));
    WorkVariMod(1).Dens(i) = str2num(char(TestInput(WereForDens+i-1)));
end

WorkVariMod(1).NbSS = str2num(char(TestInput(end)));
WorkVariMod(1).NbPP = WorkVariMod(1).NbPhases-WorkVariMod(1).NbSS;

WorkVariMod(1).VolFrac = WorkVariMod(1).VolFrac ./ sum(WorkVariMod(1).VolFrac);


[WorkVariMod] = CutMinNames(WorkVariMod);

return


function [WorkVariXMap] = ReadXMapTemp(FileName)
% ReadXMapTemp reads FileName and generates a fake variable WorkVariXMap. 
% 
% Temporary function to be remove later
%

fid = fopen(FileName);

WorkVariXMap.NbPhases = str2num(fgetl(fid));
WorkVariXMap.NbEl = str2num(fgetl(fid));

WorkVariXMap.Els = strread(fgetl(fid),'%s')';

for i=1:WorkVariXMap.NbPhases
    WorkVariXMap.COMP(i,:) = strread(fgetl(fid),'%f');
end

WorkVariXMap.Names = strread(fgetl(fid),'%s')';
WorkVariXMap.VolFrac = strread(fgetl(fid),'%f')';

WorkVariXMap.VolFrac = WorkVariXMap.VolFrac ./ sum(WorkVariXMap.VolFrac)


fclose(fid);

return


function [DefMin,DefGf] = ReadMinTrans2(MinTransFile);
DefMin = [];
DefGf = [];

fid = fopen(MinTransFile); Compt = 0;
while 1
    lalign = fgetl(fid);
    if isequal(lalign,-1)
        break
    end
    
    if length(lalign) > 1
        if isequal(lalign(1:2),'>1')
            Compt = 0;
            while 1
                lalign = fgetl(fid);
                if isequal(lalign,'')
                    break
                end
                Compt = Compt+1;
                Temp = strread(lalign,'%s');
                %eduester octave fit
                DefMin.Names{Compt} = Temp{1};
                DefMin.NbOx{Compt} = str2num(Temp{2});
                DefMin.DBNames{Compt} = Temp{3};
                DefMin.Els{Compt} = strread(Temp{4}(2:end-1),'%s','delimiter',',')';
                %DefMin{Compt,5} = str2double(strread(Temp{5}(2:end-1),'%s','delimiter',','))';
            end
        end
        
    end
    
    if length(lalign) > 1
        if isequal(lalign(1:2),'>2')
            Compt = 0;
            while 1
                lalign = fgetl(fid);
                if isequal(lalign,'')
                    break
                end
                Compt = Compt+1;
                Temp = strread(lalign,'%s');
                DefGf{Compt,1} = Temp{1};
                DefGf{Compt,2} = str2num(Temp{2});
                DefGf{Compt,3} = Temp{3};
                DefGf{Compt,4} = strread(Temp{4}(2:end-1),'%s','delimiter',',')';
                DefGf{Compt,5} = str2double(strread(Temp{5}(2:end-1),'%s','delimiter',','))';
            end
        end
    end
end
fclose(fid);
return

function [varargout] = normcdf2(x,varargin)

if nargin<1
    error(message('stats:normcdf:TooFewInputsX'));
end

if nargin>1 && strcmpi(varargin{end},'upper')
    % Compute upper tail and remove 'upper' flag
    uflag=true;
    varargin(end) = [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;  
end

[varargout{1:max(1,nargout)}] = localnormcdf(uflag,x,varargin{:});


function [p,plo,pup] = localnormcdf(uflag,x,mu,sigma,pcov,alpha)

if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>1
   if nargin<5
      error(message('stats:normcdf:TooFewInputsCovariance'));
   end
   if ~isequal(size(pcov),[2 2])
      error(message('stats:normcdf:BadCovarianceSize'));
   end
   if nargin<6
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:normcdf:BadAlpha'));
   end
end

try
    z = (x-mu) ./ sigma;
    if uflag==true
        z = -z;
    end
catch
    error(message('stats:normcdf:InputSizeMismatch'));
end

% Prepare output
p = NaN(size(z),class(z));
if nargout>=2
    plo = NaN(size(z),class(z));
    pup = NaN(size(z),class(z));
end

% Set edge case sigma=0
if uflag==true
    p(sigma==0 & x<mu) = 1;
    p(sigma==0 & x>=mu) = 0;
    if nargout>=2
        plo(sigma==0 & x<mu) = 1;
        plo(sigma==0 & x>=mu) = 0;
        pup(sigma==0 & x<mu) = 1;
        pup(sigma==0 & x>=mu) = 0;
    end
else
    p(sigma==0 & x<mu) = 0;
    p(sigma==0 & x>=mu) = 1;
    if nargout>=2
        plo(sigma==0 & x<mu) = 0;
        plo(sigma==0 & x>=mu) = 1;
        pup(sigma==0 & x<mu) = 0;
        pup(sigma==0 & x>=mu) = 1;
    end
end

% Normal cases
if isscalar(sigma)
    if sigma>0
        todo = true(size(z));
    else
        return;
    end
else
    todo = sigma>0;
end
z = z(todo);

% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
p(todo) = 0.5 * erfc(-z ./ sqrt(2));

% Compute confidence bounds if requested.
if nargout>=2
   zvar = (pcov(1,1) + 2*pcov(1,2)*z + pcov(2,2)*z.^2) ./ (sigma.^2);
   if any(zvar<0)
      error(message('stats:normcdf:BadCovarianceSymPos'));
   end
   normz = -norminv2(alpha/2);
   halfwidth = normz * sqrt(zvar);
   zlo = z - halfwidth;
   zup = z + halfwidth;

   plo(todo) = 0.5 * erfc(-zlo./sqrt(2));
   pup(todo) = 0.5 * erfc(-zup./sqrt(2));
end

function [x,xlo,xup] = norminv2(p,mu,sigma,pcov,alpha)



if nargin<1
    error(message('stats:norminv:TooFewInputsP'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>2
   if nargin<4
      error(message('stats:norminv:TooFewInputsCovariance'));
   end
   if ~isequal(size(pcov),[2 2])
      error(message('stats:norminv:BadCovarianceSize'));
   end
   if nargin<5
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:norminv:BadAlpha'));
   end
end

% Return NaN for out of range parameters or probabilities.
sigma(sigma <= 0) = NaN;
p(p < 0 | 1 < p) = NaN;

x0 = -sqrt(2).*erfcinv(2*p);
try
    x = sigma.*x0 + mu;
catch
    error(message('stats:norminv:InputSizeMismatch'));
end

% Compute confidence bounds if requested.
if nargout>=2
   xvar = pcov(1,1) + 2*pcov(1,2)*x0 + pcov(2,2)*x0.^2;
   if any(xvar<0)
      error(message('stats:norminv:BadCovarianceSymPos'));
   end
   normz = -norminv2(alpha/2);
   halfwidth = normz * sqrt(xvar);
   xlo = x - halfwidth;
   xup = x + halfwidth;
end


function [WorkVariMod] = CutMinNames(WorkVariMod)
% function to remove EM abbreviations
% 
% 
% ------------------------------------------------------------------------
% Check Theriak Names and remove EM abbreviations
% ------------------------------------------------------------------------
for i=1:length(WorkVariMod(1).Names)
    TheNameFromTheriak = WorkVariMod(1).Names{i};
    %if ~ismember(TheNameFromTheriak,ListRefMiner)
        WereD = ismember(TheNameFromTheriak,'_');
        %keyboard
        switch sum(WereD)
            case 1
                % we delete it...
                Where = find(WereD);
                WorkVariMod(1).Names{i} = TheNameFromTheriak(1:Where-1);
            case 2
                % we have to delete the first one (Compatibility with the
                % MELT model of DOUG
                Where = find(WereD);
                NameTemp = TheNameFromTheriak(1:Where(1)-1);
                if isequal(NameTemp,'LIQtc6')
                    WorkVariMod(1).Names{i} = NameTemp;
                else
                    % we delete the second one ...
                    WorkVariMod(1).Names{i} = TheNameFromTheriak(1:Where(2)-1);
                end
            case 3
                disp('Oups, too many underscores in this name contact Pierre Lanari')
                keyboard
        end
    %end
end


% ------------------------------------------------------------------------
% Check for phase demixion (not identified in ListRefMiner). 
% ------------------------------------------------------------------------
% Note this is typically caused by flat G function in complex solid 
% solution models from Roger Powell (amphiboles).  
%
% in this case we select the phase with the higher volume fraction for
% comparison with the observation and rename the other phases.
%
%                                                  Pierre Lanari (24.10.16)

for i=1:length(WorkVariMod(1).Names)
    TheName = WorkVariMod(1).Names{i};
    Ind = find(ismember(WorkVariMod(1).Names,TheName));
    
    if length(Ind)>1
        Vols = WorkVariMod.VolFrac(Ind);
        [Val,IndSort] = sort(Vols,2,'descend');
        
        for j=2:length(IndSort)
            WorkVariMod(1).Names{Ind(IndSort(j))} = [WorkVariMod(1).Names{Ind(IndSort(j))},num2str(j)];
        end    
    end
end

return

function [DefMin,DefGf] = BingoInitialization(InvMet)

BingoDefault.SelectedProgram = 0;
BingoDefault.SelectedDatabase = 0;

BingoDefault.Theriak.Path = '';                          % Not defined here
BingoDefault.Theriak.Database(1).Label = '';
BingoDefault.Theriak.Database(1).File1 = '';
BingoDefault.Theriak.Database(1).File2 = '';
BingoDefault.Theriak.Database(1).MinTrans = '';
BingoDefault.Theriak.InputBulk = [];



if ispc
    fid = fopen('BA_Dev\XTT_ConnectionConfig.txt');
else
    fid = fopen('BA_Dev/XTT_ConnectionConfig.txt','r');
end    


while 1
    LaLign = fgetl(fid);
    if isequal(LaLign,-1)
        break
    end
    if length(LaLign)
        if isequal(LaLign(1),'>')
            TheRef = str2double(LaLign(2:4));

            switch TheRef
                
                case 1      % Default Thermodynamic dataset
                    TheL = fgetl(fid);
                    Compt=0;
                    while length(TheL)
                        Compt = Compt+1;
                        %keyboard
                        TheLStr = strread(TheL,'%s');
                        BingoDefault.Theriak.Database(Compt).Label = TheLStr{1};
                        BingoDefault.Theriak.Database(Compt).File1 = TheLStr{2};
                        BingoDefault.Theriak.Database(Compt).File2 = TheLStr{3};
                        BingoDefault.Theriak.Database(Compt).MinTrans = TheLStr{4};
                        
                        TheL = fgetl(fid);
                    end  
                    
                case 2      % Default Input Compositions
                   BingoDefault.Theriak.InputBulk = fgetl(fid);
                   
            end
        end
    end
end
fclose(fid);


% READ XTT_Min_ files
InvMet.MinTrans=' ';
for i=1:length(BingoDefault.Theriak.Database)
    if ismember(BingoDefault.Theriak.Database(i).File1,InvMet.dbs)
       if ispc 
        InvMet.MinTrans=['BA_Dev\',BingoDefault.Theriak.Database(i).MinTrans];
       else
        InvMet.MinTrans=['BA_Dev/',BingoDefault.Theriak.Database(i).MinTrans];           
       end    
    end
end
if ismember(' ',InvMet.MinTrans)
    disp('ERROR: Database is not known') 
    keyboard
end
[DefMin,DefGf] = ReadMinTrans2(InvMet.MinTrans); 

return