%IMPORT DATA

% Observed MHC alleles n=27
%A =  [{'A_0201'},{'A_0205'},{'A_0101'},{'A_1101'},{'A_2402'},{'A_0301'},{'A_3101'},...
%          {'A_6801'},{'B_0401'},{'B_2702'},{'B_2705'},{'B_3701'},{'B_3801'},{'B_3901'},... 
% 		 {'B_4403'},{'B_5101'},{'B_5102'},{'B_5103'},{'B_4001'},{'B_4006'},{'B_1501'},...
% 		 {'B_0801'},{'B_5801'},{'B_5201'},{'Cw_0401'},{'Cw_0602'},{'Cw_0702'}];

[num txt raw]=xlsread('ep_p_tI.xls');
A=raw(2:end,1);

% Target antigens n=40
Q1= [{'core1a'},{'e11a'},{'e21a'},{'ns21a'},{'ns31a'},{'ns4a1a'},{'ns4b1a'},{'ns5a1a'},{'ns5b1a'},{'p71a'},...
  {'core1b'},{'e11b'},{'e21b'},{'ns21b'},{'ns31b'},{'ns4a1b'},{'ns4b1b'},{'ns5a1b'},{'ns5b1b'},{'p71b'},...
	{'core2a'},{'e12a'},{'e22a'},{'ns22a'},{'ns32a'},{'ns4a2a'},{'ns4b2a'},{'ns5a2a'},{'ns5b2a'},{'p72a'},...
	{'core3a'},{'e13a'},{'e23a'},{'ns23a'},{'ns33a'},{'ns4a3a'},{'ns4b3a'},{'ns5a3a'},{'ns5b3a'},{'p73a'}];
QQ=sort(Q1);

%***********SET***PARAMETERS*************%

% Number of epitopes to be selected
k=27;

%Minimum number of epitopes from each antigen to be included
t_A=1;

%Antigen processing threshold (minimum proteasomal cleavage score)
t_AP = -0.28328;

% Conservation threshold
t_C = 0.9;

%Minimum number of MHC alleles to be covered
t_MHC = 27; 

%probabilities of epitope occuring in target population (from database)
p=raw(2:end,2);

%immunogenicity threshold
t_I=raw(2:end,3);

% All candidate epitopes
[num txt raw]=xlsread('ep_c_sAP.xls');
E=raw(2:end,1);

% Epitope conservation
c=raw(2:end, 2);

%Antigen processing score assigned to an epitope by the proteasomal cleavage matrix given in the supplement of Vider-Shalit et al.
s_AP=raw(2:end, 3);

%%
% Pairs of overlapping epitopes
[num txt raw]=xlsread('overlapping_ep.xls');
O=raw(:,:);
overlap=reshape(O,5809*5,1);

%immunogenicities ##READ IN TEXT FILE
% fid=fopen('Text_S2.txt')
% ep=struct('core1a',[],'core1b',[],'core2a',[],'core3a',[],...
%     'e11a',[],'e11b',[],'e12a',[],'e13a',[],'e21a',[],'e22a',[],'e21b',[],'e23a',[],...
%     'ns21a',[],'ns21b',[],'ns22a',[],'ns23a',[],'ns31a',[],'ns31b',[],'ns32a',[],'ns33a',[],...
%     'ns4a1a',[],'ns4a1b',[],'ns4a2a',[],'ns4a3a',[],'ns4b1a',[], 'ns4b1b',[],'ns4b2a',[],'ns4b3a',[],...
%     'ns5a1a',[],'ns5a1b',[],'ns5a2a',[],'ns5a3a',[],...
%      'ns5b1a',[],'ns5b1b',[],'ns5b2a',[],'ns5b3a',[],'p71a',[],...
%     'p71b',[],'p72a',[],'p73a',[]);
% 	
% j=1;
% for j=1:length(Q)
%     for i=1:500
%         epitope=(fscanf(fid,'%s',1));
%         ep.(Q{j}){i}=mat2cell(epitope);
%         if isequal(epitope,';')
%          break
%         end
%     end
% end
% 
% fclose(fid)
load ep.mat 

%%
%Immunogenicity of an epitope with respect to an allele
fid2=fopen('i_values.txt')
% allele=struct('A_0201',[],'A_0205',[],'A_0101',[],'A_1101',[],'A_2402',[],'A_0301',[],...
%     'A_3101',[],'A_6801',[],'B_0401',[],'B_2702',[],'B_2705',[],'B_3701',[],'B_3801',[],...
%     'B_3901',[],'B_4403',[],'B_5101',[],'B_5102',[],'B_5103',[],'B_4001',[],'B_4006',[],...
%     'B_1501',[],'B_0801',[],'B_5801',[],'B_5201',[],'Cw_0401',[],'Cw_0602',[],'Cw_0702',[]);

for j=1:4461
    imm=fscanf(fid2, '%s', 28);
    pep(j).seq=imm(1:9);
end
fclose(fid2)

fid2=fopen('i_values.txt')
for j=1:4461
imm=fscanf(fid2, '%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 27);
    pep(j).im=imm;
end
fclose(fid2)

for j=1:1339
    for i=1:27
    A{i,j+1}=pep(ind(j)).im(i);
    end
end
    
for j=1:1339
    f(j)=-sum(cell2mat(p(:))' * cell2mat(A(:,j+1)));
end

% B is a 27row (alleles) by 1340 column (epitopes) matrix that shows
% immunogenicity of the epitope at that specific allele. this is I

%%

%maximize Overall_Immunogenicity: 

clear f j k i 

AP=cell2mat(s_AP);

j=0
for i=1:4461
    % CONSTRAINT: SUBJECT TO ANTIGEN PROCESSING THRESHOLD
    if AP(i)<t_AP %select only epitopes with a proteasomal cleavage score greater or equal to t_AP
        continue;
    end
    j=j+1
    ind(j)=i; %indexes all peptides of 4461 that meet criteria
end

%% Set of epitopes, which when bound to an MHC allele a,
% display an immunogenicity greater than or equal to a given
% threshold t_I[a]
% set I_allele { a in A } := { e in E: i[e,a] >= t_I[a] }; 
%sort through all the epitopes and find the ones that are equal or larger
%to the t_I for the allele they bind to

for i=1:27
    for j=2:1340
        if A{i,j}<=t_I{i};
            A{i,j}=0;
        end
    end
end
        
%% 
%further constraints
 
C=1 - cell2mat(c); %CONSTRAINT: EPITOPE CONSERVATION
E_e=ones(4461,1) %CONSTRAINT: NUMBER OF EPITOPES =27
%%
%CONSTRAINT: EPITOPE OVERLAPPING
ov=cell(length(overlap)-2,2);
j=0;
for i=1:(length(overlap)-1)
    if isnan(overlap{i})
        continue
    end
    j=j+1;
    ov{j,1}=cellstr(overlap{i}(2:10));
    ov{j,2}=cellstr(overlap{i}(13:21));
end

ovlp=[1;1];
h = waitbar(0,'Please wait...');
for i=1:length(ov)
    for j=1:length(E)
        if strcmpi(ov{i,1},E{j}(2:end))
            ovlp(i,1)=j;
        end
        if strcmpi(ov{i,2},E{j}(2:end))
            ovlp(i,2)=j;
        end
    end
    waitbar(i/length(ov))
end
close(h)      

%only include those who meet antigen processing threshold
ovlap=[];
j=0;
for i=1:length(ovlp)
        if ismember(ovlp(i,1),ind) || ismember (ovlp(i,2),ind)
            j=j+1;
            ovlap(j,:)=ovlp(i,:);
        else
            continue
        end
end
%%
%CONSTRAINT: if there are overlapping epitopes, make sure only one is
%selected in the final group of 27
        

for i=1:length(ovlap) % get rid of repeats (replace index with 0)
    if ismember(ovlap(i,1), ovlap(:,2))
        ovlap(i,1)=0
    end    
end
    
nonov=ovlap(:,1);
nonov=nonov(nonov~=0);
nonov2=ovlap(:,2);
nonov2=unique(nonov2);
nonoverlap=unique([nonov;nonov2]); %now down to 4205 epitopes

%which of these epitopes are contained in the 1339 we are working with?
j=0;
for i=1:length(ind)
    if ~ismember(ind(i), nonoverlap)
        j=j+1;
        getridof(j)=ind(i)
    end
end

for i=1:5 %now get rid of these in the 1339
    f(getridof(i))=0;
end

%left with 1334...well that didn't do much
%% ILP
%f(j)=-27*sum(cell2mat(p)' * cell2mat(A(:,2:end)));

[x fval]=bintprog(f(1:1339)',[C(ind)'], [1-t_C], [E_e(ind)'],[27])

k=0;
for i=1:1339
if x(i)~=0;
    k=k+1;
match(k)=ind(i);
end
end
    
for i=1:27
    index=match(i);
    disp(E(index))
end
