% load GRN_fp16.mat %This loads the fixed point constraints.[]]=
% load('con1.mat')%This loads the connection matrix which consists of 1's and 0's. 1 means that a particular J_ij exists 
% and 0 means that the J_ij of interest is vanishing or in other words
% those gene classes do not interact.
% typically con1.mat is a matrix of the following form
% con1=ones(nvar,nvar)if we assume that all the interactions are allowed. 
% load('fp.mat')%This loads the matrix consisting of the fixed point conditions. The rows are the gene modules and columns are the different cell states. 
%The (i,j)-th element of the matrix 
phi =0.01; 
nfp = size(fp,2);%fp is the matrix which consists of the fixed point constraints.nfp is the number of fixed points
nvar = size(fp,1);%number of gene classes
CON=con1(:);%connection matrix 
CON1=CON;
CON1(CON1==0)=[]; %Delete all the elements i.e. the coupling constants of the connection matrix that are zero
% sym2=find(CON1>1.5);
%sym1=find(CON1>0);
%neg=find(CON1<0);
Jij=find(~CON);%the nonzero coupling constants
A = zeros(nfp*nvar+1,nvar^2);%Want to solve the linear programming problem of the form Ax (less, greater or equal B)
b = zeros(nfp*nvar+1,1);%(Setting up A and B)
for i=1:nfp
    for k=1:nvar
        p = nvar*(i-1) + k;
        if fp(k,i) == 1
            b(p) = -phi; 
            A(p,(nvar*(k-1)+1):k*nvar) = -fp(:,i)';
        elseif fp(k,i) == 0
            b(p) = phi;
            A(p,(nvar*(k-1)+1):k*nvar) = fp(:,i)';
        end
    end
end
A(:,Jij)=[];
ntrials = 10000;%Number of models that we simulate
con_num=nnz(con1);%number of non-zero coupling constants
X = zeros(con_num,ntrials);
lb = -ones(con_num,1);
% diag_ele = (nvar*(0:(nvar-1)))+(1:nvar);
% lb(diag_ele) = 0;
ub = ones(con_num,1);
lb = lb; ub = ub;
%ub(neg)=0;%I impose the the symmetric to asymmetric gene interaction to be positive and the asymmetric to asymmetric interaction to be negative
%lb(sym2)=0;

opts = optimoptions(@linprog,'Algorithm','interior-point', 'Display', 'off');%Set options to use the 'interior-point' algorithm.

for nt = 1:ntrials
    if mod(nt,500)==0
        nt
    end
    f = randi(3,1,con_num)-2;
    A(nvar*nfp+1,:) = f;
    x1 = linprog(f', A, b, [], [], lb, ub, [], opts);% As set before I use linprog function to do the solve the problem
    X(:,nt) = x1;
end
M=mean(X,2);%compute the meand and std deviation of the J_ijs
S=std(X,[],2);
C=abs(S./M);
J_1=find(CON);
X2=zeros(nvar^2,ntrials);%This matrix will consist of all the solutions for the 10000 different models. 
X2(J_1,:)=X;%J_1 consists of all the non zero coupling constants
M1=mean(X2,2);
S1=std(X2,[],2);
C1=abs((S1./M1)); %Get the mean and CVs of the distributions of the coupling constants
C1(isnan(C1))=0;
%Defines the 
tf25={'klf4';'Hes6';'Hmga1';'Tead1';'Sp5';'Baz1a';'Msx2';'Snai1';'Ciao1';'churc1';'bmp';'lif';'fgf';'pou5f1';'sox2';'atf2';'otx2';'Smarce1';'Ets2';'T';'apex1';'Hes1';'Pax6';'Xab2';'Gm13051';'Brd7';'Etv5';'Fhl1';'Hmgn2'};
meanJij = zeros(nvar,nvar);%define the initial mean, CV, and connection matrices. Connection matrix consists of the signs of the interaction
conJij = zeros(nvar,nvar);
cvJij = zeros(nvar,nvar);
fracJij = zeros(nvar,nvar);
for ai=1:nvar
for kii=1:nvar
  meanJij(kii,ai)=M1(nvar*(ai-1)+kii,1);
  cvJij(kii,ai)=C1(nvar*(ai-1)+kii,1);
end
end
conJij=sign(meanJij);%Carries the sign of the different J_ijs

plotgraph_barcodes(tf25, conJij.*(cvJij<1000), 'barcodegraph_29_allgene1.sif');
%% Creating the file to use in cytoscape for visualizing the network. 
con = conJij.*(cvJij<10);% con contains the interaction with a cv of less than 10 only. Based on the network this can be changed. 
fid2 = fopen('means29_allgene1.txt', 'w'); %should be .txt file which contains the mean and cv's of the interactions.
fprintf(fid2, 'name\t means\t cvs\t frac\n');
for i=1:nvar
    for j=1:nvar
        if con(i,j)== -1
            fprintf(fid2, '%s %s %s\t %f\t %f\t %f\n', tf25{i}, '(neg)', tf25{j}, abs(meanJij(i,j)), cvJij(i,j), fracJij(i,j));
        elseif con(i,j) == 1
            fprintf(fid2, '%s %s %s\t %f\t %f\t %f\n', tf25{i}, '(pos)', tf25{j}, abs(meanJij(i,j)), cvJij(i,j), fracJij(i,j));
        end
    end
end
fclose(fid2);




