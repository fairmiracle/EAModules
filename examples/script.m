addpath('../EAModules/utils/')
addpath('../EAModules/functions/')
addpath('../data/')
load('GSE35103FromString.mat')
%addpath('../../jActivModule/')
% 
% load('GSE35103FromBioGRID.mat');
% randomscore = zeros(length(G),2);
% iter_monte = 10000;
% tic;
% for i = 1:length(G)
%     [mu,sigma] = monte(G,i,array_basic_z,iter_monte);
%     randomscore(i,1) =  mu;
%     randomscore(i,2) =  sigma;
% end
% toc
% save GSE35103FromBioGRID G GENE array_basic_z randomscore comps;
% 
% popsize = 100;
% crossrate = 0.9;
% lsrate = 0.5;
% iteration = 10000;
% lsiter = 10;
% minmodulesize = 100;
% maxmodulesize = 200;
% % T_start = 1;
% % T_end = 0.001;
% %[corrected_subnet_score_SA, fsubset_SA,func_SA] = SA(G, array_basic_z, randomscore,minmodulesize,maxmodulesize,T_start,T_end,iteration);
% 
% for i = 1:10
%     tic;
%     [corrected_subnet_score_MA, fsubset_MA,func_MA] = MA(G, array_basic_z, randomscore,minmodulesize,maxmodulesize,popsize,crossrate,lsrate,lsiter,iteration);
%     toc
%     fid = fopen('resultGSE35103BioGRID_MX.txt', 'a+');
%     fprintf(fid, '%f \n', corrected_subnet_score_MA);
%     fprintf(fid, '%s \n', num2str(fsubset_MA));
%     for j = 1:length(fsubset_MA)
%         fprintf(fid, '%s\n', GENE{fsubset_MA(j)});
%     end
%     fprintf(fid, '\n\n');
%     fclose(fid);
%     
%     file_name = ['func_MA' num2str(i)];
%     save(file_name,'func_MA');
%     
%     G(fsubset_MA,:)=[];
%     G(:,fsubset_MA)=[];
%     array_basic_z(fsubset_MA)=[];
%     randomscore(fsubset_MA,:)=[];
%     GENE(fsubset_MA) = [];
% end
% 
% clear all
% addpath('../EAModules/utils/')
% addpath('../EAModules/functions/')
% addpath('../data/')
% addpath('mActive/')
% load('GSE3635ToGSE5283FromBioGRID.mat');
% randomscore = zeros(length(G),2);
% iter_monte = 10000;
% tic;
% for i = 1:length(G)
%     [mu,sigma] = monte(G,i,array_basic_z,iter_monte);
%     randomscore(i,1) =  mu;
%     randomscore(i,2) =  sigma;
% end
% toc
% save GSE3635ToGSE5283FromBioGRID G GENE array_basic_z randomscore comps;

popsize = 100;
crossrate = 0.9;
lsrate = 0.5;
iteration = 1000;
lsiter = 10;
minmodulesize = 30;
maxmodulesize = 100;
mutrate = 0.5;
% T_start = 1;
% T_end = 0.001;
%[corrected_subnet_score_SA, fsubset_SA,func_SA] = SA(G, array_basic_z, randomscore,minmodulesize,maxmodulesize,T_start,T_end,iteration);

for i = 1:20
    tic;
    [corrected_subnet_score_MA, fsubset_MA,func_MA] = MA(G, array_basic_z, randomscore,minmodulesize,maxmodulesize,popsize,crossrate,lsrate,lsiter,iteration);
    toc
    fid = fopen('Result/GSE35103FromString_MA.txt', 'a+');
    fprintf(fid, '%f \n', corrected_subnet_score_MA);
    fprintf(fid, '%s \n', num2str(fsubset_MA));
    for j = 1:length(fsubset_MA)
        fprintf(fid, '%s\n', GENE{fsubset_MA(j)});
    end
    fprintf(fid, '\n\n');
    fclose(fid);
    
    file_name = ['func_MA' num2str(i)];
    save(file_name,'func_MA');
    
    G(fsubset_MA,:)=[];
    G(:,fsubset_MA)=[];
    array_basic_z(fsubset_MA)=[];
    randomscore(fsubset_MA,:)=[];
    GENE(fsubset_MA) = [];
end

load('GSE35103FromString.mat')
for i = 1:20
    [corrected_subnet_score_GA,fsubset_GA,func_GA] = GA(G, array_basic_z, randomscore,minmodulesize,maxmodulesize,popsize,crossrate,mutrate,iteration);
    
    fid = fopen('Result/GSE35103FromString_GA.txt', 'a+');
    fprintf(fid, '%f \n', corrected_subnet_score_GA);
    fprintf(fid, '%s \n', num2str(fsubset_GA));
    for j = 1:length(fsubset_GA)
        fprintf(fid, '%s\n', GENE{fsubset_GA(j)});
    end
    fprintf(fid, '\n\n');
    fclose(fid);
    
    G(fsubset_GA,:)=[];
    G(:,fsubset_GA)=[];
    array_basic_z(fsubset_GA)=[];
    randomscore(fsubset_GA,:)=[];
    GENE(fsubset_GA) = [];
end