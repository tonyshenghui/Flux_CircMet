% calculate direct contributions (and the asscoiated errors) to a
% circulating nutrient from all other circulating nutrients;
%
% values of direct contribution are constrained to be non-negative
%
% use bootstrapping to estimate errors
%
% Hui et al.


input_file = 'normalizedLabelingData_fastedCD.xlsx'; %see the example input file for format
output_file = 'result_directContributions_fastedCD.xlsx';


[values,text,raw] = xlsread(input_file);

ind_nan = find(isnan(values(:,1)));

M = values(1:ind_nan(1)-1,1:2:end);%normalized labeling of circulating nutrients
dM = values(1:ind_nan(1)-1,2:2:end);%errors for normalized labeling of circulating nutrients
M = M';
dM = dM';

matrix = M;
dmatrix = dM;
N = size(matrix,1);

% number of rounds of bootstrapping
n_iter = 100;

%normalized labeling (and the errors) of tissue metabolite
metabolite = values(1:ind_nan(1)-1,:);
output = zeros(size(metabolite));

for i = 1:size(metabolite,1)%for each tissue
    
    zind = setdiff(1:N,i);
    M = matrix(zind,zind);
    dM = dmatrix(zind,zind);
    
    metabolite_i = metabolite;
    metabolite_i(:,i*2-1:i*2) = [];
    
    output_i = NaN(n_iter,size(metabolite_i,2)/2);
    Y_i = metabolite_i(i,1:2:end);
    dY_i = metabolite_i(i,2:2:end);
    for j = 1:n_iter
        Y_ij = Y_i + randn(size(dY_i)).*dY_i;
        M_ij = M + randn(size(dM)).*dM;
        
        %optimization
        fun=@(X)myvariance(X,M_ij,Y_ij');
        lb=zeros(N-1,1);
        ub=ones(N-1,1);
        X0 = M_ij\Y_ij';
        X = fmincon(fun,X0,[],[],[],[],lb,ub);
        
        output_i(j,:) = X';
    end

    output(i,zind*2-1) = mean(output_i);
    output(i,zind*2) = std(output_i);
end

%output results
ind_empty = find(cellfun('isempty',text(:,1)));
tissues = text(4:ind_empty(end)-1,1);

nutrients = strrep(text(1,:),'infusion','');
nutrients = strrep(nutrients,'Infusion','');
nutrients = strtrim(nutrients);

output = [tissues num2cell(output)];
output = [nutrients;output];
writetable(table(output),output_file,'WriteVariableNames',false);

