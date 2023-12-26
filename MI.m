    function M = MI(OP1,OP2)

% w1 = (max(OP1) - min(OP1)) / (1 + log2(length(OP1)));
% w2 = (max(OP2) - min(OP2)) / (1 + log2(length(OP2)));
% 
% % w1 = 3.49 * std(OP1) * length(OP1)^(-1/3);
% % w2 = 3.49 * std(OP2) * length(OP2)^(-1/3);
% % 
% int_OP1 = min(OP1):w1:max(OP1);
% if max(int_OP1) ~= max(OP1)
%    int_OP1(1,end+1) = max(OP1); 
% end
% t1 = (length(int_OP1) - 1);
% int_OP2 = min(OP2):w2:max(OP2);
% if max(int_OP2) ~= max(OP2)
%    int_OP2(1,end+1) = max(OP2); 
% end
% t2 = (length(int_OP2) - 1);
% clear int_OP1 int_OP2
t1 = 10;
t2 = 10;
int_OP1 = linspace(min(OP1),max(OP1),t1+1);
int_OP2 = linspace(min(OP2),max(OP2),t2+1);

%% i0j0
A4 = histcounts2(OP2,OP1,linspace(min(OP2),max(OP2),t2+1),linspace(min(OP1),max(OP1),t1+1));
p4 = A4 ./ sum(sum(A4));

%% i0
A2 = hist(OP1,t1);
p2 = A2 ./ sum(A2);

%% j0
A3 = hist(OP2,t2);
p3 = A3 ./ sum(A3);

sumM = 0;
for i = 1:t1
   for j = 1:t2
      if p4(j,i) ~= 0 
         M(j,i) = p4(j,i) * log(p4(j,i) / (p2(i) * p3(j)));
%          sumM = sumM + M;
      end
   end
end
M = sum(sum(M)) + ((nnz(p4) - nnz(p2) - nnz(p3) + 1) / (2*length(OP1)));
end