% clc
% clear all


N = size(x_dev,2); %number of amino acid
t = size(x_dev,1); %number of frames

j = 1; % R is n tx3N matrix. Each i:i+2 column where rem(i,3)=1 corresponds to fluctuations of xyz coordinates of ith residue.
for i = 1:3:3*N
    R(:,i:i+2) = [x_dev(:,j) y_dev(:,j) z_dev(:,j)];
    j = j + 1;
end

%% Anisotropic
%Calculation of det<deltaRdeltaRT>:
%Forms a tx6 matrix containing fluctuations of residue i and j. Multiply
%its transpose with itself, which corresponds to multiplication of each row
%with itself that forms t 6x6 matrices and summing all resulting matrices.
%Take its time average. Find its determinant.
for i = 1:3:3*N
    for j = i:3:3*N
        C = [R(:,i:i+2) R(:,j:j+2)]' * [R(:,i:i+2) R(:,j:j+2)];
        avgC(ceil(i/3),ceil(j/3)) = det(C/t);
    end
end
%Finally, sum upper and lower triangular matrices. Each entry of diagonal of avgC
%resulting from the loop above must be equal to zero but it doesn't happen
%probably due to matlab's limited capacity for storing number of decimal
%point. (-diag(diag(avgC))) term serve to make the diagonal zero.
avgC = avgC + avgC' - diag(diag(avgC)); 

%Calculation of det<deltaRideltaRiT>:
%Do the same thing as above for i and i.
j = 1;
for i = 1:3:3*N
    sumCi = 0;
    for k = 1:t
        Ci = [R(k,i:i+2)]' * [R(k,i:i+2)];
        sumCi = sumCi + Ci;
    end
    avgCi(j) = det(sumCi/t);
    j = j + 1;
end

%Calculation of mutual information:
for i = 1:N
    for j = i:N
        MI(i,j) = -0.5 * log(avgC(i,j)/(avgCi(i)*avgCi(j)));
    end
end
MIan = MI + MI'; 

%% Isotropic
% Take dot product of transposes of matrices of residue i and j, that
% equals to <deltaRi . deltaTj>
for i = 1:3:3*N
    for j = i:3:3*N
        RidotRj = dot(R(:,i:i+2)',R(:,j:j+2)');
        avgR(ceil(i/3),ceil(j/3)) = mean(RidotRj);
    end
end
avgR = avgR + avgR' - diag(diag(avgR)); % (diag(diag(avgR))) terms is required for the same reason as above.

%Calculation of mutual information:
for i = 1:N
    for j = i:N
        MIiso(i,j) = -0.5 * log(1 - (avgR(i,j)^2/(avgR(i,i)*avgR(j,j))));
    end
end
MIiso = MIiso + MIiso';
