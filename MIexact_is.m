for i = 1:size(x_dev,1) %estimate magnitude of deviations
    for j = 1:size(x_dev,2)
        R(i,j) = sqrt(x_dev(i,j)^2 + y_dev(i,j)^2 + z_dev(i,j)^2)^0.5;
    end
end

for i = 1:320
    for j = i:320
        M(i,j) = MI(R(:,i),R(:,j));
    end
end

M = M + M' - diag(diag(M));