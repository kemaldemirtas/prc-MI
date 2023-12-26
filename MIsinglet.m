N = 320; %number of residues

for x = 1:320
    X1 = x_dev(:,x);
    Y1 = y_dev(:,x);
    Z1 = z_dev(:,x);
    
    t = 6; %number of bin

    int_X1 = linspace(min(X1),max(X1),t+1);
    int_Y1 = linspace(min(Y1),max(Y1),t+1);
    int_Z1 = linspace(min(Z1),max(Z1),t+1);
    A = hist(Z1,t);

    F = cell(t,1);
    E = cell(t,1);
    D = cell(t,1);
    Ps = cell(t,1);
    
    for i = 1:t
        F{i} = 0;
    end
    for i = 1:t
        E{i} = F;
    end
    for i = 1:t
        D{i} = E;
    end
    
    for i = 1:size(X1,1)
        n = 0;
        for d = 1:t
            if X1(i) >= int_X1(d)-10^-9 && X1(i) <= int_X1(d+1)
                for e = 1:t
                    if Y1(i) >= int_Y1(e)-10^-9 && Y1(i) <= int_Y1(e+1)
                        for f = 1:t
                            if Z1(i) >= int_Z1(f)-10^-9 && Z1(i) <= int_Z1(f+1)
                                D{d}{e}{f} = D{d}{e}{f} + 1;
                                n = 1;
                            end
                            if n == 1
                                break
                            end
                        end
                    end
                    if n == 1
                        break
                    end
                end
            end
            if n == 1
                break
            end 
        end
    end
    
    for i = 1:t
        for j = 1:t
            for k = 1:t
                Ps{i}{j}{k} = D{i}{j}{k} / size(X1,1);
            end
        end
    end
    K{x}{1} = D;
    K{x}{2} = Ps;
end