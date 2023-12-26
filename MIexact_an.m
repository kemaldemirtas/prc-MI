% load('singlet MIan.mat')

for x = 1:320
    for y = x:320
        X1 = x_dev(:,x);
        Y1 = y_dev(:,x);
        Z1 = z_dev(:,x);
        X2 = x_dev(:,y);
        Y2 = y_dev(:,y);
        Z2 = z_dev(:,y);
        
%         X1 = [0.4917; 0.5165; 0.1946; 0.4962; 0.7382];
%         Y1 = [0.5516; 0.6046; 0.7160; 0.6835; 0.7534];
%         Z1 = [0.2060; 0.5779; 0.4144; 0.8062; 0.3452];
%         X2 = [0.9470; 0.6101; 0.5458; 0.7516; 0.0436];
%         Y2 = [0.6173; 0.4129; 0.1417; 0.1948; 0.1009];
%         Z2 = [0.5858; 0.0656; 0.7150; 0.4471; 0.0162];
        
        t = 6; %number of bins
        
        F = cell(t,1);
        E = cell(t,1);
        D = cell(t,1);
        C = cell(t,1);
        B = cell(t,1);
        A = cell(t,1);
        for i = 1:t
            F{i} = 0;
        end
        for i = 1:t
            E{i} = F;
        end
        for i = 1:t
            D{i} = E;
        end
        for i = 1:t
            C{i} = D;
        end
        for i = 1:t
            B{i} = C;
        end
        for i = 1:t
            A{i} = B;
        end
        
        int_X1 = linspace(min(X1),max(X1),t+1);
        int_Y1 = linspace(min(Y1),max(Y1),t+1);
        int_Z1 = linspace(min(Z1),max(Z1),t+1);
        int_X2 = linspace(min(X2),max(X2),t+1);
        int_Y2 = linspace(min(Y2),max(Y2),t+1);
        int_Z2 = linspace(min(Z2),max(Z2),t+1);
        
        for i = 1:size(X1,1)
            n = 0;
            for a = 1:t
                if X1(i) >= int_X1(a)-10^-9 && X1(i) <= int_X1(a+1)
                    for b = 1:t
                        if Y1(i) >= int_Y1(b)-10^-9 && Y1(i) <= int_Y1(b+1)
                            for c = 1:t
                                if Z1(i) >= int_Z1(c)-10^-9 && Z1(i) <= int_Z1(c+1)
                                    for d = 1:t
                                        if X2(i) >= int_X2(d)-10^-9 && X2(i) <= int_X2(d+1)
                                            for e = 1:t
                                                if Y2(i) >= int_Y2(e)-10^-9 && Y2(i) <= int_Y2(e+1)
                                                    for f = 1:t
                                                        if Z2(i) >= int_Z2(f)-10^-9 && Z2(i) <= int_Z2(f+1)
                                                            A{a}{b}{c}{d}{e}{f} = A{a}{b}{c}{d}{e}{f} + 1;
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

        for a = 1:t
            for b = 1:t
                for c = 1:t
                    for d = 1:t
                        for e = 1:t
                            for f = 1:t
                                P{a}{b}{c}{d}{e}{f} = A{a}{b}{c}{d}{e}{f} / size(X1,1);
                            end
                        end
                    end
                end
            end
        end
        
        Ps1 = K{x}{2};
        Ps2 = K{y}{2};
                                                                                                                                                                                                                
        sumMI = 0;
        for a = 1:t
            for b = 1:t
                for c = 1:t
                    for d = 1:t
                        for e = 1:t
                            for f = 1:t
                                if P{a}{b}{c}{d}{e}{f} ~= 0 && Ps1{a}{b}{c} ~= 0 && Ps2{d}{e}{f} ~= 0
                                    sumMI = sumMI + P{a}{b}{c}{d}{e}{f} * log(P{a}{b}{c}{d}{e}{f} / (Ps1{a}{b}{c} * Ps2{d}{e}{f}));
                                end
                            end
                        end
                    end
                end
            end
        end
        MI(x,y) = sumMI;

%         T{x}{y} = A;
%         sumA = 0;
%         for r = 1:t
%             for j = 1:t
%                 for k = 1:t
%                     for l = 1:t
%                         for m = 1:t
%                             for n = 1:t
%                                 sumA = A{r}{j}{k}{l}{m}{n} + sumA;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
    end
end
MIex = MI + MI' - diag(diag(MI));
MI = MIex;