%% main
function main()
    close all;
    clc;
    addpath(genpath(pwd));
    
    %% initalization
    T = 0.05;
    F=[ 1,  T,  T*T/2.,  0,  0,  0;          ...
        0,  1,  T,       0,  0,  0;         ...
        0,  0,  1,       0,  0,  0;         ...
        0,  0,  0,       1,  T,  T*T/2.;  ...
        0,  0,  0,       0,  1,  T;        ...
        0,  0,  0,       0,  0,  1];

    H = [1  0   0   0   0   0;                          % 测量模型
         0  0   0   1   0   0];

    Q = [0.02    0       0       0       0       0;      % 过程噪声协方差矩阵
         0      0.002     0       0       0       0;
         0      0       0.0002    0       0       0;
         0      0       0       0.02    0       0;
         0      0       0       0       0.002   0;
         0      0       0       0       0       0.0002];

    R = [1.5  0;                                      % 测量噪声协方差矩阵
         0  1.5];

    P = [10  0   0   0   0   0;                          % 初始误差协方差矩阵
         0  10   0   0   0   0;
         0  0   10   0   0   0;
         0  0   0   10   0   0;
         0  0   0   0   10   0;
         0  0   0   0   0   10];

    I = eye(size(F,1));
    N = size(Q,1);                                          % 状态维度
    M = size(R,1);                                          % 测量维度

    %% Monte Carlo
    monteCarlo.count = 10000;
    monteCarlo.means = zeros(4, monteCarlo.count);
    monteCarlo.RMSEs = zeros(4, monteCarlo.count);
    monteCarlo.time  = zeros(4, monteCarlo.count);
    for mc=1:monteCarlo.count
       %% get measurements(with noise)
        [time, Z] = S_Trajectory();
        n = length(time);
        Y0(1,:) = Z(1,:) .* sin(Z(2,:));
        Y0(2,:) = Z(1,:) .* cos(Z(2,:));
        if(1)
            Y(1,:) = Y0(1,:) + normrnd(0, R(1,1), 1, n);
            Y(2,:) = Y0(2,:) + normrnd(0, R(2,2), 1, n);
            save('Y');
        else
            raw = load('Y.mat');
            Y = raw.Y;
        end
        
        X0 = [Y(1,1) 0 0 Y(2,1) 0 0 ]';  % 用测量值初始化状态值

       %% Singer-KF 
        a = 1/10;
        F_sg = [ 1  T  (a*T-1+exp(-a*T)/a*a) 0  0  0;...
                 0  1  (1-exp(-a*T)/a)       0  0  0;...
                 0  0  exp(-a*T)             0  0  0;...
                 0  0  0                     1  T  (a*T-1+exp(-a*T)/a*a);...
                 0  0  0                     0  1  (1-exp(-a*T)/a);...
                 0  0  0                     0  0  exp(-a*T); ];
        X_sg = X0;
        P_sg = P;
        Q_sg = Q;
        R_sg = R;
        XX_sg = zeros(N,n);
        PP_sg = zeros(N,N,n);
        t1 = clock;
        for i = 1 : n
            Xp_sg = F_sg * X_sg;
            Pp_sg = F_sg * P_sg * F_sg' + Q_sg;
            innovY = Y(:,i) - H * Xp_sg;
            S_sg = H * Pp_sg * H' + R_sg;
            K_sg = Pp_sg * H' / S_sg;
            X_sg = Xp_sg + K_sg * innovY;
            P_sg = (I - K_sg * H) * Pp_sg;
            XX_sg(:,i) = X_sg;
            PP_sg(:,:,i) = P_sg;
        end
        time_sg = etime(clock, t1);

       %% CA-KF
        F_ca = [ 1,  T,  T*T/2.,  0,  0,  0;          ...
                 0,  1,  T,       0,  0,  0;         ...
                 0,  0,  1,       0,  0,  0;         ...
                 0,  0,  0,       1,  T,  T*T/2.;  ...
                 0,  0,  0,       0,  1,  T;        ...
                 0,  0,  0,       0,  0,  1];
        X_ca = X0;
        P_ca = P;
        Q_ca = 3*Q;
        R_ca = R;
        XX_ca = zeros(N,n);
        PP_ca = zeros(N,N,n);
        t1 = clock;
        for i = 1 : n
            Xp_ca = F_ca * X_ca;
            Pp_ca = F_ca * P_ca * F' + Q_ca;
            innovY = Y(:,i) - H * Xp_ca;
            S_ca = H * Pp_ca * H' + R_ca;
            K_ca = Pp_ca * H' / S_ca;
            X_ca = Xp_ca + K_ca * innovY;
            P_ca = (I - K_ca * H) * Pp_ca; 
            XX_ca(:,i) = X_ca;
            PP_ca(:,:,i) = P_ca;
        end
        time_ca = etime(clock, t1);

       %% SVR-KF
%         % 训练集为ds对角列向量
%         F_svr = [ 1,  T,  T*T/2.,  0,  0,  0;          ...
%                 0,  1,  T,       0,  0,  0;         ...
%                 0,  0,  1,       0,  0,  0;         ...
%                 0,  0,  0,       1,  T,  T*T/2.;  ...
%                 0,  0,  0,       0,  1,  T;        ...
%                 0,  0,  0,       0,  0,  1];
%         X_svr = X0;
%         P_svr = P;
%         Q_svr = Q;
%         R_svr = R;
%         XX_svr = zeros(N,n);
%         PP_svr = zeros(N,N,n);
%         winSize = 5;
%         innovationY = zeros(M, winSize);
%         innovationX = zeros(N, winSize);
%         train = zeros(M, winSize)';
%         label = zeros(1, winSize)';
%         qk = 1;
%         rk = 1; 
%         qqk = zeros(n,1);
%         dds = zeros(n,2);
%         t1 = clock;
%         for i = 1 : n
%             qqk(i) = qk;
%             Qa = qk*Q_svr;
%             Ra = rk*R_svr;
%             Xp_svr = F_svr * X_svr;
%             Pp_svr = F_svr * P_svr * F' + Qa;
%             innovY = Y(:,i) - H * Xp_svr;
%             S_svr = H * P_svr * H' + Ra;
%             K_svr = Pp_svr * H' / S_svr;
%             X_svr = Xp_svr + K_svr * innovY;
%             P_svr = (I - K_svr * H) * Pp_svr;
%             XX_svr(:,i) = X_svr;
%             PP_svr(:,:,i) = P_svr;
% 
%             if(i > winSize)
%                 S1 = innovationY * innovationY' / winSize;
%                 ds = diag(S1 - S_svr)';
%                 dds(i,:) = ds;
%                 if(i > 2*winSize)
%                     model = svmtrain(label, train,'-s 3 -t 0 -c 15.2 -g 0.01 -p 0.01');
%                     qk= svmpredict(0, ds, model);
%                     label(mod(i-1,winSize)+1,:) = qk;
%                 else
%                     P1 = innovationX * innovationX' / winSize;
%                     dp = diag(P1 - P_svr)';
%                     label(mod(i-1,winSize)+1,:) = 2*sqrt(dp*dp');
%                 end
%                 train(mod(i-1,winSize)+1,:) = ds;
%             end
%             innovationX(:,mod(i-1, winSize)+1) = X_svr - Xp_svr;
%             innovationY(:,mod(i-1, winSize)+1) = innovY;
%         end
%         time_svr = etime(clock, t1);
%         if(monteCarlo.count==1)
%             figure;
%             plot(qqk);
%             hold on;
%             plot(dds(:,1));
%             grid on;
%         end
        
       %% SVR-KF-1
    %     F_svr = [ 1,  T,  T*T/2.,  0,  0,  0;          ...
    %             0,  1,  T,       0,  0,  0;         ...
    %             0,  0,  1,       0,  0,  0;         ...
    %             0,  0,  0,       1,  T,  T*T/2.;  ...
    %             0,  0,  0,       0,  1,  T;        ...
    %             0,  0,  0,       0,  0,  1];
    %     X_svr = X0;
    %     P_svr = P;
    %     Q_svr = Q;
    %     R_svr = R;
    %     XX_svr = zeros(N,n);
    %     PP_svr = zeros(N,N,n);
    %     winSize = 20;
    %     innovationY = zeros(M, winSize);
    %     innovationX = zeros(N, winSize);
    %     train = zeros(1, winSize)';
    %     label = zeros(1, winSize)';
    %     for i = 1:winSize
    %         train(i) = i*5;
    %         label(i) = i;
    %     end
    %     qk = 1;
    %     rk = 1;
    %     qqk = label;
    %     dds = train;
    %     t1 = clock;
    %     for i = 1 : n
    %         Qa = qk*Q_svr;
    %         Ra = rk*R_svr;
    %         Xp_svr = F_svr * X_svr;
    %         Pp_svr = F_svr * P_svr * F' + Qa;
    %         innovY = Y(:,i) - H * Xp_svr;
    %         S_svr = H * P_svr * H' + Ra;
    %         K_svr = Pp_svr * H' / S_svr;
    %         X_svr = Xp_svr + K_svr * innovY;
    %         P_svr = (I - K_svr * H) * Pp_svr;
    %         XX_svr(:,i) = X_svr;
    %         PP_svr(:,:,i) = P_svr;
    %         
    %         if (i > winSize)
    %             model = svmtrain(label, train,'-s 3 -t 2 -c 5.2 -g 1 -p 0.01');
    %             S1 = innovationY * innovationY' / winSize;
    %             ds = norm(diag(S1-S_svr), 'fro');
    %             qk = svmpredict(0, ds, model);
    %             train(mod(i-1,winSize)+1,:) = ds;
    %             label(mod(i-1,winSize)+1,:) = qk;
    %             dds = [dds; ds];
    %             qqk = [qqk; qk];
    %         end
    %         innovationY(:,mod(i-1, winSize)+1) = innovY;
    %     end
    %     time_svr = etime(clock, t1);
    % %     figure;
    % %     subplot(2,1,1);
    % %     plot(dds(:,1));
    % %     title('新息协方差误差X分量')
    % %     grid on;
    % %     axis([0 600 -50 100]);
    % %     subplot(2,1,2);
    % %     plot(dds(:,2));
    % %     title('新息协方差误差Y分量')
    % %     grid on;
    % %     axis([0 600 -20 40]);
    %     figure;
    %     plot(qqk);
    %     hold on;
    %     plot(dds);
    %     grid on;
    
       %% SVR-KF-2
        % 训练集为ds F范数
        F_svr = [ 1,  T,  T*T/2.,  0,  0,  0;          ...
                0,  1,  T,       0,  0,  0;         ...
                0,  0,  1,       0,  0,  0;         ...
                0,  0,  0,       1,  T,  T*T/2.;  ...
                0,  0,  0,       0,  1,  T;        ...
                0,  0,  0,       0,  0,  1];
        X_svr = X0;
        P_svr = P;
        Q_svr = Q;
        R_svr = R;
        XX_svr = zeros(N,n);
        PP_svr = zeros(N,N,n);
        winSize = 5;
        innovationY = zeros(M, winSize);
        innovationX = zeros(N, winSize);
        train = zeros(winSize,1);
        label = zeros(winSize,1);
        qk = 1;
        qqk = zeros(n,1);
        dds = zeros(n,1);
        t1 = clock;
        for i = 1 : n
            Qa = qk*Q_svr;
            Ra = R_svr;
            Xp_svr = F_svr * X_svr;
            Pp_svr = F_svr * P_svr * F' + Qa;
            innovY = Y(:,i) - H * Xp_svr;
            S_svr = H * P_svr * H' + Ra;
            K_svr = Pp_svr * H' / S_svr;
            X_svr = Xp_svr + K_svr * innovY;
            P_svr = (I - K_svr * H) * Pp_svr;
            XX_svr(:,i) = X_svr;
            PP_svr(:,:,i) = P_svr;

            if(i > winSize)
                S1 = innovationY * innovationY' / winSize;
                ds = norm(S1 - S_svr,'fro');
                dds(i) = ds;
                qqk(i) = qk;
                if(i > 2*winSize)
                    model = svmtrain(label, train,'-s 3 -t 0 -p 0.01');
                    qk0= svmpredict(0, ds, model);
                    if qk0 < 40 && qk0>0
                        qk = qk0;
                    end
                    label(mod(i-1,winSize)+1,:) = qk;
                else
                    P1 = innovationX * innovationX' / winSize;
                    dp = norm(P1 - P_svr,'fro');
                    label(mod(i-1,winSize)+1,:) = 2*dp;
                    qqk(i) = dp;
                end
                train(mod(i-1,winSize)+1,:) = ds;
            end
            innovationX(:,mod(i-1, winSize)+1) = X_svr - Xp_svr;
            innovationY(:,mod(i-1, winSize)+1) = innovY;
        end
        time_svr = etime(clock, t1);
        if(monteCarlo.count==1)
            figure;
            plot(qqk);
            hold on;
            plot(dds(:,1));
            grid on;
        end
        
       %% 误差统计
        meanPos_Y   = mean( sqrt((Y0(1,:)-Y(1,:)).^2 + (Y0(2,:)-Y(2,:)).^2));
        meanPos_CA  = mean( sqrt((Y0(1,:)-XX_ca(1,:)).^2 + (Y0(2,:)-XX_ca(4,:)).^2));
        meanPos_SG  = mean( sqrt((Y0(1,:)-XX_sg(1,:)).^2 + (Y0(2,:)-XX_sg(4,:)).^2));
        meanPos_SVR = mean( sqrt((Y0(1,:)-XX_svr(1,:)).^2 + (Y0(2,:)-XX_svr(4,:)).^2));
        rmsePos_Y    = var( sqrt((Y0(1,:)-Y(1,:)).^2 + (Y0(2,:)-Y(2,:)).^2));
        rmsePos_CA   = var( sqrt((Y0(1,:)-XX_ca(1,:)).^2 + (Y0(2,:)-XX_ca(4,:)).^2));
        rmsePos_SG   = var( sqrt((Y0(1,:)-XX_sg(1,:)).^2 + (Y0(2,:)-XX_sg(4,:)).^2));
        rmsePos_SVR  = var( sqrt((Y0(1,:)-XX_svr(1,:)).^2 + (Y0(2,:)-XX_svr(4,:)).^2));
        monteCarlo.means(:, mc) = [meanPos_Y; meanPos_CA; meanPos_SG; meanPos_SVR];
        monteCarlo.RMSEs(:, mc) = [rmsePos_Y; rmsePos_CA; rmsePos_SG; rmsePos_SVR];
        monteCarlo.time(:, mc) = [0; time_ca; time_sg; time_svr];
    end
    error.means = zeros(size(monteCarlo.means,1),1);
    error.RMSE = zeros(size(monteCarlo.RMSEs,1),1);
    error.time = zeros(size(monteCarlo.time,1),1);
    if monteCarlo.count == 1
        error.means = monteCarlo.means;
        error.RMSE = monteCarlo.RMSEs;
        error.time = monteCarlo.time;
    else
        for i=1:size(monteCarlo.means,1)
            error.means(i) = mean(monteCarlo.means(i,:));
            error.RMSE(i) = mean(monteCarlo.RMSEs(i,:));
            error.time(i) = mean(monteCarlo.time(i,:));
        end
    end
    
    %% print
    fprintf('\n========================================\n');
    fprintf('Monte Carlo Simulation Count: %d\n', monteCarlo.count);
    fprintf('          mean(m) \t RMSE(m) \t time(sec)\n');
    fprintf('Meausre:  %+.4f \t %+.4f \t %+.4f\n', error.means(1),error.RMSE(1),error.time(1));
    fprintf('CA:       %+.4f \t %+.4f \t %+.4f\n', error.means(2),error.RMSE(2),error.time(2));
    fprintf('Singer:   %+.4f \t %+.4f \t %+.4f\n', error.means(3),error.RMSE(3),error.time(3));
    fprintf('SVR:      %+.4f \t %+.4f \t %+.4f\n', error.means(4),error.RMSE(4),error.time(4));
    fprintf('========================================\n\n');
    
    %% plot
    color = {[0,0,0],...
             [0.4,0.4,0.4],...
             [0.467,0.674,0.188],...
             [0,0.353,0.749],...
             [1,0,0]};
    figure;
    set(gcf,'color','white')
    % Monte Carlo means
    subplot(2,1,1)
    plot(monteCarlo.means(1,:), '--', 'color', color{1}); hold on;
    plot(monteCarlo.means(2,:), '-', 'color', color{3}); hold on;
    plot(monteCarlo.means(3,:), '-.', 'color', color{4}); hold on;
    plot(monteCarlo.means(4,:), ':','color',color{5},'linewidth',1.5); hold on;
    xlabel('tick');
    ylabel('Position Mean/(m)');
    grid on;
    title('1,0000 Monte Carlo Simulation Mean in scenario 1');
    legend('Measurement','CA','Singer','SVR-CA');
    
    % Monte Carlo RMSE
    subplot(2,1,2);
    plot(monteCarlo.RMSEs(1,:), '--', 'color', color{1}); hold on;
    plot(monteCarlo.RMSEs(2,:), '-', 'color', color{3}); hold on;
    plot(monteCarlo.RMSEs(3,:), '-.', 'color', color{4}); hold on;
    plot(monteCarlo.RMSEs(4,:), ':','color',color{5},'linewidth',1.5); hold on;
    xlabel('tick');
    ylabel('Position RMSE/(m)');
    grid on;
    title('1,0000 Monte Carlo Simulation RMSE  in scenario 1');
    legend('Measurement','CA','Singer','SVR-CA');
    
    % Trajectory
    figure;
    set(gcf,'color','white')
    plot(Y0(1,:),Y0(2,:), '-', 'color', color{1}, 'linewidth',1.3); hold on;
    plot(Y(1,:),Y(2,:),   '*',  'color', color{2},   'linewidth',0.2, 'markersize',5); hold on;
    plot(XX_ca(1,:),XX_ca(4,:), '-.', 'color', color{3} ,'linewidth',1.5); hold on;
    plot(XX_sg(1,:),XX_sg(4,:), '.-', 'color', color{4}, 'linewidth',1); hold on;
    plot(XX_svr(1,:),XX_svr(4,:), '-','color',color{5}, 'linewidth',1.2);
    xlabel('x/m');
    ylabel('y/m');
    grid on, axis equal;
    axis([-10, 60, 5, 95]);
    title('Trajectory in scenario 1');
    legend('Theoretical','Noisy','CA','Singer','SVR-CA');
end


