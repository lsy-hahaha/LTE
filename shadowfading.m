function s=shadowfading(roi_x,roi_y,resolution)
%roi_x=[-250,250];
%roi_y=[-216.5064 216.5064];
%resolution=20;
roi_maximum_pixels = LTE_common_pos_to_pixel( [roi_x(2) roi_y(2)], [roi_x(1) roi_y(1)], resolution);
n_neighbors=8;
switch n_neighbors
                case 4
                    offsets = [
                        -1  -1
                        -0  -1
                        1  -1
                        -1   0
                        ];
                case 8 % offsets与注释中提到的文献有关，它体现的是文献中Figure5的产生s_n图样相对位置的坐标偏差，第一列为横坐标，第二列为纵坐标。
                       % 例如，offsets的第一行-1 -1代表s_n向左移动一格向上移动一格，则得到s_1
                    offsets = [
                        -1  -1
                        -0  -1
                        1  -1
                        -1   0
                        -1  -2
                        1  -2
                        -2  -1
                        2  -1
                        ];
                otherwise
                    error('Only 4 and 8 values supported');
end
            alpha = 1/20;
            r = @(d) exp(-alpha*d);% 文献中的公式（5）
            d = resolution; % How many meters a hop is %也就是所谓的网格分辨率

            % Map size in pixels
            N_y = roi_maximum_pixels(1);
            N_x = roi_maximum_pixels(2);
num_eNodeBs=1;
r_eNodeBs = 0.5;
mean = 0;
std = 10;
a_n_random_matrix = mean + std*randn(N_y,N_x,num_eNodeBs);
a_n_original_map = mean + std*randn(N_y,N_x);
a_n_matrix = zeros(N_y,N_x,num_eNodeBs);
                for i_=1:num_eNodeBs
                    % Generate i gaussian maps based on an 'original' one
                    a_n_matrix(:,:,i_) = sqrt(r_eNodeBs)*a_n_original_map + sqrt(1-r_eNodeBs)*a_n_random_matrix(:,:,i_);% ？？？
                end
switch n_neighbors %参见文献的conclusions部分
                case 4
                    R = [
                        1              r(d)           r(2*d)         r(d)           r(sqrt(2)*d)
                        r(d)           1              r(d)           r(sqrt(2)*d)   r(d)
                        r(2*d)         r(d)           1              r(sqrt(5)*d)   r(sqrt(2)*d)
                        r(d)           r(sqrt(2)*d)   r(sqrt(5)*d)   1              r(d)
                        r(sqrt(2)*d)   r(d)           r(sqrt(2)*d)   r(d)           1
                        ];
                case 8
                    R = [
                        1              r(d)           r(2*d)         r(d)           r(d)          r(sqrt(5)*d)   r(d)           r(3*d)         r(sqrt(2)*d)
                        r(d)           1              r(d)           r(sqrt(2)*d)   r(sqrt(2)*d)  r(sqrt(2)*d)   r(2*d)         r(2*d)         r(d)
                        r(2*d)         r(d)           1              r(sqrt(5)*d)   r(sqrt(5)*d)  r(d)           r(3*d)         r(d)           r(sqrt(2)*d)
                        r(d)           r(sqrt(2)*d)   r(sqrt(5)*d)   1              r(2*d)        r(sqrt(8)*d)   r(sqrt(2)*d)   r(sqrt(10)*d)  r(d)
                        r(d)           r(sqrt(2)*d)   r(sqrt(5)*d)   r(2*d)         1             r(2*d)         r(sqrt(2)*d)   r(sqrt(10)*d)  r(sqrt(5)*d)
                        r(sqrt(5)*d)   r(sqrt(2)*d)   r(d)           r(sqrt(8)*d)   r(2*d)        1              r(sqrt(10)*d)  r(sqrt(2)*d)   r(sqrt(5)*d)
                        r(d)           r(2*d)         r(3*d)         r(sqrt(2)*d)   r(sqrt(2)*d)  r(sqrt(10)*d)  1              r(4*d)         r(sqrt(5)*d)
                        r(3*d)         r(2*d)         r(d)           r(sqrt(10)*d)  r(sqrt(10)*d) r(sqrt(2)*d)   r(4*d)         1              r(sqrt(5)*d)
                        r(sqrt(2)*d)   r(d)           r(sqrt(2)*d)   r(d)           r(sqrt(5)*d)  r(sqrt(5)*d)   r(sqrt(5)*d)   r(sqrt(5)*d)   1
                        ];
                otherwise
                    error('Only 4 and 8 values supported');
end
L= chol(R,'lower');%利用Cholesky分解得到一个下三角矩阵L满足矩阵方程L*L'=R.
            lambda_n_T  = L(end,:);%取出L矩阵的最后一行
            R_tilde     = R(1:end-1,1:end-1);%将R矩阵降一维保存成为R_tilde
            L_tilde     = chol(R_tilde,'lower');%利用Cholesky分解得到一个下三角矩阵L_tilde满足矩阵方程L_tilde*L_tilde'=R_tilde.
            inv_L_tilde = inv(L_tilde);%求L_tilde矩阵的逆inv_L_tilde

            N_x = size(a_n_matrix,2);
            N_y = size(a_n_matrix,1);
            s = zeros(N_y,N_x,num_eNodeBs);
            for y_=1:N_y
                for x_=1:N_x
                    % Substitutes the LTE_aux_shadowFadingMapClaussen_get_neighbors function
                    % LTE_aux_shadowFadingMapClaussen_get_neighbors(s,position,values)
                    s_tilde = zeros(n_neighbors,num_eNodeBs);
                    positions        = [x_+offsets(:,1) y_+offsets(:,2)];% 计算出理论图样的坐标
                    positions_geq0   = positions>0;% 选出不超过ROI下界的网格图样(0,0)
                    positions_leqroi = [positions(:,1)<=N_x positions(:,2)<=N_y];% 选出不超过ROI上界的网格图样(N_x,N_y)
                    positions_valid  = (positions_geq0(:,1)&positions_geq0(:,2)) & (positions_leqroi(:,1)&positions_leqroi(:,2));% 取以上两个集合的交集
                    for i_=1:n_neighbors
                        if positions_valid(i_)
                            neighbor_position = positions(i_,:);% 选出合法的相邻网格存储到neighbor_position当中
                            s_tilde(i_,:) = s(neighbor_position(2),neighbor_position(1),:);% 与s = zeros(N_y,N_x,num_eNodeBs)相对应
                        end
                    end
                    % end of LTE_aux_shadowFadingMapClaussen_get_neighbors
                    
                    % 参照论文的公式（13）进行思考
                    inv_L_tilde_s_tilde = inv_L_tilde * s_tilde;
                    a_n_all = reshape(a_n_matrix(y_,x_,:),1,[]);
                    inv_L_tilde_s_tilde_a_n = [inv_L_tilde_s_tilde;a_n_all];
                    s(y_,x_,:) = lambda_n_T*inv_L_tilde_s_tilde_a_n;
                end
            end













