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
                case 8 % offsets��ע�����ᵽ�������йأ������ֵ���������Figure5�Ĳ���s_nͼ�����λ�õ�����ƫ���һ��Ϊ�����꣬�ڶ���Ϊ�����ꡣ
                       % ���磬offsets�ĵ�һ��-1 -1����s_n�����ƶ�һ�������ƶ�һ����õ�s_1
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
            r = @(d) exp(-alpha*d);% �����еĹ�ʽ��5��
            d = resolution; % How many meters a hop is %Ҳ������ν������ֱ���

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
                    a_n_matrix(:,:,i_) = sqrt(r_eNodeBs)*a_n_original_map + sqrt(1-r_eNodeBs)*a_n_random_matrix(:,:,i_);% ������
                end
switch n_neighbors %�μ����׵�conclusions����
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
L= chol(R,'lower');%����Cholesky�ֽ�õ�һ�������Ǿ���L������󷽳�L*L'=R.
            lambda_n_T  = L(end,:);%ȡ��L��������һ��
            R_tilde     = R(1:end-1,1:end-1);%��R����һά�����ΪR_tilde
            L_tilde     = chol(R_tilde,'lower');%����Cholesky�ֽ�õ�һ�������Ǿ���L_tilde������󷽳�L_tilde*L_tilde'=R_tilde.
            inv_L_tilde = inv(L_tilde);%��L_tilde�������inv_L_tilde

            N_x = size(a_n_matrix,2);
            N_y = size(a_n_matrix,1);
            s = zeros(N_y,N_x,num_eNodeBs);
            for y_=1:N_y
                for x_=1:N_x
                    % Substitutes the LTE_aux_shadowFadingMapClaussen_get_neighbors function
                    % LTE_aux_shadowFadingMapClaussen_get_neighbors(s,position,values)
                    s_tilde = zeros(n_neighbors,num_eNodeBs);
                    positions        = [x_+offsets(:,1) y_+offsets(:,2)];% ���������ͼ��������
                    positions_geq0   = positions>0;% ѡ��������ROI�½������ͼ��(0,0)
                    positions_leqroi = [positions(:,1)<=N_x positions(:,2)<=N_y];% ѡ��������ROI�Ͻ������ͼ��(N_x,N_y)
                    positions_valid  = (positions_geq0(:,1)&positions_geq0(:,2)) & (positions_leqroi(:,1)&positions_leqroi(:,2));% ȡ�����������ϵĽ���
                    for i_=1:n_neighbors
                        if positions_valid(i_)
                            neighbor_position = positions(i_,:);% ѡ���Ϸ�����������洢��neighbor_position����
                            s_tilde(i_,:) = s(neighbor_position(2),neighbor_position(1),:);% ��s = zeros(N_y,N_x,num_eNodeBs)���Ӧ
                        end
                    end
                    % end of LTE_aux_shadowFadingMapClaussen_get_neighbors
                    
                    % �������ĵĹ�ʽ��13������˼��
                    inv_L_tilde_s_tilde = inv_L_tilde * s_tilde;
                    a_n_all = reshape(a_n_matrix(y_,x_,:),1,[]);
                    inv_L_tilde_s_tilde_a_n = [inv_L_tilde_s_tilde;a_n_all];
                    s(y_,x_,:) = lambda_n_T*inv_L_tilde_s_tilde_a_n;
                end
            end













