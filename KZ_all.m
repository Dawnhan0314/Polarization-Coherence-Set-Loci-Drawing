% % 计算kz计算用的相关参数
% % kz_para = struct('base_perp',[],'c',[],'p',[],'R',[],'incidence',[]);
% 
% % %% 0413-0424
% % kz_para.base_perp = 79.4813; %垂直基线，使用gamma软件的base_calc命令进行计算
% % kz_para.c = 299792458; %微波速度
% % kz_para.p = 9.65*10^9; %微波频率，参照主影像par文件中的radar_frequency这一项
% % kz_para.R = 593428.2010; %斜距，参照主影像par文件中的center_range_slc这一项
% % kz_para.incidence = 31.5839; %入射角，参照主影像par文件中的incidence_angle这一项
% % save('0413_0424_kz_para','kz_para');
% % 
% % %% 0505-0424
% % kz_para.base_perp = 255.6263; %垂直基线，使用gamma软件的base_calc命令进行计算
% % kz_para.c = 299792458; %微波速度
% % kz_para.p = 9.65*10^9; %微波频率，参照主影像par文件中的radar_frequency这一项
% % kz_para.R = 593428.2010; %斜距，参照主影像par文件中的center_range_slc这一项
% % kz_para.incidence = 31.5839; %入射角，参照主影像par文件中的incidence_angle这一项
% % save('0505_0424_kz_para','kz_para');
% % 
% % %% puer 20150206TSX_TDX
% % kz_para.base_perp = 618.70010; %垂直基线，使用gamma软件的base_calc命令进行计算
% % kz_para.c = 299792458; %微波速度
% % kz_para.p = 9.6499993e+09; %微波频率，参照主影像par文件中的radar_frequency这一项
% % kz_para.R = 646029.1959; %斜距，参照主影像par文件中的center_range_slc这一项
% % kz_para.incidence = 39.7289; %入射角，参照主影像par文件中的incidence_angle这一项
% % save('puer_20150206TSX_TDX_kz_para','kz_para');
% 
% 
% 
% 
% KZ_PARA = zeros(91,9);
% % 初始化kz计算参数的列表，共有91个基线对，每对都有8个参数
% % 分别是：主影像名称、从影像名称、垂直基线长度、光速、波段频率、斜距、入射角和由前面参数计算的kz和模糊高。
% 
% TEMP = importdata('J:\data\TSX_TDX_yunnan_puer\bperp_file.xlsx');
% KZ_PARA(:,1:3) = TEMP.data(:,2:4);% 将主影像名称、从影像名称和垂直基线长度从berpfile中读进去。
% KZ_PARA(:,3) = abs(KZ_PARA(:,3));
% 
% KZ_PARA(:,4) = 299792458 * ones(91,1);% 读取光速
% 
% KZ_PARA(:,5) = 9.6499993e+09 * ones(91,1);% 读取波段频率
% 
% TEMP1 = importdata('J:\data\TSX_TDX_yunnan_puer\master_para.xlsx');% 读取主影像相关参数excel
% 
% for i = 1:91
%     [itemp,jtemp] = find(KZ_PARA(i,1) == TEMP1.data(:,1));
%     KZ_PARA(i,6) = TEMP1.data(itemp,3);% 读取斜距
%     KZ_PARA(i,7) = TEMP1.data(itemp,4);% 读取入射角
%     % 计算kz
%     KZ_PARA(i,8) = (4 * pi * KZ_PARA(i,3) * KZ_PARA(i,5)) / (KZ_PARA(i,4) *KZ_PARA(i,6) * sind(KZ_PARA(i,7)));
%     % 计算模糊高
%     KZ_PARA(i,9) = 2 * pi / KZ_PARA(i,8);
% end
% % 
% % kzparatemp = KZ_PARA;
% % for i = 1:91
% %     if(KZ_PARA(i,1) == KZ_PARA(i,2))
% %         continue;
% %     else
% %         kzparatemp(i,:) = nan;
% %     end
% % end

%% kz_all
clc;
clear;

cd 'F:\PUER\slc_sub';

KZ = cell(32,9);

temp = importdata('bperp_file.xlsx');

% 将主,从影像名字,垂直基线，光速，波段频率存进去
KZ(:,1:3) = temp.textdata.Sheet2;
KZ{1,4} = 'C';
KZ{1,5} = 'radar_frequency';
for i = 2:32
    KZ{i,3} = temp.data.Sheet2(i-1  );
    KZ{i,4} = 299792458;% 光速
    KZ{i,5} = 9.6499993e+09; %　波段频率
end

temp2 = importdata('master_para.xlsx');

% 读取斜距和入射角
KZ(1,6:7) = temp2.textdata(1,3:4);

for i = 2:32
for j = 1:7
    if(num2str(temp2.data(j,1)) == KZ{i,1}(1:8))
        KZ{i,6} = temp2.data(j,3);
        KZ{i,7} = temp2.data(j,4);
    else
        continue;
    end
end
end
    
% 计算kz
KZ{1,8} = 'KZ';
KZ{1,9} = 'amHeight';
for i = 2:32
    KZ{i,8} = (4 * pi * KZ{i,3} * KZ{i,5}) / (KZ{i,4} *KZ{i,6} * sind(KZ{i,7}));
    KZ{i,9} = (2 * pi) / KZ{i,8};
end

save('KZ.mat','KZ');


bperp = 230.72480;
f = 9.6499993e+09;
c = 299792458;
R = 635060.6925;
incidence = 37.8036;
kz = (4 * pi * bperp * f) / (c * R * sind(incidence));

















    



