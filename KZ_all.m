% % ����kz�����õ���ز���
% % kz_para = struct('base_perp',[],'c',[],'p',[],'R',[],'incidence',[]);
% 
% % %% 0413-0424
% % kz_para.base_perp = 79.4813; %��ֱ���ߣ�ʹ��gamma�����base_calc������м���
% % kz_para.c = 299792458; %΢���ٶ�
% % kz_para.p = 9.65*10^9; %΢��Ƶ�ʣ�������Ӱ��par�ļ��е�radar_frequency��һ��
% % kz_para.R = 593428.2010; %б�࣬������Ӱ��par�ļ��е�center_range_slc��һ��
% % kz_para.incidence = 31.5839; %����ǣ�������Ӱ��par�ļ��е�incidence_angle��һ��
% % save('0413_0424_kz_para','kz_para');
% % 
% % %% 0505-0424
% % kz_para.base_perp = 255.6263; %��ֱ���ߣ�ʹ��gamma�����base_calc������м���
% % kz_para.c = 299792458; %΢���ٶ�
% % kz_para.p = 9.65*10^9; %΢��Ƶ�ʣ�������Ӱ��par�ļ��е�radar_frequency��һ��
% % kz_para.R = 593428.2010; %б�࣬������Ӱ��par�ļ��е�center_range_slc��һ��
% % kz_para.incidence = 31.5839; %����ǣ�������Ӱ��par�ļ��е�incidence_angle��һ��
% % save('0505_0424_kz_para','kz_para');
% % 
% % %% puer 20150206TSX_TDX
% % kz_para.base_perp = 618.70010; %��ֱ���ߣ�ʹ��gamma�����base_calc������м���
% % kz_para.c = 299792458; %΢���ٶ�
% % kz_para.p = 9.6499993e+09; %΢��Ƶ�ʣ�������Ӱ��par�ļ��е�radar_frequency��һ��
% % kz_para.R = 646029.1959; %б�࣬������Ӱ��par�ļ��е�center_range_slc��һ��
% % kz_para.incidence = 39.7289; %����ǣ�������Ӱ��par�ļ��е�incidence_angle��һ��
% % save('puer_20150206TSX_TDX_kz_para','kz_para');
% 
% 
% 
% 
% KZ_PARA = zeros(91,9);
% % ��ʼ��kz����������б�����91�����߶ԣ�ÿ�Զ���8������
% % �ֱ��ǣ���Ӱ�����ơ���Ӱ�����ơ���ֱ���߳��ȡ����١�����Ƶ�ʡ�б�ࡢ����Ǻ���ǰ����������kz��ģ���ߡ�
% 
% TEMP = importdata('J:\data\TSX_TDX_yunnan_puer\bperp_file.xlsx');
% KZ_PARA(:,1:3) = TEMP.data(:,2:4);% ����Ӱ�����ơ���Ӱ�����ƺʹ�ֱ���߳��ȴ�berpfile�ж���ȥ��
% KZ_PARA(:,3) = abs(KZ_PARA(:,3));
% 
% KZ_PARA(:,4) = 299792458 * ones(91,1);% ��ȡ����
% 
% KZ_PARA(:,5) = 9.6499993e+09 * ones(91,1);% ��ȡ����Ƶ��
% 
% TEMP1 = importdata('J:\data\TSX_TDX_yunnan_puer\master_para.xlsx');% ��ȡ��Ӱ����ز���excel
% 
% for i = 1:91
%     [itemp,jtemp] = find(KZ_PARA(i,1) == TEMP1.data(:,1));
%     KZ_PARA(i,6) = TEMP1.data(itemp,3);% ��ȡб��
%     KZ_PARA(i,7) = TEMP1.data(itemp,4);% ��ȡ�����
%     % ����kz
%     KZ_PARA(i,8) = (4 * pi * KZ_PARA(i,3) * KZ_PARA(i,5)) / (KZ_PARA(i,4) *KZ_PARA(i,6) * sind(KZ_PARA(i,7)));
%     % ����ģ����
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

% ����,��Ӱ������,��ֱ���ߣ����٣�����Ƶ�ʴ��ȥ
KZ(:,1:3) = temp.textdata.Sheet2;
KZ{1,4} = 'C';
KZ{1,5} = 'radar_frequency';
for i = 2:32
    KZ{i,3} = temp.data.Sheet2(i-1  );
    KZ{i,4} = 299792458;% ����
    KZ{i,5} = 9.6499993e+09; %������Ƶ��
end

temp2 = importdata('master_para.xlsx');

% ��ȡб��������
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
    
% ����kz
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

















    



