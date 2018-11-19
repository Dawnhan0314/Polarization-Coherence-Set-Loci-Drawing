

%% 读取T6

T6_path = 'J:\GAOHAN\EXPERIMENT\8.NR_Puer\slc_sub\Experiment_data\long_Bperp\20141213_TSX_20141224_TDX\20141213_TSX_20141224_TDX_FER_LEE_7_7\T6';
[T6,row,col] = readT6PUER(T6_path);

%%%%%%%%%%%%%%%% 单点NR %%%%%%%%%%%%%%%%%

%% 单点NR
AW1 = freadbk([T6_path,'\cmplx_coh_Opt_NR1.bin'],row,'cpxfloat32');
AW2 = freadbk([T6_path,'\cmplx_coh_Opt_NR2.bin'],row,'cpxfloat32');
AW3 = freadbk([T6_path,'\cmplx_coh_Opt_NR3.bin'],row,'cpxfloat32');
HH = freadbk([T6_path,'\cmplx_coh_HH.bin'],row,'cpxfloat32');
VV = freadbk([T6_path,'\cmplx_coh_VV.bin'],row,'cpxfloat32');
HV = freadbk([T6_path,'\cmplx_coh_HV.bin'],row,'cpxfloat32');

spauli = imread('J:\GAOHAN\EXPERIMENT\8.NR_Puer\slc_sub\Experiment_data\long_Bperp\20141213_TSX_20141224_TDX\20141213_TSX_20141224_TDX_FER_LEE_7_7\T6\cmplx_coh_Opt_NR1_mod.bmp');

[ROIrow,ROIcol] = NR_ROI_Point_Building(spauli);% 选择感兴趣点

% [SelectAW1,SelectAW2,SelectAW3] = NR_extraction(0,pi,SelectT6);
SelectAW1 = AW1(ROIrow,ROIcol);
SelectAW2 = AW2(ROIrow,ROIcol);
SelectAW3 = AW3(ROIrow,ROIcol);

figure;polar(0,1);hold on;
[CoPoint4] = CoherenceSetLoci(T6,ROIrow,ROIcol,'Point',0);
polar(angle(CoPoint4),abs(CoPoint4),'*');hold on;
NR1 = polar(angle(SelectAW1),abs(SelectAW1),'ro');hold on;
NR2 = polar(angle(SelectAW2),abs(SelectAW2),'r*');hold on;
NR3 = polar(angle(SelectAW3),abs(SelectAW3),'r^');hold on;




%%%%%%%%%%%%%%%% 窗口内平均方法 %%%%%%%%%%%%%%%%%
% 窗口内平均
MEAN = cell(2,3);

%% 勾画要画相干集区域的像素

AW1 = freadbk([T6_path,'\cmplx_coh_Opt_NR1.bin'],row,'cpxfloat32');
AW2 = freadbk([T6_path,'\cmplx_coh_Opt_NR2.bin'],row,'cpxfloat32');
AW3 = freadbk([T6_path,'\cmplx_coh_Opt_NR3.bin'],row,'cpxfloat32');
HH = freadbk([T6_path,'\cmplx_coh_HH.bin'],row,'cpxfloat32');
VV = freadbk([T6_path,'\cmplx_coh_VV.bin'],row,'cpxfloat32');
HV = freadbk([T6_path,'\cmplx_coh_HV.bin'],row,'cpxfloat32');
dif = angle(AW1) - angle(AW2);
kz = 0.8769;
% spauli = imread('F:\PUER\slc_sub\puer_930_760_T3NL\C3\SinclairRGB.bmp');
spauli = imread('F:\PUER\slc_sub\Experiment_data\long_Bperp\20141213_TDX_20141224_TSX\20141213_TDX_20141224_TSX_FER_LEE_7_7\T6\PauliRGB_T1.bmp');
ROIregion = NR_ROI_Building(spauli);% 勾画叠掩区域
%% 读取对应极化通道下的相干
j =1;
meanT6 = zeros(6,6);
for i = 1:36
    meanT6(i) = mean(T6{i}(ROIregion==1));
end
[meanAW1,meanAW2,meanAW3] = NR_extraction(0,pi*0.5,meanT6);
MEAN{j,1} = meanT6;MEAN{j,2} = meanAW1;MEAN{j,3} = meanAW2;MEAN{j,4}=meanAW3;

figure;polar(0,1);hold on;


[CoPoint4] = CoherenceSetLoci(MEAN{j,1},0,0,'NR',0);
polar(angle(CoPoint4),abs(CoPoint4));hold on;
NR1 = polar(angle(MEAN{j,2}),abs(MEAN{j,2}),'ro');hold on;
NR2 = polar(angle(MEAN{j,3}),abs(MEAN{j,3}),'r*');hold on;
NR3 = polar(angle(MEAN{j,4}),abs(MEAN{j,4}),'r^');hold on;

phaseHH = angle(HH);
meanHH = mean(phaseHH(ROIregion==1 ));


