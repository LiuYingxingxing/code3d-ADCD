clear;
clc;
h=2060; w=4000;C=[];D=[];E=[];F=[];cm=0;mm_num=0;
bs=load([pwd '\小基线对组合列表_init_37.gdh']);
nimg= size(bs,1); %小基线对个数
bs=double(bs);%数值-字符串
A=bs(:,1);
B=bs(:,2);
workpath=pwd;
bs_name=num2str(bs); 
%先把uintp复制为cp
 headpath1=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_uintP.gdh',nimg,1)];
 headpath11=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_cp.gdh',nimg,1)];
 headpath2=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_uintP.gd',nimg,1)];
 headpath22=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_cp.gd',nimg,1)];
 headpath3=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_uintP.hdr',nimg,1)];
 headpath33=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_cp.hdr',nimg,1)];
 for hp=1:size(bs,1) 
     oldname1=headpath1(hp,:);
     newname1=headpath11(hp,:);
     oldname2=headpath2(hp,:);
     newname2=headpath22(hp,:);
     oldname3=headpath3(hp,:);
     newname3=headpath33(hp,:);    
     copyfile(oldname1,newname1)
     copyfile(oldname2,newname2)
     copyfile(oldname3,newname3)
 end
 phasepath=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_cp.gd',nimg,1)];  
 wrap_phasepath=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_diffF.gd',nimg,1)];
 name=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_cp.gd',nimg,1)];
%找到闭合环
parfor i=1:nimg
    id1=find(A==B(i));
   for j=1:size(id1)
       id2=find(B==B(id1(j)));
    for k=1:size(id2)
       if(A(i)==A(id2(k)))
        cm=cm+1;
        f=[i,id1(j),id2(k)];
        C=[C;i];
        D=[D;id1(j)];
        E=[E;id2(k)];
        F=[F;f];
         break;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
       end
    end
   end
end
%相位梯度图
diff_total=zeros(h,w,nimg);
phase_total=zeros(h,w,nimg);
for y=1:nimg
    phase= multibandread(phasepath(y,:),[h,w,1],'float',0,'bsq','l');
    phase_total(:,:,y)=phase;
    %列坐标相位差
%     col=zeros(h,w-1);
    phase_mask=phase;
    phase_mask(phase==0)=NaN;
    diff_col=diff(phase_mask,1,2);
%     col(abs(diff_col)>=4.7)=1;
    diff_col(isnan(diff_col)) = 0;
    col=round(diff_col/(2*pi));
    col_int=[zeros(size(col,1),1) col];
    col_trans_init=[col zeros(size(col,1),1)];
    %行坐标相位差
%     row=zeros(h-1,w);
    diff_row=diff(phase_mask,1,1);
    diff_row(isnan(diff_row)) = 0;
%     row(abs(diff_row)>=4.7)=1;
    row=round(diff_row/(2*pi));
    row_int=[zeros(1,size(row,2));row];
    row_trans_init=[row;zeros(1,size(row,2))];
    diff_phase=row_int+col_int+col_trans_init+row_trans_init;
    diff_total(:,:,y)=diff_phase;
end
%计算闭合相位并对闭合相位分块及探测每个块状的干涉图
Storage_interferogram_Error=zeros(h,w,nimg);
close_block=zeros(h,w,cm);
for x=1:cm
     closenew=zeros(h,w);
     block_storage=zeros(h,w,3);
     A_diff=diff_total(:,:,C(x,:));
     B_diff=diff_total(:,:,D(x,:));
     C_diff=diff_total(:,:,E(x,:));
     [row_Adiff,col_Adiff]=find(A_diff~=0);
     nozeros_Adiff=[row_Adiff col_Adiff];
     [row_Bdiff,col_Bdiff]=find(B_diff~=0);
     nozeros_Bdiff=[row_Bdiff col_Bdiff];  
     [row_Cdiff,col_Cdiff]=find(C_diff~=0);
     nozeros_Cdiff=[row_Cdiff col_Cdiff];     
     phase1 = multibandread(phasepath(F(x,1),:),[h,w,1],'float',0,'bsq','l');
     phase2 = multibandread(phasepath(F(x,2),:),[h,w,1],'float',0,'bsq','l');
     phase3 = multibandread(phasepath(F(x,3),:),[h,w,1],'float',0,'bsq','l');
     wrap_phase1 = multibandread(wrap_phasepath(F(x,1),:),[h,w,2],'float',0,'bip','l');
     wrap_phase2 = multibandread(wrap_phasepath(F(x,2),:),[h,w,2],'float',0,'bip','l');
     wrap_phase3 = multibandread(wrap_phasepath(F(x,3),:),[h,w,2],'float',0,'bip','l');
     wrap_phase1 =(complex(wrap_phase1(:,:,1),wrap_phase1(:,:,2)));
     wrap_phase2 =(complex(wrap_phase2(:,:,1),wrap_phase2(:,:,2)));
     wrap_phase3 =(complex(wrap_phase3(:,:,1),wrap_phase3(:,:,2)));
     close_phase=phase1+phase2-phase3; 
     close_wrap=angle(wrap_phase1.*wrap_phase2 .*conj(wrap_phase3));
     close=round((close_phase-close_wrap)/(2*pi));  
     mask=phase1==0|phase2==0|phase3==0;
     %与掩膜边界相邻的点，换句话就是找到块边界的有效点
     se = strel('arbitrary', [0 1 0; 1 1 1; 0 1 0]);
     dilatedMask = imdilate(mask, se);
     adjacentPixelsMask = dilatedMask & ~mask; 
     [row_adjacent,col_adjacent]=find(adjacentPixelsMask~=0);
     adjacent_boundary_point=[row_adjacent col_adjacent];
     %解缠相对坐标问题，让解缠错误区域，为小区域
     close_zero=size(find(close==0))-size(find(mask==1));
     close_fzero=size(find(close~=0));
     if close_fzero(1,1)>close_zero(1,1)
        close_A=(close(:))';
        close_A(close_A(:)==0)=[];
        close=close-ones(h,w)*mode(close_A);
     end  
    close(mask==1)=0;  
    %不同close值
    nonzero_nums=close(close~=0);
    unique_nonzero_nums=unique(nonzero_nums);
    %分块
    for q=1:numel(unique_nonzero_nums)
        error_values=(close==unique_nonzero_nums(q,1));
        SZ=size(error_values);%获取图像的大小
        [bwlabel_all,bwlabel_num]=bwlabel(error_values);
        bwlabel_area=regionprops(bwlabel_all,'Area');
        Area=[bwlabel_area.Area];
        %解缠错误块的大小
        Label=find(Area>800);
        bw=ismember(bwlabel_all,Label);
        close_error=close;
        close_error((bw==0))=0;
        close_error(closenew~=0)=0;
        [bwlabel_new,bwlabelnew_num]=bwlabel(close_error);
        %块状探测
        for mm=1:bwlabelnew_num 
            mm_num=mm_num+1;
            error=close_error;
            error(bwlabel_new~=mm)=0;
            boundary=bwboundaries(error,'noholes');
            boundary=cell2mat(boundary);
            boundary_display=zeros(SZ);
            for vv=1:length(boundary)%可以替换成boundaryCW
                boundary_display(sub2ind(SZ,boundary(:,1),boundary(:,2)))=1;
            end
            %这个地方其实和error的东西是一样的
            boundary_into=poly2mask(boundary(:,2),boundary(:,1),SZ(1),SZ(2));   
            result=(boundary_into|boundary_display);
            result(close_block(:,:,x)~=0)=0;
            close_block(:,:,x)=close_block(:,:,x)+result*mm_num; 
            %总边界与掩膜边界相邻的点
            [is_adjacentPixelsMask, idx_adjacentPixelsMask] = ismember(boundary, adjacent_boundary_point, 'rows');
            valid_point=boundary(~is_adjacentPixelsMask, :);
            out_M=(size(boundary,1)-size(valid_point,1))/size(boundary,1);  
            %阈值M
            threshold_M=0.4;
            %边界与相位梯度非零值的匹配情况
            %干涉图A的匹配情况
            A_nonum=intersect(nozeros_Adiff,valid_point,'rows');
            A_vaidnum=size(A_nonum,1)/size(valid_point,1);%有效值情况
            %干涉图B的匹配情况
            B_nonum=intersect(nozeros_Bdiff,valid_point,'rows');
            B_vaidnum=size(B_nonum,1)/size(valid_point,1);%有效值情况
            %干涉图C的匹配情况
            C_nonum=intersect(nozeros_Cdiff,valid_point,'rows');
            C_vaidnum=size(C_nonum,1)/size(valid_point,1);%有效值情况  
            %改正值准备
             [N,~,ix]=unique(boundary(:,1));
             [M,~,iy]=unique(boundary(:,2));
             G=[];H=[];GG=[];HH=[];G_gd_A=[];G_gd_B=[];G_gd_C=[];H_gd_A=[];H_gd_B=[];H_gd_C=[];
             for n=1:length(N)
                 rows=find(ix==n);
                 max_val=max(boundary(rows,2));
                 min_val=min(boundary(rows,2));
                 %每个区域的边界最大值和最小值
                 if(mask(boundary(rows(1,1),1),min_val-1)==0 && close(boundary(rows(1,1),1),min_val-1)==0) 
                     G=[G;boundary(rows(1,1),1),min_val];
                     GG=[GG;boundary(rows(1,1),1),min_val-1];
                 end
                 if(mask(boundary(rows(1,1),1),max_val+1)==0 && close(boundary(rows(1,1),1),max_val+1)==0)
                     H=[H;boundary(rows(1,1),1),max_val];
                     HH=[HH;boundary(rows(1,1),1),max_val+1];
                 end          
             end   
            %最大的且大于阈值
            if A_vaidnum>=B_vaidnum && A_vaidnum>=C_vaidnum && A_vaidnum>= threshold_M
             %改正值A
             for g=1:size(G,1)
                 G_grent_A=round((phase1(G(g,1),G(g,2))-phase1(GG(g,1),GG(g,2)))/(2*pi));
                 G_gd_A=[G_gd_A;G_grent_A];
             end
             for hh=1:size(H,1)
                 H_grent_A=round((phase1(H(hh,1),H(hh,2))-phase1(HH(hh,1),HH(hh,2)))/(2*pi));
                 H_gd_A=[H_gd_A;H_grent_A];
             end 
             gd_A=[G_gd_A;H_gd_A];
             gd_A=gd_A(gd_A~=0);
             gd_mode_A=mode(gd_A);
             if isnan(gd_mode_A)
                gd_mode_A=0;
             end
             result(block_storage(:,:,1)~=0)=0;
             block_storage(:,:,1)=block_storage(:,:,1)+result*gd_mode_A;          
            elseif B_vaidnum>=A_vaidnum && B_vaidnum>=C_vaidnum && B_vaidnum>= threshold_M
             %改正值B
             for g=1:size(G,1)
                 G_grent_B=round((phase2(G(g,1),G(g,2))-phase2(GG(g,1),GG(g,2)))/(2*pi));
                 G_gd_B=[G_gd_B;G_grent_B];
             end
             for hh=1:size(H,1)
                 H_grent_B=round((phase2(H(hh,1),H(hh,2))-phase2(HH(hh,1),HH(hh,2)))/(2*pi));
                 H_gd_B=[H_gd_B;H_grent_B];
             end 
             gd_B=[G_gd_B;H_gd_B];
             gd_B=gd_B(gd_B~=0);
             gd_mode_B=mode(gd_B);
             if isnan(gd_mode_B)
                gd_mode_B=0;
             end
             result(block_storage(:,:,2)~=0)=0;
             block_storage(:,:,2)=block_storage(:,:,2)+result*gd_mode_B;
            elseif C_vaidnum>=A_vaidnum && C_vaidnum>=B_vaidnum && C_vaidnum>= threshold_M
             %改正值C
             for g=1:size(G,1)
                 G_grent_C=round((phase3(G(g,1),G(g,2))-phase3(GG(g,1),GG(g,2)))/(2*pi));
                 G_gd_C=[G_gd_C;G_grent_C];
                 G_gd_C=[];
             end
             for hh=1:size(H,1)
                 H_grent_C=round((phase3(H(hh,1),H(hh,2))-phase3(HH(hh,1),HH(hh,2)))/(2*pi));
                 H_gd_C=[H_gd_C;H_grent_C];
             end 
             gd_C=[G_gd_C;H_gd_C];
             gd_C=gd_C(gd_C~=0);     
             gd_mode_C=mode(gd_C); 
             if isnan(gd_mode_C)
                gd_mode_C=0;
             end
             result(block_storage(:,:,3)~=0)=0;
             block_storage(:,:,3)=block_storage(:,:,3)+result*gd_mode_C;    
            end
        end
    end
    phase_a=phase1-2*pi*block_storage(:,:,1);
    phase_b=phase2-2*pi*block_storage(:,:,2);
    phase_c=phase3-2*pi*block_storage(:,:,3);
    multibandwrite(phase_a,name(C(x,:),:),'bsq','precision','float');
    multibandwrite(phase_b,name(D(x,:),:),'bsq','precision','float');
    multibandwrite(phase_c,name(E(x,:),:),'bsq','precision','float');
end  


        
        
        
        
        
        