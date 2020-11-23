function ACOfirst()
clear;
clc;
%author Xing Peng
%@copyright reserved
%% 测试矩阵D
D=[0 1 3 0;   1 0 1 2;   3 1 0 5;   0 2 5 0;  ];
%% get 邻接矩阵
%D=getD();
%% 将不连通赋值为0
N=size(D,1);
for i=1:N
    for j=1:N
        if D(i,j)==inf
            D(i,j)=0;
        end
    end
end
%%
%% 起点与终点
S=1 ;         % S 起始点（最短路径的起始点）
E=N;         % E 终止点（最短路径的目的点）
%% 
Tau=ones(N,N);% Tau 初始信息素矩阵（认为前面的觅食活动中有残留的信息素）
Tau=8.*Tau;
K=10;          % K 迭代次数（指蚂蚁出动多少波）
M=5;         % M 蚂蚁个数（每一波蚂蚁有多少个）
Alpha=2;    % Alpha 表征信息素重要程度的参数
Beta=5;     % Beta 表征启发式因子重要程度的参数
Rho=0.1 ;     % Rho 信息素蒸发系数
Q=1.2;          % Q 信息素增加强度系数
minkl=inf;
mink=0;
minl=0;
dlmwrite('distance.txt',D,'delimiter',' ') 
%xlswrite('d',D,sheet,range);%输出excel表
%% 启发式信息，取为当前点至最终目标终点的直线距离的倒数
Eta=zeros(N,1);
for i=1:N
   if i~=E
       if D(i,N)~=0 
           Eta(i)=1/D(i,N);
       end
   else
       Eta(i)=10;
   end
end
ROUTES=cell(K,M);%用细胞结构存储每一代的每一只蚂蚁的爬行路线
PL=zeros(K,M);%用矩阵存储每一代的每一只蚂蚁的爬行路线长度
%% -----------启动K轮蚂蚁觅食活动，每轮派出M只蚂蚁--------------------
for k=1:K
    for m=1:M
        %% 第一步：状态初始化
        W=S;%当前节点初始化为起始点，W为蚂蚁所在节点，
        Path=S;%爬行路线初始化
        PLkm=0;%爬行路线长度初始化
        TABUkm=ones(N,1);%禁忌表初始化，所有都为1，表示没有走过
        TABUkm(S,1)=0;%已经在初始点了，因此要排除
        DD=D;%邻接矩阵初始化
        %% 第二步：下一步可以前往的节点
        DW=DD(W,:); % 选定第W个点，W为地图编码，把邻接矩阵中第W行给DW
                DW1=find(DW);% 找出DW这一行中，不为零的数，其所在的列序列，比如，W=1，则DW1=2 4，就是说第2,4与1有联系
        for j=1:length(DW1)    
            if TABUkm(DW1(j),1)==0 % DW1(1)=2, TABUkm(2)=1, TABUkm(2)表示第2行的第1个数
                DW(DW1(j))=0;
            end
        end        
        LJD=find(DW);% 2,4
        Len_LJD=length(LJD);%可选节点的个数,2
        iter=0;
        %% 觅食停止条件：蚂蚁未遇到食物或者陷入死胡同，第一代第一只蚂蚁
        while W~=E&&Len_LJD>=1 %终止条件为，W为终点，或者，可选节点个数大于等于1
            iter=iter+1;
            %% 第三步：转轮赌法选择下一步怎么走
            PP=zeros(Len_LJD,1); % 根据蚂蚁当前所在节点的下一个可选节点数，建立0矩阵行数为可选节点数，列数为1
            for i=1:Len_LJD
                PP(i)=(Tau(W,LJD(i))^Alpha)*((Eta(LJD(i)))^Beta);
            end           
            sumpp=sum(PP); %求PP这一列的和，
            PP=PP/sumpp;%建立概率分布，即这一列中，每个元素代表可选节点被选中的概率            
            Pcum(1)=PP(1);%第一个可行节点被选中的概率
            for i=2:Len_LJD
                Pcum(i)=Pcum(i-1)+PP(i);
            end
            Pcum(1:Len_LJD); %概率累加，最后一个概率是1
            Select=find(Pcum>=rand);           
            to_visit=LJD(Select(1));% 确定下一步要走的节点            
            %% 第四步：状态更新和记录            
            Path=[Path,to_visit];%蚂蚁从W移动到to_visit，节点序列增加，路径长度增加路径节点序列拼接，矩阵拼接
            PLkm=PLkm+DD(W,to_visit);%路径长度增加
            W=to_visit;%蚂蚁移到下一个节点            
            for kk=1:N %搜索每一个节点
                if TABUkm(kk,1)==0 %如果这个节点走过了，为0表示走过了
                    DD(W,kk)=0; % 把与走过点的长度赋值为0
                    DD(kk,W)=0;
                end
            end
            TABUkm(W,1)=0;%已访问过的节点从禁忌表中删除重复再下一步可以前往的节点
            DW=DD(W,:);
            DW1=find(DW);
            for j=1:length(DW1)
                if TABUkm(DW1(j),1)==0
                    DW(j)=0;
                end
            end
            LJD=find(DW);
            Len_LJD=length(LJD);%可选节点的个数            
        end
        %Path,PLkm,%一代一只蚂蚁的信息
        %% 第五步：记下每一代每一只蚂蚁的觅食路线和路线长度
        ROUTES{k,m}=Path;%ROUTES,细胞结构，每一个元素都是一个数组
        if Path(end)==E %如果路径序列的末尾是终点，即该序列有效
            PL(k,m)=PLkm; %把长度信息给第k行，m列，代表第k代第m只蚂蚁的路径长度
            if PLkm<minkl %初始minkl为inf,第一只蚂蚁后，更新
                minkl=PLkm; %把路径长度给minkl
                mink=k; % 代数
                minl=m; % 第几只蚂蚁               
            end
        else
            PL(k,m)=0;
        end
    end    
    %% 第六步：更新信息素
    Delta_Tau=zeros(N,N);%更新量初始化
    for m=1:M
        if PL(k,m)~=0 %如果该路径长度有效
            ROUT=ROUTES{k,m}; %序列给路径
            TS=length(ROUT)-1;%跳数
            PL_km=PL(k,m);%长度给out
            for s=1:TS
                x=ROUT(s);
                y=ROUT(s+1);
                Delta_Tau(x,y)=Delta_Tau(x,y)+Q/PL_km;
                Delta_Tau(y,x)=Delta_Tau(y,x)+Q/PL_km;
            end
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau;%信息素挥发一部分，新增加一部分
end
%% 输出结果
Path
PLkm

%% 绘收敛曲线
minPL=zeros(K);%K是迭代次数，
for i=1:K
     PLK=PL(i,:);%PL是用矩阵存储每一代的每一只蚂蚁的爬行路线长度，把第i次迭代的蚂蚁路线长度给PLK
     Nonzero=find(PLK);%找出PLK中不为0的数的位置
     PLKPLK=PLK(Nonzero);
     minPL(i)=min(PLKPLK);
end
figure(1)
plot(minPL);
hold on
grid on
title('收敛曲线（最小路径长度）');
xlabel('迭代次数');
ylabel('路径长度');



%% ---------------------------绘图--------------------------------

%% 邻接矩阵D的生成
function D=getD()
%% 虚拟城市节点距离邻接矩阵citydistance
% 扩展城市坐标
city=[
    % 汽车站  港口    火车站     机场
    16,88;   9,84;   19,82;    17,98;
    42,98;  38,92;   44,95;    32,97;
    68,96;  65,93;   70,92;    58,98;
    96,78;  93,75;   98,70;    97,87;
    
    8,70;   3,59;   11,61;     9,80;
    36,68;  30,61;   39,62;    36,76;
    69,58;  63,55;   68,50;    68,67;
    
    24,39;  20,37;   26,36;    24,47;
    64,37;  60,34;   66,35;    63,46;
    96,38;  93,35;   98,30;    97,48;
    
    44,22;  40,19;   46,20;    44,30;
    84,14;  80,8;    86,9;     84,23;];
[city_number,city_dim]=size(city);

% 初始城市距离distance矩阵,即每两个节点之间有距离
citydistance=[city_number,city_number];
for i=1:city_number
    for j=1:city_number   
        citydistance(i,j)=((city(i,1)-city(j,1))^2+(city(i,2)-city(j,2))^2)^0.5;
        citydistance(j,i)=citydistance(i,j);                        %对称矩阵
    end
end

% 默认不连通的城市节点,illusion_city矩阵表示虚拟城市所在模型中像素编码，暂定为，(公路，水路，铁路，机场)，城市修改，增加，只需对矩阵进行增改
% 即在此处，某一个城市的公路，默认无法直接连通到另一个城市的铁路或者水路
illusion_city_code=[ 1,2,3,4;
    5,6,7,8;
    9,10,11,12;
    13,14,15,16;
    17,18,19,20;
    21,22,23,24;
    25,26,27,28;
    29,30,31,32;
    33,34,35,36;
    37,38,39,40;
    41,42,43,44;
    45,46,47,48;];
[illusion_city_array,illusion_city_col]=size(illusion_city_code);%取矩阵的大小，行数，列数
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            if any(b_array==m)||any(b_col==m)           
                citydistance(m,illusion_city_code(i,j))=citydistance(illusion_city_code(i,j),m);
            else
                citydistance(illusion_city_code(i,j),m)=0;
                citydistance(m,illusion_city_code(i,j))=citydistance(illusion_city_code(i,j),m);
            end
        end
    end
end
%至此，得到一个 citydistance矩阵，网络中的三棱柱模型，没有距离，为0

% 根据实际情况，根据需求自定义不连通的城市节点，只需对cut矩阵进行增改
      cut=[];
[cut_array,cut_col]=size(cut);
for ray=1:cut_array
    cutstart=cut(ray,1);
    cutdest=cut(ray,2);
    cutstart_array=illusion_city_code(cutstart,:);
    cutdest_array=illusion_city_code(cutdest,:);%把不连通的实际城市编号取出来
    for i=1:illusion_city_col-1%只前三种方式，公路，水路，铁路，但是航空却是任意两个城市之间均有连接
        for j=1:illusion_city_col-1
            citydistance(cutstart_array(i),cutdest_array(j))=0;
            citydistance(cutdest_array(j),cutstart_array(i))=citydistance(cutstart_array(i),cutdest_array(j));
        end
    end
end
% 距离邻接矩阵去零化
for i=1:city_number
    for j=1:city_number
        if citydistance(i,j)==0
            citydistance(i,j)=inf;
            citydistance(j,i)=citydistance(i,j);
        end
    end
end
citydistance_D=citydistance;
%%


%%
%% 成本 
%邻接矩阵c
cc=eps*ones(city_number,city_number);
fprintf('基本成本是否采用默认参数，是请输入1，否请输入0\n');
c_data_choice=input('');
if c_data_choice==0    
    cost_road=input(); 
    cost_river=input(); 
    cost_railway=input(); 
    cost_fly=input(); 
    LL=input();
    WW=input();
elseif c_data_choice==1
    cost_road=40; %成本因子RMB/km
    cost_river=18; %成本因子RMB/km
    cost_railway=21; %成本因子RMB/km
    cost_fly=140;
    LL=7.5;%5;%10;%7.5;
    WW=10;   
end
    cc_truck=cost_road*LL/WW;% RMB/km，负载比重下公路成本因子
    cc_river=cost_river*LL/WW; % RMB/km，负载比重下水路成本因子
    cc_railway=cost_railway*LL/WW; % RMB/km，负载比重下铁路成本因子
    cc_fly=cost_fly*LL/WW; % RMB/km，负载比重下铁路成本因子
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            % 运输部分
            if any(b_col==m)&&(illusion_city_code(i,j)~=m)   %如果m在这一列里                
                if j==1
                    cc(m,illusion_city_code(i,j))=cc_truck;%计算两点之间的成本因子
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));
                elseif j==2
                    cc(m,illusion_city_code(i,j))=cc_river;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));
                elseif j==3
                    cc(m,illusion_city_code(i,j))=cc_railway;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));  
                elseif j==4
                    cc(m,illusion_city_code(i,j))=cc_fly;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));  
                end
             % 中转部分
            elseif any(b_array==m)&&(illusion_city_code(i,j)~=m)
                if (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==2)
                    cc(m,illusion_city_code(i,j))=cc_truck;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==3)
                    cc(m,illusion_city_code(i,j))=cc_truck;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==0)
                    cc(m,illusion_city_code(i,j))=cc_truck;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==3)
                    cc(m,illusion_city_code(i,j))=cc_truck;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==0)
                    cc(m,illusion_city_code(i,j))=cc_truck;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==3)&&(mod(m,illusion_city_col)==0)
                    cc(m,illusion_city_code(i,j))=cc_truck;
                    cc(illusion_city_code(i,j),m)=cc(m,illusion_city_code(i,j));
                end
            end
        end
    end
end

%
%cost_D=citydistance.*cc;%%****基本成本矩阵
%% ************************************************************************************
%% 时间
%平均速度 邻接矩阵v
v=eps*ones(city_number,city_number);
fprintf('速度是否采用默认参数，是请输入1，否请输入0\n');
v_data_choice=input('');
if v_data_choice==0    
    v_truck =input('请输入公路平均速度km/h：');
    v_river =input('请输入水路平均速度km/h：');
    v_railway =input('请输入铁路平均速度km/h：');
    v_fly=input('请输入航空平均速度km/h：');
    
    v_road2river = input('请输入水路-公路平均速度km/h：');
    v_road2railway =input('请输入公路-铁路平均速度km/h：');
    v_river2railway =input('请输入水路-铁路平均速度km/h：');
    v_road2fly=input('请输入公路-航空平均速度km/h：');
    v_river2fly =input('请输入水路-航空平均速度km/h：');
    v_railway2fly =input('请输入铁路-航空平均速度km/h：');
elseif v_data_choice==1    
    v_truck =90; %input('请输入公路平均速度km/h：');
    v_river =40; %input('请输入水路平均速度km/h：');
    v_railway =60; % input('请输入铁路平均速度km/h：');
    v_fly=600;
    v_road2river = 40;%input('请输入水路-公路平均速度km/h：');
    v_road2railway = 30;%input('请输入公路-铁路平均速度km/h：');
    v_river2railway =50;%input('请输入水路-铁路平均速度km/h：'); 
    v_road2fly=50;
    v_river2fly =50;
    v_railway2fly =50;
end
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            % 运输部分
            if any(b_col==m)&&(illusion_city_code(i,j)~=m)   %如果m在这一列里                
                if j==1
                    v(m,illusion_city_code(i,j))=v_truck;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                elseif j==2
                    v(m,illusion_city_code(i,j))=v_river;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                elseif j==3
                    v(m,illusion_city_code(i,j))=v_railway;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                elseif j==4
                    v(m,illusion_city_code(i,j))=v_fly;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                end
             % 中转部分
            elseif any(b_array==m)&&(illusion_city_code(i,j)~=m)
                if (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==2)
                    v(m,illusion_city_code(i,j))=v_road2river;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==3)
                    v(m,illusion_city_code(i,j))=v_road2railway;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==3)
                    v(m,illusion_city_code(i,j))=v_river2railway;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==0)
                    v(m,illusion_city_code(i,j))=v_road2fly;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==0)
                    v(m,illusion_city_code(i,j))=v_river2fly;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==3)&&(mod(m,illusion_city_col)==0)
                    v(m,illusion_city_code(i,j))=v_railway2fly;
                    v(illusion_city_code(i,j),m)=v(m,illusion_city_code(i,j));
                end
            end
        end
    end
end
%time_D=citydistance./v;%基本时间矩阵
%% ************************************************************************************************
%% 排放 
emission=inf*ones(city_number,city_number);% 任何两个节点之间某种运输方式的排放因子，g/km，或者g/h
e=inf*ones(city_number,city_number);% 单位g，任意两个城市节点之间某种运输方式的排放总量
fprintf('排放是否采用默认参数，是请输入1，否请输入0\n');
p_data_choice=input('');
if p_data_choice==0    
    ef=input('请输入公路排放因子g/L：'); %g/L
    fc=input('请输入公路燃料消耗L/km：');
    LL=input('请输入实际负载t吨：');
    WW=input('请输入理论负载t吨：');
    cf=input('请输入实际每千万时功率排放g/kwh：');
    ep=input('请输入动力机车实际功率kw：');
    ep_fly=input('请输入飞机引擎功率kw：');
    eprice =input('请输入每克排放单价RMB/g ：'); 
elseif p_data_choice==1
    ef=500;%input('请输入公路排放因子：'); g/L
    fc=1;%input('请输入公路燃料消耗：'); L/km
    LL=7.5;%5;%10;%7.5; %input('请输入实际负载：');% t吨
    WW=10;  %input('请输入理论负载：');% t吨
    cf=10;%input('请输入实际每千万时功率排放：'); g/kwh
    ep=1000;%input('请输入实际功率：');kw
    ep_fly=10000;%input('请输入航空引擎实际功率：');kw
    eprice=0.1;%input('请输入每克排放单价：')排放单价 RMB/g    
end
    ee_truck=fc*ef*LL/WW;% g/km，公路排放因子
    ee_river=cf*ep*LL/WW; % g/h，水路排放因子
    ee_railway=cf*ep*LL/WW; % g/h，铁路排放因子
    ee_fly=cf*ep_fly*LL/WW;% g/h ,航空排放因子
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            % 运输部分
            if any(b_col==m)&&(illusion_city_code(i,j)~=m)   %如果m在这一列里                
                if j==1
                    emission(m,illusion_city_code(i,j))=ee_truck;%计算两点之间的排放因子
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%计算两点之间的排放
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif j==2
                    emission(m,illusion_city_code(i,j))=ee_river;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_river.*(citydistance(m,illusion_city_code(i,j))./v(m,illusion_city_code(i,j)));
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif j==3
                    emission(m,illusion_city_code(i,j))=ee_railway;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_railway.*(citydistance(m,illusion_city_code(i,j))./v(m,illusion_city_code(i,j)));
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif j==4
                    emission(m,illusion_city_code(i,j))=ee_fly;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_fly.*(citydistance(m,illusion_city_code(i,j))./v(m,illusion_city_code(i,j)));
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                end
             % 中转部分
            elseif any(b_array==m)&&(illusion_city_code(i,j)~=m)
                if (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==2)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%计算两点之间的排放
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==3)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%计算两点之间的排放
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==3)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%计算两点之间的排放
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==0)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%计算两点之间的排放
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==0)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%计算两点之间的排放
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==3)&&(mod(m,illusion_city_col)==0)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%计算两点之间的排放
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                end
            end
        end
    end
end
%e_D=e;% 基本排放矩阵
%cost_e_D=eprice.*e_D;% 基本排放成本矩阵
%% **************************************************************************************************
%% 能耗
oil=inf*ones(city_number,city_number);
fprintf('能耗是否采用默认参数，是请输入1，否请输入0\n');
p_data_choice=input('');
if p_data_choice==0    
    fc=input('请输入公路燃料消耗L/km：');
    fc_river=input('请输入水路燃料消耗L/km：');
    fc_railway=input('请输入铁路燃料消耗L/km：');
    fc_fly=input('请输入航空燃料消耗L/km：');
    LL=input('请输入实际负载t吨：');
    WW=input('请输入理论负载t吨：');
    oprice =input('请输入每升单价RMB/L ：'); 
elseif p_data_choice==1
    fc=1;%input('请输入公路燃料消耗：'); L/km
    fc_river=0.5;%input('请输入水路燃料消耗：'); L/km
    fc_railway=0.5;%input('请输入铁路燃料消耗：'); L/km
    fc_fly=10;%input('请输入航空燃料消耗：'); L/km
    LL=7.5;%5;%10;%7.5; %input('请输入实际负载：');% t吨
    WW=10;  %input('请输入理论负载：');% t吨
    oprice =7.8; %RMB/L
end
    oil_truck=fc*LL/WW;% L/km，公路能耗因子
    oil_river=fc_river*LL/WW; % L/km，水路能耗因子
    oil_railway=fc_railway*LL/WW; % L/km，铁路能耗因子
    oil_fly=fc_fly*LL/WW; % L/km，铁路能耗因子
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            % 运输部分
            if any(b_col==m)&&(illusion_city_code(i,j)~=m)   %如果m在这一列里                
                if j==1
                    oil(m,illusion_city_code(i,j))=oil_truck;%计算两点之间的能耗因子
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));
                elseif j==2
                    oil(m,illusion_city_code(i,j))=oil_river;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));
                elseif j==3
                    oil(m,illusion_city_code(i,j))=oil_railway;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));        
                elseif j==4
                    oil(m,illusion_city_code(i,j))=oil_fly;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j)); 
                end
             % 中转部分
            elseif any(b_array==m)&&(illusion_city_code(i,j)~=m)
                if (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==2)
                    oil(m,illusion_city_code(i,j))=oil_truck;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==3)
                    oil(m,illusion_city_code(i,j))=oil_truck;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==3)
                    oil(m,illusion_city_code(i,j))=oil_truck;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==0)
                    oil(m,illusion_city_code(i,j))=oil_truck;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==0)
                    oil(m,illusion_city_code(i,j))=oil_truck;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==3)&&(mod(m,illusion_city_col)==0)
                    oil(m,illusion_city_code(i,j))=oil_truck;
                    oil(illusion_city_code(i,j),m)=oil(m,illusion_city_code(i,j));
                end
            end
        end
    end
end
%o_D=citydistance.*oil;%基本油耗矩阵
%cost_o_D=oprice.*o_D;%基本油耗成本

%% 菜单选择
choice = input('welcome!\nplease choose\n1.min cost path without considering emission\n2.min time\n3.min emission\n4.min oil\n5.min cost of basic cost,oil and emission\n6.min cost,oil and emission and analyse of eprice\n7.block to time\n8.velocity to time\n9.velocity to emission\n10.dynamicn\n');

% 成本最少，不考虑 排放，能耗
if choice==1 
    fprintf('min cost path without considering emission\n');
    cost_D=citydistance.*cc;%%****基本成本矩阵
    time_D=citydistance./v;%基本时间矩阵
    emission_D=e;% 基本排放矩阵
    cost_e_D=eprice.*emission_D;% 基本排放成本矩阵
    oil_D=citydistance.*oil;%基本油耗矩阵
    cost_o_D=oprice.*oil_D;%基本油耗成本
        cost_c_e_o_D=cost_D;  %基本成本+能耗成本
        D=cost_c_e_o_D;
      

end



