function ACOfirst()
clear;
clc;
%author Xing Peng
%@copyright reserved
%% ���Ծ���D
D=[0 1 3 0;   1 0 1 2;   3 1 0 5;   0 2 5 0;  ];
%% get �ڽӾ���
%D=getD();
%% ������ͨ��ֵΪ0
N=size(D,1);
for i=1:N
    for j=1:N
        if D(i,j)==inf
            D(i,j)=0;
        end
    end
end
%%
%% ������յ�
S=1 ;         % S ��ʼ�㣨���·������ʼ�㣩
E=N;         % E ��ֹ�㣨���·����Ŀ�ĵ㣩
%% 
Tau=ones(N,N);% Tau ��ʼ��Ϣ�ؾ�����Ϊǰ�����ʳ����в�������Ϣ�أ�
Tau=8.*Tau;
K=10;          % K ����������ָ���ϳ������ٲ���
M=5;         % M ���ϸ�����ÿһ�������ж��ٸ���
Alpha=2;    % Alpha ������Ϣ����Ҫ�̶ȵĲ���
Beta=5;     % Beta ��������ʽ������Ҫ�̶ȵĲ���
Rho=0.1 ;     % Rho ��Ϣ������ϵ��
Q=1.2;          % Q ��Ϣ������ǿ��ϵ��
minkl=inf;
mink=0;
minl=0;
dlmwrite('distance.txt',D,'delimiter',' ') 
%xlswrite('d',D,sheet,range);%���excel��
%% ����ʽ��Ϣ��ȡΪ��ǰ��������Ŀ���յ��ֱ�߾���ĵ���
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
ROUTES=cell(K,M);%��ϸ���ṹ�洢ÿһ����ÿһֻ���ϵ�����·��
PL=zeros(K,M);%�þ���洢ÿһ����ÿһֻ���ϵ�����·�߳���
%% -----------����K��������ʳ���ÿ���ɳ�Mֻ����--------------------
for k=1:K
    for m=1:M
        %% ��һ����״̬��ʼ��
        W=S;%��ǰ�ڵ��ʼ��Ϊ��ʼ�㣬WΪ�������ڽڵ㣬
        Path=S;%����·�߳�ʼ��
        PLkm=0;%����·�߳��ȳ�ʼ��
        TABUkm=ones(N,1);%���ɱ��ʼ�������ж�Ϊ1����ʾû���߹�
        TABUkm(S,1)=0;%�Ѿ��ڳ�ʼ���ˣ����Ҫ�ų�
        DD=D;%�ڽӾ����ʼ��
        %% �ڶ�������һ������ǰ���Ľڵ�
        DW=DD(W,:); % ѡ����W���㣬WΪ��ͼ���룬���ڽӾ����е�W�и�DW
                DW1=find(DW);% �ҳ�DW��һ���У���Ϊ������������ڵ������У����磬W=1����DW1=2 4������˵��2,4��1����ϵ
        for j=1:length(DW1)    
            if TABUkm(DW1(j),1)==0 % DW1(1)=2, TABUkm(2)=1, TABUkm(2)��ʾ��2�еĵ�1����
                DW(DW1(j))=0;
            end
        end        
        LJD=find(DW);% 2,4
        Len_LJD=length(LJD);%��ѡ�ڵ�ĸ���,2
        iter=0;
        %% ��ʳֹͣ����������δ����ʳ�������������ͬ����һ����һֻ����
        while W~=E&&Len_LJD>=1 %��ֹ����Ϊ��WΪ�յ㣬���ߣ���ѡ�ڵ�������ڵ���1
            iter=iter+1;
            %% ��������ת�ֶķ�ѡ����һ����ô��
            PP=zeros(Len_LJD,1); % �������ϵ�ǰ���ڽڵ����һ����ѡ�ڵ���������0��������Ϊ��ѡ�ڵ���������Ϊ1
            for i=1:Len_LJD
                PP(i)=(Tau(W,LJD(i))^Alpha)*((Eta(LJD(i)))^Beta);
            end           
            sumpp=sum(PP); %��PP��һ�еĺͣ�
            PP=PP/sumpp;%�������ʷֲ�������һ���У�ÿ��Ԫ�ش����ѡ�ڵ㱻ѡ�еĸ���            
            Pcum(1)=PP(1);%��һ�����нڵ㱻ѡ�еĸ���
            for i=2:Len_LJD
                Pcum(i)=Pcum(i-1)+PP(i);
            end
            Pcum(1:Len_LJD); %�����ۼӣ����һ��������1
            Select=find(Pcum>=rand);           
            to_visit=LJD(Select(1));% ȷ����һ��Ҫ�ߵĽڵ�            
            %% ���Ĳ���״̬���ºͼ�¼            
            Path=[Path,to_visit];%���ϴ�W�ƶ���to_visit���ڵ��������ӣ�·����������·���ڵ�����ƴ�ӣ�����ƴ��
            PLkm=PLkm+DD(W,to_visit);%·����������
            W=to_visit;%�����Ƶ���һ���ڵ�            
            for kk=1:N %����ÿһ���ڵ�
                if TABUkm(kk,1)==0 %�������ڵ��߹��ˣ�Ϊ0��ʾ�߹���
                    DD(W,kk)=0; % �����߹���ĳ��ȸ�ֵΪ0
                    DD(kk,W)=0;
                end
            end
            TABUkm(W,1)=0;%�ѷ��ʹ��Ľڵ�ӽ��ɱ���ɾ���ظ�����һ������ǰ���Ľڵ�
            DW=DD(W,:);
            DW1=find(DW);
            for j=1:length(DW1)
                if TABUkm(DW1(j),1)==0
                    DW(j)=0;
                end
            end
            LJD=find(DW);
            Len_LJD=length(LJD);%��ѡ�ڵ�ĸ���            
        end
        %Path,PLkm,%һ��һֻ���ϵ���Ϣ
        %% ���岽������ÿһ��ÿһֻ���ϵ���ʳ·�ߺ�·�߳���
        ROUTES{k,m}=Path;%ROUTES,ϸ���ṹ��ÿһ��Ԫ�ض���һ������
        if Path(end)==E %���·�����е�ĩβ���յ㣬����������Ч
            PL(k,m)=PLkm; %�ѳ�����Ϣ����k�У�m�У������k����mֻ���ϵ�·������
            if PLkm<minkl %��ʼminklΪinf,��һֻ���Ϻ󣬸���
                minkl=PLkm; %��·�����ȸ�minkl
                mink=k; % ����
                minl=m; % �ڼ�ֻ����               
            end
        else
            PL(k,m)=0;
        end
    end    
    %% ��������������Ϣ��
    Delta_Tau=zeros(N,N);%��������ʼ��
    for m=1:M
        if PL(k,m)~=0 %�����·��������Ч
            ROUT=ROUTES{k,m}; %���и�·��
            TS=length(ROUT)-1;%����
            PL_km=PL(k,m);%���ȸ�out
            for s=1:TS
                x=ROUT(s);
                y=ROUT(s+1);
                Delta_Tau(x,y)=Delta_Tau(x,y)+Q/PL_km;
                Delta_Tau(y,x)=Delta_Tau(y,x)+Q/PL_km;
            end
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau;%��Ϣ�ػӷ�һ���֣�������һ����
end
%% ������
Path
PLkm

%% ����������
minPL=zeros(K);%K�ǵ���������
for i=1:K
     PLK=PL(i,:);%PL���þ���洢ÿһ����ÿһֻ���ϵ�����·�߳��ȣ��ѵ�i�ε���������·�߳��ȸ�PLK
     Nonzero=find(PLK);%�ҳ�PLK�в�Ϊ0������λ��
     PLKPLK=PLK(Nonzero);
     minPL(i)=min(PLKPLK);
end
figure(1)
plot(minPL);
hold on
grid on
title('�������ߣ���С·�����ȣ�');
xlabel('��������');
ylabel('·������');



%% ---------------------------��ͼ--------------------------------

%% �ڽӾ���D������
function D=getD()
%% ������нڵ�����ڽӾ���citydistance
% ��չ��������
city=[
    % ����վ  �ۿ�    ��վ     ����
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

% ��ʼ���о���distance����,��ÿ�����ڵ�֮���о���
citydistance=[city_number,city_number];
for i=1:city_number
    for j=1:city_number   
        citydistance(i,j)=((city(i,1)-city(j,1))^2+(city(i,2)-city(j,2))^2)^0.5;
        citydistance(j,i)=citydistance(i,j);                        %�Գƾ���
    end
end

% Ĭ�ϲ���ͨ�ĳ��нڵ�,illusion_city�����ʾ�����������ģ�������ر��룬�ݶ�Ϊ��(��·��ˮ·����·������)�������޸ģ����ӣ�ֻ��Ծ����������
% ���ڴ˴���ĳһ�����еĹ�·��Ĭ���޷�ֱ����ͨ����һ�����е���·����ˮ·
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
[illusion_city_array,illusion_city_col]=size(illusion_city_code);%ȡ����Ĵ�С������������
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
%���ˣ��õ�һ�� citydistance���������е�������ģ�ͣ�û�о��룬Ϊ0

% ����ʵ����������������Զ��岻��ͨ�ĳ��нڵ㣬ֻ���cut�����������
      cut=[];
[cut_array,cut_col]=size(cut);
for ray=1:cut_array
    cutstart=cut(ray,1);
    cutdest=cut(ray,2);
    cutstart_array=illusion_city_code(cutstart,:);
    cutdest_array=illusion_city_code(cutdest,:);%�Ѳ���ͨ��ʵ�ʳ��б��ȡ����
    for i=1:illusion_city_col-1%ֻǰ���ַ�ʽ����·��ˮ·����·�����Ǻ���ȴ��������������֮���������
        for j=1:illusion_city_col-1
            citydistance(cutstart_array(i),cutdest_array(j))=0;
            citydistance(cutdest_array(j),cutstart_array(i))=citydistance(cutstart_array(i),cutdest_array(j));
        end
    end
end
% �����ڽӾ���ȥ�㻯
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
%% �ɱ� 
%�ڽӾ���c
cc=eps*ones(city_number,city_number);
fprintf('�����ɱ��Ƿ����Ĭ�ϲ�������������1����������0\n');
c_data_choice=input('');
if c_data_choice==0    
    cost_road=input(); 
    cost_river=input(); 
    cost_railway=input(); 
    cost_fly=input(); 
    LL=input();
    WW=input();
elseif c_data_choice==1
    cost_road=40; %�ɱ�����RMB/km
    cost_river=18; %�ɱ�����RMB/km
    cost_railway=21; %�ɱ�����RMB/km
    cost_fly=140;
    LL=7.5;%5;%10;%7.5;
    WW=10;   
end
    cc_truck=cost_road*LL/WW;% RMB/km�����ر����¹�·�ɱ�����
    cc_river=cost_river*LL/WW; % RMB/km�����ر�����ˮ·�ɱ�����
    cc_railway=cost_railway*LL/WW; % RMB/km�����ر�������·�ɱ�����
    cc_fly=cost_fly*LL/WW; % RMB/km�����ر�������·�ɱ�����
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            % ���䲿��
            if any(b_col==m)&&(illusion_city_code(i,j)~=m)   %���m����һ����                
                if j==1
                    cc(m,illusion_city_code(i,j))=cc_truck;%��������֮��ĳɱ�����
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
             % ��ת����
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
%cost_D=citydistance.*cc;%%****�����ɱ�����
%% ************************************************************************************
%% ʱ��
%ƽ���ٶ� �ڽӾ���v
v=eps*ones(city_number,city_number);
fprintf('�ٶ��Ƿ����Ĭ�ϲ�������������1����������0\n');
v_data_choice=input('');
if v_data_choice==0    
    v_truck =input('�����빫·ƽ���ٶ�km/h��');
    v_river =input('������ˮ·ƽ���ٶ�km/h��');
    v_railway =input('��������·ƽ���ٶ�km/h��');
    v_fly=input('�����뺽��ƽ���ٶ�km/h��');
    
    v_road2river = input('������ˮ·-��·ƽ���ٶ�km/h��');
    v_road2railway =input('�����빫·-��·ƽ���ٶ�km/h��');
    v_river2railway =input('������ˮ·-��·ƽ���ٶ�km/h��');
    v_road2fly=input('�����빫·-����ƽ���ٶ�km/h��');
    v_river2fly =input('������ˮ·-����ƽ���ٶ�km/h��');
    v_railway2fly =input('��������·-����ƽ���ٶ�km/h��');
elseif v_data_choice==1    
    v_truck =90; %input('�����빫·ƽ���ٶ�km/h��');
    v_river =40; %input('������ˮ·ƽ���ٶ�km/h��');
    v_railway =60; % input('��������·ƽ���ٶ�km/h��');
    v_fly=600;
    v_road2river = 40;%input('������ˮ·-��·ƽ���ٶ�km/h��');
    v_road2railway = 30;%input('�����빫·-��·ƽ���ٶ�km/h��');
    v_river2railway =50;%input('������ˮ·-��·ƽ���ٶ�km/h��'); 
    v_road2fly=50;
    v_river2fly =50;
    v_railway2fly =50;
end
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            % ���䲿��
            if any(b_col==m)&&(illusion_city_code(i,j)~=m)   %���m����һ����                
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
             % ��ת����
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
%time_D=citydistance./v;%����ʱ�����
%% ************************************************************************************************
%% �ŷ� 
emission=inf*ones(city_number,city_number);% �κ������ڵ�֮��ĳ�����䷽ʽ���ŷ����ӣ�g/km������g/h
e=inf*ones(city_number,city_number);% ��λg�������������нڵ�֮��ĳ�����䷽ʽ���ŷ�����
fprintf('�ŷ��Ƿ����Ĭ�ϲ�������������1����������0\n');
p_data_choice=input('');
if p_data_choice==0    
    ef=input('�����빫·�ŷ�����g/L��'); %g/L
    fc=input('�����빫·ȼ������L/km��');
    LL=input('������ʵ�ʸ���t�֣�');
    WW=input('���������۸���t�֣�');
    cf=input('������ʵ��ÿǧ��ʱ�����ŷ�g/kwh��');
    ep=input('�����붯������ʵ�ʹ���kw��');
    ep_fly=input('������ɻ����湦��kw��');
    eprice =input('������ÿ���ŷŵ���RMB/g ��'); 
elseif p_data_choice==1
    ef=500;%input('�����빫·�ŷ����ӣ�'); g/L
    fc=1;%input('�����빫·ȼ�����ģ�'); L/km
    LL=7.5;%5;%10;%7.5; %input('������ʵ�ʸ��أ�');% t��
    WW=10;  %input('���������۸��أ�');% t��
    cf=10;%input('������ʵ��ÿǧ��ʱ�����ŷţ�'); g/kwh
    ep=1000;%input('������ʵ�ʹ��ʣ�');kw
    ep_fly=10000;%input('�����뺽������ʵ�ʹ��ʣ�');kw
    eprice=0.1;%input('������ÿ���ŷŵ��ۣ�')�ŷŵ��� RMB/g    
end
    ee_truck=fc*ef*LL/WW;% g/km����·�ŷ�����
    ee_river=cf*ep*LL/WW; % g/h��ˮ·�ŷ�����
    ee_railway=cf*ep*LL/WW; % g/h����·�ŷ�����
    ee_fly=cf*ep_fly*LL/WW;% g/h ,�����ŷ�����
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            % ���䲿��
            if any(b_col==m)&&(illusion_city_code(i,j)~=m)   %���m����һ����                
                if j==1
                    emission(m,illusion_city_code(i,j))=ee_truck;%��������֮����ŷ�����
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%��������֮����ŷ�
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
             % ��ת����
            elseif any(b_array==m)&&(illusion_city_code(i,j)~=m)
                if (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==2)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%��������֮����ŷ�
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==3)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%��������֮����ŷ�
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==3)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%��������֮����ŷ�
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==1)&&(mod(m,illusion_city_col)==0)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%��������֮����ŷ�
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==2)&&(mod(m,illusion_city_col)==0)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%��������֮����ŷ�
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                elseif (mod(illusion_city_code(i,j),illusion_city_col)==3)&&(mod(m,illusion_city_col)==0)
                    emission(m,illusion_city_code(i,j))=ee_truck;
                    emission(illusion_city_code(i,j),m)=emission(m,illusion_city_code(i,j));
                    e(m,illusion_city_code(i,j))=ee_truck.*citydistance(m,illusion_city_code(i,j));%��������֮����ŷ�
                    e(illusion_city_code(i,j),m)=e(m,illusion_city_code(i,j));
                end
            end
        end
    end
end
%e_D=e;% �����ŷž���
%cost_e_D=eprice.*e_D;% �����ŷųɱ�����
%% **************************************************************************************************
%% �ܺ�
oil=inf*ones(city_number,city_number);
fprintf('�ܺ��Ƿ����Ĭ�ϲ�������������1����������0\n');
p_data_choice=input('');
if p_data_choice==0    
    fc=input('�����빫·ȼ������L/km��');
    fc_river=input('������ˮ·ȼ������L/km��');
    fc_railway=input('��������·ȼ������L/km��');
    fc_fly=input('�����뺽��ȼ������L/km��');
    LL=input('������ʵ�ʸ���t�֣�');
    WW=input('���������۸���t�֣�');
    oprice =input('������ÿ������RMB/L ��'); 
elseif p_data_choice==1
    fc=1;%input('�����빫·ȼ�����ģ�'); L/km
    fc_river=0.5;%input('������ˮ·ȼ�����ģ�'); L/km
    fc_railway=0.5;%input('��������·ȼ�����ģ�'); L/km
    fc_fly=10;%input('�����뺽��ȼ�����ģ�'); L/km
    LL=7.5;%5;%10;%7.5; %input('������ʵ�ʸ��أ�');% t��
    WW=10;  %input('���������۸��أ�');% t��
    oprice =7.8; %RMB/L
end
    oil_truck=fc*LL/WW;% L/km����·�ܺ�����
    oil_river=fc_river*LL/WW; % L/km��ˮ·�ܺ�����
    oil_railway=fc_railway*LL/WW; % L/km����·�ܺ�����
    oil_fly=fc_fly*LL/WW; % L/km����·�ܺ�����
for i=1:illusion_city_array
    for j=1:illusion_city_col
        b_array=illusion_city_code(i,:);
        b_col=illusion_city_code(:,j);
        for m=1:city_number
            % ���䲿��
            if any(b_col==m)&&(illusion_city_code(i,j)~=m)   %���m����һ����                
                if j==1
                    oil(m,illusion_city_code(i,j))=oil_truck;%��������֮����ܺ�����
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
             % ��ת����
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
%o_D=citydistance.*oil;%�����ͺľ���
%cost_o_D=oprice.*o_D;%�����ͺĳɱ�

%% �˵�ѡ��
choice = input('welcome!\nplease choose\n1.min cost path without considering emission\n2.min time\n3.min emission\n4.min oil\n5.min cost of basic cost,oil and emission\n6.min cost,oil and emission and analyse of eprice\n7.block to time\n8.velocity to time\n9.velocity to emission\n10.dynamicn\n');

% �ɱ����٣������� �ŷţ��ܺ�
if choice==1 
    fprintf('min cost path without considering emission\n');
    cost_D=citydistance.*cc;%%****�����ɱ�����
    time_D=citydistance./v;%����ʱ�����
    emission_D=e;% �����ŷž���
    cost_e_D=eprice.*emission_D;% �����ŷųɱ�����
    oil_D=citydistance.*oil;%�����ͺľ���
    cost_o_D=oprice.*oil_D;%�����ͺĳɱ�
        cost_c_e_o_D=cost_D;  %�����ɱ�+�ܺĳɱ�
        D=cost_c_e_o_D;
      

end



