2.Ln=200;        %格点边长
3.L=zeros(Ln);   %格点矩阵
4.Q=120;         %总取向数
5.step_num=500;  %MC总步数
6.interval_save_jpg=20;         %图形存储间隔
7.interval_stastics=2;          %晶粒平均参数和相对密度统计间隔
8.stastics_data=zeros(step_num/interval_stastics,5);           %存储每interval_stastics次MCS后的平均晶粒尺寸和相对密度，存储格式为(MCS,grain count,average area,average diameter,relative density)
9.
10.
11.%烧结模拟过程参数赋值
12.T=1;    %温度参数
13.J1=1;   %晶界能
14.
15.
16.%初始结构的格点赋值
17.rand_l=randperm(Ln^2);
18.for i=1:Ln^2*V_pore
19.    if rem(rand_l(i),Ln)==0
20.        x=Ln;
21.        y=fix(rand_l(i)/Ln);
22.      else
23.        x=rem(rand_l(i),Ln);
24.        y=fix(rand_l(i)/Ln)+1;
25.    end
26.    L(x,y)=-1;  %标识孔洞区域
27.end
28.for i=(Ln^2*V_pore+1):Ln^2
29.    if rem(rand_l(i),Ln)==0
30.        x=Ln;
31.        y=fix(rand_l(i)/Ln);
32.      else
33.        x=rem(rand_l(i),Ln);
34.        y=fix(rand_l(i)/Ln)+1;
35.    end
36.    rand_Q=randperm(Q);
37.    L(x,y)=rand_Q(1);  %标识晶粒区域
38.end
39.
40.temp_L=zeros(Ln+2);    %标识边界区域标示为0，便于后续处理
41.temp_L(2:Ln+1,2:Ln+1)=L;
42.L=temp_L;
43.Ln=Ln+2;               %此时L边长Ln=Ln+2
44.
45.s=[-1 -1
46.   -1  0
47.   -1  1
48.    0 -1
49.    0  1
50.    1 -1
51.    1  0
52.    1  1];     %便于随即选取所选格点周围相邻的一个格点
53.
54.
55.%开始CAS模拟
56.for step=1:step_num
57.
58.    step   %显示MCS进程
59.
60.    rand_l=randperm(Ln^2);
61.
62.    for i=1:Ln^2    %随机选取格点
63.
64.        %if rem(i,1000)==0  %显示选取格点进程
65.        %    i
66.        %end
67.
68.        if rem(rand_l(i),Ln)==0   %把rand_l(i)转换成实际坐标L(x,y)
69.            x=Ln;
70.            y=fix(rand_l(i)/Ln);
71.          else
72.            x=rem(rand_l(i),Ln);
73.            y=fix(rand_l(i)/Ln)+1;
74.        end
75.
76.        if L(x,y)~=0  %如果格点不在边界区域则开始模拟
77.
78.            %---------------------如果选取点为晶粒格点---------------------
79.            if L(x,y)~=-1  %如果选取点为晶粒格点
80.
81.                ss=s+[x y; x y; x y; x y; x y; x y; x y; x y];      %存储所选晶粒格点L(x,y)周围格点坐标
82.
83.                diff_grain_count=0;      %计数所选晶粒格点L(x,y)周围与之不同的晶粒格点数
84.                same_grain_count=0;      %计数所选晶粒格点L(x,y)周围与之相同的晶粒格点数
85.                diff_ss=zeros(8,3);      %存储所选晶粒格点L(x,y)周围与之不同的晶粒格点坐标及其取向度值
86.                for ii=1:8
87.                    if L(ss(ii,1), ss(ii,2))~=0
88.                        if L(ss(ii,1), ss(ii,2))~=L(x,y) & L(ss(ii,1), ss(ii,2))~=-1
89.                            diff_grain_count=diff_grain_count+1;
90.                            diff_ss(diff_grain_count,1)=ss(ii,1);
91.                            diff_ss(diff_grain_count,2)=ss(ii,2);
92.                            diff_ss(diff_grain_count,3)=L(ss(ii,1),ss(ii,2));
93.                        end
94.                        if L(ss(ii,1), ss(ii,2))==L(x,y)
95.                            same_grain_count=same_grain_count+1;
96.                        end
97.                        if L(ss(ii,1), ss(ii,2))==-1
98.                            pore_count=pore_count+1;
99.                        end
100.                    end
101.                end
102.                total_grain_count=diff_grain_count+same_grain_count;
103.
104.                if diff_grain_count~=0  %如果所选晶粒格点L(x,y)周围有与之不同的晶粒格点
105.
106.                    BG_energy=J1*diff_grain_count;    %BG_energy为格点所处晶界能
107.
108.                    diff_ss_1=diff_ss(1:diff_grain_count,3);      %diff_ss_1存储被选格点周围取向度值与之不同的格点取向度值
109.                    diff_ss_2=unique(diff_ss_1);                  %去除diff_ss_1中的重复元素，并存储到diff_ss_2
110.                    rand_ll=randperm(length(diff_ss_2));
111.                    temp_Q=diff_ss_2(rand_ll(1));
112.                    change_BG_energy=J1*(total_grain_count-length(find(diff_ss_1==temp_Q)));
113.
114.                    if change_BG_energy<=BG_energy
115.                        L(x,y)=temp_Q;
116.                    end
117.
118.                    if change_BG_energy>BG_energy
119.                        set_probability=rand();
120.                        if exp(-(change_BG_energy-BG_energy)/T)>=set_probability
121.                            L(x,y)=temp_Q;
122.                        end
123.                    end
124.
125.                end  %对应于如果所选晶粒格点L(x,y)周围有与之不同的晶粒格点
126.
127.            end  %对应于如果选取点为晶粒格点
128.            %---------------------如果选取点为晶粒格点---------------------
129.
130.        end  %对应于如果格点不在边界区域则开始模拟
131.
132.    end    %对应于随机选取格点
133.
134.    %========================================================后处理过程========================================================
135.
136.    %后处理1开始---------------每interval_save_jpg次MCS后存储图形矩阵---------------%
137.    if rem(step,interval_save_jpg)==0
138.        figure1=figure('visible','off','PaperPosition',[3.067 9.28 14.81 11.1],'PaperSize',[20.98 29.68]);
139.        axes1 = axes('Layer','top','YDir','reverse','Parent',figure1);
140.        axis(axes1,[0.5 Ln-2+0.5 0.5 Ln-2+0.5]);
141.        image1 = image('CData',L(2:Ln-2+1,2:Ln-2+1),'CDataMapping','scaled','XData',[1 Ln],'YData',[1 Ln],'Parent',axes1);
142.        if V_pore==0
143.            jpg_name=strcat(num2str(step),'_','Ln=',num2str(Ln-2),'_','Q=',num2str(Q),'_','T=',num2str(T),'_','J1=',num2str(J1),'.jpg');
144.         else
145.            jpg_name=strcat(num2str(step),'_','Ln=',num2str(Ln-2),'_','Q=',num2str(Q),'_','Vpore=',num2str(V_pore),'_','T=',num2str(T),'_','J1=',num2str(J1),'_','J2=',num2str(J2),'.jpg');
146.        end
147.        saveas(image1,jpg_name,'jpg');
148.    end
149.    %后处理1结束---------------每interval_save_jpg次MCS后存储图形矩阵---------------%
150.
151.
152.    %后处理2开始---------------每interval_stastics次MCS后统计平均晶粒尺寸，并存入stastics_data中---------------%
153.    if rem(step,interval_stastics)==0
154.
155.        stastics_data(step/interval_stastics,1)=step;       %存储此次step
156.
157.        temp_L=L(2:Ln-1,2:Ln-1);      %去掉标示为0的边界
158.        Q_exist=unique(temp_L);       %Q_exist存储L矩阵中仍存在的晶粒取向度Q
159.        if Q_exist(1)==-1
160.            Q_exist=Q_exist(2:length(Q_exist));
161.        end
162.        Q_length=length(Q_exist);     %Q_length为L矩阵中仍存在的晶粒取向度个数（不包括孔洞）
163.
164.        for qq=1:Q_length             %只统计具有在L矩阵中仍存在的取向度的晶粒的个数
165.            temp_L=L;
166.            temp_L(temp_L~=Q_exist(qq))=0;
167.            temp_L=bwlabel(temp_L,8);
168.            now_grainnum=max(max(temp_L));       %返回目前取向度为Q_exist(qq)的晶粒数目
169.            stastics_data(step/interval_stastics,2)=stastics_data(step/interval_stastics,2)+now_grainnum;       %累加晶粒个数
170.        end
171.
172.        stastics_data(step/interval_stastics,3)=(Ln-2)^2*(1-V_pore)/stastics_data(step/interval_stastics,2);       %统计晶粒平均面积
173.
174.        stastics_data(step/interval_stastics,4)=sqrt(stastics_data(step/interval_stastics,3));                     %统计晶粒平均直径
175.
176.    end
177.    %后处理2结束---------------每interval_stastics次MCS后统计平均晶粒尺寸，并存入stastics_data中---------------%
178.
179.
180.    %========================================================后处理过程========================================================
181.
182.end
183.%结束MCS模拟
184.
185.
186.%写入文件
187.if V_pore==0
188.    xls_name1=strcat('stastics_Ln=',num2str(Ln-2),'_','Q=',num2str(Q),'_','T=',num2str(T),'_','J1=',num2str(J1),'.xls');
189.    xls_name2=strcat('grainsize_Ln=',num2str(Ln-2),'_','Q=',num2str(Q),'_','T=',num2str(T),'_','J1=',num2str(J1),'.mat');
190.else
191.    xls_name1=strcat('stastics_Ln=',num2str(Ln-2),'_','Q=',num2str(Q),'_','Vpore=',num2str(V_pore),'_','T=',num2str(T),'_','J1=',num2str(J1),'_','J2=',num2str(J2),'.xls');
192.    xls_name2=strcat('grainsize_Ln=',num2str(Ln-2),'_','Q=',num2str(Q),'_','Vpore=',num2str(V_pore),'_','T=',num2str(T),'_','J1=',num2str(J1),'_','J2=',num2str(J2),'.mat');
193.end
194.xlswrite(xls_name1,stastics_data);
195.save(xls_name2,'stastics_grainsize');
