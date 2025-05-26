
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.optim.lr_scheduler import CosineAnnealingLR
import warnings
warnings.filterwarnings('ignore')
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
# ����GPU�����Ƿ��ʹ��
print(torch.__version__) # pytorch�汾
print(torch.version.cuda) # cuda�汾
print(torch.cuda.is_available()) # �鿴cuda�Ƿ����

Numm = 1
Num = Numm * 2 #��������
Epoch = 200
N_samples = 10000
m = 1  # ����
a = 1  # ʱ�䲽��
omega = 1  # ��Ƶ��

N = Num
mean2 = np.array(np.zeros(N))
Sigma_inv2 = np.zeros([N,N])
for i in range(N):
     for j in range(N):
        diff = abs(i-j)
        if diff==0:
            Sigma_inv2[i][j]=2/a+a
        elif diff==1 or diff==N-1:
            Sigma_inv2[i][j]=(-2/a)
#print(Sigma_inv2)
Sigma2 = np.linalg.inv(Sigma_inv2)
vec = np.linalg.eigvals(Sigma_inv2)
#samples2 = np.random.multivariate_normal(mean2, Sigma2, N_samples)
print(Sigma2[0, 0])
Z = -0.5*np.log(np.linalg.det(Sigma2))
print(Z)
#����ֵ
#loss����ֵ
class flow_unit(nn.Module):
    def __init__(self,dim_channels) : 
        super().__init__() #�̳�

        #��������
        self.scale = nn.Sequential( #ģ���װ
            nn.Linear(1, dim_channels),
            nn.ReLU(), #ȡ��
            nn.Linear(dim_channels, dim_channels),#ȫ���Ӳ㣬ʵ�����Ա任
            nn.ReLU(),
            nn.Linear(dim_channels, 1),
            torch.nn.Tanh() #˫�����У������Ա任
        )

        #����ƽ��
        self.trans = nn.Sequential(
            nn.Linear(1, dim_channels),
            nn.ReLU(),
            nn.Linear(dim_channels, dim_channels),
            nn.ReLU(),
            nn.Linear(dim_channels, 1),
            torch.nn.Tanh()

        )

    def forward(self, x, update_type):
        update = update_type % 2
        x_update = x[:, update::2].unsqueeze(1)
        x_used = x[:, (1 - update)::2].unsqueeze(1)
        
        scale = self.scale(x_used)
        trans = self.trans(x_used)
        x_update = x_update * torch.exp(scale) + trans


        y = torch.tensor(x)
        y[:,update::2] = x_update.squeeze(1)
        y[:,(1 - update)::2] = x_used.squeeze(1)
        log_det = -scale

        return y, log_det
# ������ģ���б����Ż����б�
flow_list = []
optim_list = []
#lrs_list = []
for i in range(30):  
    flow = flow_unit(dim_channels = 10)
    flow.cuda()
    optimizer = torch.optim.Adam(flow.parameters())
    flow_list.append(flow)
    optim_list.append(optimizer)    
    
 #   lrs = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=Epoch)
 #   lrs_list.append(lrs)

loss_list= []
# ѵ��ѭ��
for epoch in range(Epoch): #����ÿ��ѵ�����ڣ�Epoch ��һ��Ԥ�ȶ������������ʾѵ������������
    samples0 = torch.randn(size=[N_samples , Num]).cuda()  #[������������������]<-�����������
    log_det = 0 #��������ʽ
    samples = samples0
    for i in range(Num):  
        samples, log_det0 =  flow_list[i](samples, i)
        log_det = log_det + log_det0

# ������ʧ����������ģ�Ͳ���...
    energy = m * (((samples - torch.roll(samples, -1, dims = 1)).pow(2)).sum(1)) / (2 * a) + a * m * (omega**2) * (samples * samples).sum(1) / 2
#����ʱ�䲽����ô����

# ������ʧ����
    loss = (-0.5 * (samples0 * samples0).sum(1) + log_det + energy).mean()


    for i in range(Num):
        optim_list[i].zero_grad()
    loss.backward() #������ʧ���������ģ�Ͳ������ݶ�
    for i in range(Num):
        optim_list[i].step() #�����ݶȺ�ѧϰ�ʵȲ�������ģ�͵Ĳ���
#        lrs_list[i].step()
        
    #optimizer.zero_grad()
    loss_list.append(loss.cpu().detach().numpy() )
    #loss.backward()
    #optimizer.step()
print(loss_list)
plt.plot(np.arange(np.array(loss_list).shape[0]),np.array(loss_list))
(samples*samples).mean()
print(samples)