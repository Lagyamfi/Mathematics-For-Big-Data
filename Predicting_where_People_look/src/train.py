################
##### CUDA #####
################

import torch

use_cuda = torch.cuda.is_available()
device = torch.device("cuda" if use_cuda else "cpu")
Tensor = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor

print('\nGPU: ' + str(use_cuda) + '\n')

################
##### DATA #####
################

from torchvision import transforms
from torch.utils.data import DataLoader
from dataset import SaliencyDataset

resolution = (128, 128)
batch_size = 32

data = {
    'train': SaliencyDataset(
        '../images/train',
        '../maps/train',
        transform=transforms.Compose([
            transforms.Resize(resolution),
            transforms.ToTensor(),
        ])
    ),
    'val': SaliencyDataset(
        '../images/val',
        '../maps/val',
        transform=transforms.Compose([
            transforms.Resize(resolution),
            transforms.ToTensor(),
        ])
    )
}

loader = {
    'train': DataLoader(
        data['train'],
        batch_size=batch_size,
        shuffle=True,
        num_workers=4
    ),
    'val': DataLoader(
        data['val'],
        batch_size=batch_size,
        shuffle=True,
        num_workers=4
    )
}

#############
### MODEL ###
#############

from torch import nn, optim
from cnn import ConvNet

num_epochs = 10
learning_rate = 0.001

net = ConvNet()
net.to(device)

criterion = nn.MSELoss()
optimizer = optim.Adam(net.parameters(), lr=learning_rate, betas=(0.5, 0.999))

################
### TRAINING ###
################

import matplotlib.pyplot as plt


n_batches = len(loader['train'])
for epoch in range(num_epochs):
    for i, sample in enumerate(loader['train']):

        sample['image'] = sample['image'].to(device)
        sample['map'] = sample['map'].to(device)

        net.zero_grad()

        output = net(sample['image'])
        loss = criterion(output, sample['map'])

        loss.backward()
        optimizer.step()

        if i % 10 == 10 - 1:
            print('batch [{}/{}], epoch [{}/{}], loss:{:.4f}'.format(
                i + 1, n_batches, epoch + 1, num_epochs, loss.data
            ))

            plt.figure()

            plt.subplot(131)
            plt.imshow(sample['image'][0].cpu().permute(1, 2, 0))
            plt.title('Image')
            plt.axis('off')

            plt.subplot(132)
            plt.imshow(sample['map'][0].cpu().view(*resolution))
            plt.title('Saliency Map')
            plt.axis('off')

            plt.subplot(133)
            plt.imshow(output[0].detach().cpu().view(*resolution))
            plt.title('Prediction')
            plt.axis('off')

            plt.show()

