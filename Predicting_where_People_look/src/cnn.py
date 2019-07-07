import torch.nn as nn


class ConvLayer(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size, stride, padding):
        super(ConvLayer, self).__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size, stride, padding),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True)
        )

    def forward(self, x):
        return self.conv(x)


class ConvNet(nn.Module):
    def __init__(self):
        super(ConvNet, self).__init__()
        self.net = nn.Sequential(
            ConvLayer(3, 64, 5, 1, 3),
            nn.MaxPool2d(3, 2, 1),

            ConvLayer(64, 128, 3, 1, 1),
            nn.MaxPool2d(3, 2, 1),

            ConvLayer(128, 256, 3, 1, 1),
            nn.MaxPool2d(3, 2, 1),

            ConvLayer(256, 128, 3, 1, 1),
            ConvLayer(128, 64, 3, 1, 1),
            nn.Upsample((64, 64), mode='bilinear'),

            nn.Conv2d(64, 1, 3, 1, 1),
            nn.Upsample((128, 128), mode='bilinear'),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)
