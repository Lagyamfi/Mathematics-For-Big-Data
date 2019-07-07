import os
from PIL import Image
from torch.utils.data import Dataset
import torch

class SaliencyDataset(Dataset):
    def __init__(self, images_dir, maps_dir, transform=None):
        self.images_dir = images_dir
        self.maps_dir = maps_dir
        self.IDs = os.listdir(images_dir)
        self.transform = transform

    def __len__(self):
        return len(self.IDs)

    def __getitem__(self, idx):
        image = Image.open(
            os.path.join(self.images_dir, self.IDs[idx])
        ).convert('RGB')

        map = Image.open(
            os.path.join(self.maps_dir, self.IDs[idx][:-4] + '.png')
        )

        if self.transform:
            image = self.transform(image)
            map = self.transform(map)

        sample = {
            'image': image,
            'map': map #torch.round(map - 0.1)
        }

        return sample