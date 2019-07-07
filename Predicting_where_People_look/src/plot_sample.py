from dataset import SaliencyDataset

data = {
    'train': SaliencyDataset(
        '../images/train',
        '../maps/train'
    ),
    'val': SaliencyDataset(
        '../images/val',
        '../maps/val'
    ),
}

import matplotlib.pyplot as plt
import random

def plot_sample(dataset):
    idx = random.randint(0, len(dataset) - 1)

    plt.figure()

    plt.subplot(121)
    plt.imshow(dataset[idx]['image'])
    plt.title('Image')
    plt.axis('off')

    plt.subplot(122)
    plt.imshow(dataset[idx]['map'])
    plt.title('Saliency Map')
    plt.axis('off')

    plt.show()


if __name__ == '__main__':
    plot_sample(data['train'])