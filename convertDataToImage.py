from PIL import Image
import numpy
import os
space = numpy.load('result.space.np.npy')
x_size = space.shape[0]
if not os.path.exists('img/'):
    os.makedirs('img/')
for i in range(x_size):
    im = Image.fromarray(space[i]*255, 'L')
    im.save(f'img/{i}.png')
