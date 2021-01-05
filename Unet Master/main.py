from model import *
from data import *
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0"


data_gen_args = dict(rotation_range=0.2,
                    width_shift_range=0.05,
                    height_shift_range=0.05,
                    shear_range=0.05,
                    zoom_range=0.05,
                    horizontal_flip=True,
                    fill_mode='nearest')

train_folder = 'data/hela_boundaries/train'
test_folder = "data/hela_boundaries/test"
steps_per_epoch = 1
epochs = 5
batch_size = 2
num_test_image = 3
nr_color_channels = 3

myGene = trainGenerator(batch_size,train_folder,'image','label',data_gen_args,
                        image_color_mode = 'rgb',
                        save_to_dir = False)

#model = unet(pretrained_weights='unet_membrane.hdf5') #uncomment to use pretrained weights
model = unet(pretrained_weights = 'hela_boundaries_network.hdf5', input_size = (512,512,nr_color_channels))

model_checkpoint = ModelCheckpoint('hela_boundaries_network.hdf5', monitor='loss',verbose=1, save_best_only=True)
model.fit_generator(myGene,steps_per_epoch=steps_per_epoch,epochs=epochs,callbacks=[model_checkpoint])

testGene = testGenerator(test_folder,num_image=num_test_image, as_gray=False)
results = model.predict_generator(testGene, num_test_image, verbose=1)
saveResult(test_folder,results)