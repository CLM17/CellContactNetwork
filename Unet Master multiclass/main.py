from model import unet
from data import trainGenerator, ModelCheckpoint, testGenerator, saveResult
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0"


data_gen_args = dict(rotation_range=0.2,
                    width_shift_range=0.05,
                    height_shift_range=0.05,
                    shear_range=0.05,
                    zoom_range=0.05,
                    horizontal_flip=True,
                    fill_mode='nearest')

train_folder = 'data/hela/train'
test_folder = "data/hela/test"
steps_per_epoch = 1
epochs = 300
batch_size = 2
num_test_image = 1

myGene = trainGenerator(batch_size,train_folder,'image','label',data_gen_args, flag_multi_class = True, num_class = 3, save_to_dir = False)

#model = unet(pretrained_weights='unet_membrane.hdf5') #uncomment to use pretrained weights
model = unet()

model_checkpoint = ModelCheckpoint('hela_test_network.hdf5', monitor='loss',verbose=1, save_best_only=True)
model.fit_generator(myGene,steps_per_epoch=steps_per_epoch,epochs=epochs,callbacks=[model_checkpoint])

testGene = testGenerator(test_folder,num_image=num_test_image)
results = model.predict_generator(testGene,num_test_image,verbose=1)
saveResult(test_folder,results)