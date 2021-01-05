from model import *
from data import *
import os
import tensorflow as tf
import pickle

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

#---------------------------SPECIFY PARAMETERS---------------------------------

# Folder specifications
train_folder = 'data/cos_soma/train'
test_folder = 'data/cos_soma/test'

# Pretrained weights? Set to None if you don't want to use pretrained weights.
pretrained_weights = None
#pretrained_weights = 'hela_boundaries.hdf5'
output = 'cos_soma_final'
output_filename = 'models/' +  output + '.hdf5'

# input image specifications
input_xsize = 512
input_ysize = 512
color_mode = 'rgb'
nr_color_channels = 3
nb_classes = 2
flag_multi_class = False

# Train and test specifications
steps_per_epoch = 1
epochs = 500
batch_size = 3
learning_rate = 5e-5
num_test_image = 10

#---------------------------------START CODE-----------------------------------

data_gen_args = dict(rotation_range=0.2,
                    width_shift_range=0.05,
                    height_shift_range=0.05,
                    shear_range=0.05,
                    zoom_range=0.05,
                    horizontal_flip=True,
                    fill_mode='nearest')

myGene = trainGenerator(batch_size,train_folder,'image','label',data_gen_args,
                        image_color_mode = color_mode,
                        flag_multi_class = flag_multi_class, num_class = nb_classes, 
                        save_to_dir = False,
                        target_size = (input_xsize, input_ysize))

#model = unet(nb_classes = nb_classes, pretrained_weights='unet_membrane.hdf5') #uncomment to use pretrained weights
model = unet(learning_rate = learning_rate, nb_classes = nb_classes, input_size = (input_xsize,input_ysize,nr_color_channels))

model_checkpoint = ModelCheckpoint(output_filename, monitor='loss',verbose=1, save_best_only=True)
history = model.fit_generator(myGene,steps_per_epoch=steps_per_epoch,epochs=epochs,callbacks=[model_checkpoint])

as_gray = True
if nr_color_channels > 1:
    as_gray = False

testGene = testGenerator(test_folder,num_image=num_test_image,
                         target_size = (input_xsize, input_ysize),
                         flag_multi_class=flag_multi_class,as_gray=as_gray)

results = model.predict_generator(testGene,num_test_image,verbose=1)
saveResult(test_folder, results,
           target_size = (input_xsize, input_ysize),
           flag_multi_class = flag_multi_class, num_class = nb_classes)

#plt.figure()
#plt.plot(history.history['loss'])
#plt.xlabel('epoch')
#plt.ylabel('loss')
#plt.savefig(output+'_loss.png')

#plt.figure()
#plt.plot(history.history['acc'])
#plt.xlabel('epoch')
#plt.ylabel('accuracy')
#plt.savefig(output+'_accuracy.png')

summary = history.history
summary['name'] = output
summary['steps_per_epoch'] = steps_per_epoch
summary['epochs'] = epochs
summary['batch_size'] = batch_size
summary['learning_rate'] = learning_rate

summary_file = 'models/' + output + "_summary.pkl"
f = open(summary_file,'wb')
pickle.dump(summary,f)
f.close()
