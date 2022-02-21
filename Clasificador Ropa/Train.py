# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 21:07:40 2020

@author: Eduardo
"""

# Clasificador Ropa
import tensorflow as tf
from tensorflow import keras

import numpy as np
import matplotlib.pyplot as plt

fashion_mnist = keras.datasets.fashion_mnist

(train_images, train_labels), (test_images, test_labels) = fashion_mnist.load_data()

class_names = ['T-shirt/top', 'Trouser', 'Pullover', 'Dress', 'Coat',
               'Sandal', 'Shirt', 'Sneaker', 'Bag', 'Ankle boot']

# Normalize

train_images = train_images / 255.0
test_images = test_images / 255.0

model = keras.Sequential(
        [keras.layers.Flatten(input_shape = (28, 28)),
         keras.layers.Dense(4096, activation = 'relu'),
         keras.layers.Dense(2048, activation = 'relu'),
         keras.layers.Dense(1024, activation = 'relu'),
         keras.layers.Dense(512, activation = 'relu'),
         keras.layers.Dense(256, activation = 'relu'),
         keras.layers.Dense(128, activation = 'relu'),
         keras.layers.Dense(64, activation = 'relu'),
         keras.layers.Dense(32, activation = 'relu'),
         keras.layers.Dense(16, activation = 'relu'),
         keras.layers.Dropout(0.5),
         keras.layers.Dense(10, activation = 'softmax')])

model.compile(optimizer = 'adam',
              loss = 'sparse_categorical_crossentropy',
              metrics = ['accuracy'])

model.fit(train_images, train_labels, epochs = 10)

target_dir = './modelo/'
if not os.path.exists(target_dir):
  os.mkdir(target_dir)
model.save('./modelo/modelo.h5')
model.save_weights('./modelo/pesos.h5')

test_loss, test_acc = model.evaluate(test_images,  test_labels, verbose=2)

print('\nTest accuracy:', test_acc)

predictions = model.predict(test_images)

res_idx = np.argmax(predictions[0])

print('The first object is a: ', class_names[res_idx])
