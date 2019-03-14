#devtools::use_package('rPython')
setwd("/Users/Albi/Dropbox/Roth Lab/projects/twas_git/scripts")
devtools::load_all('../packages/twasAnalysis')
devtools::document('../packages/twasAnalysis')


rPython::python.load(system.file('python/fitness_plot.py',package='twasAnalysis'))
#rPython::python.exec('import csv')
#rPython::python.exec('import matplotlib.pyplot as plt')
#rPython::python.exec('import numpy as np')
#rPython::python.exec('import os')
#rPython::python.exec('import pandas as pd')
#stop()

rPython:::python.call('fitness_plot')