# -*- coding: utf-8 -*-
"""
Created on Thu Feb 09 14:39:41 2017
This is free software, under the terms of the GNU General Public License
Version 3 (www.gnu.org/licenses) without any implied warranty of
merchantability or fitness for a particular purpose.
@author: king.r.paul@gmail.com
"""

import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
#import datetime
#import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
#import matplotlib.cbook as cbook
#import pylab as pl

# pylint: disable=W0142
# pylint: disable=W0702
# pylint: disable=W0122
# pylint: disable=W0406

class GeigieLog():
    """ 
    Analyze bGeigie measurement log file. 
    https://github.com/Safecast/bGeigieNanoKit/wiki/Nano-Operation-Manual
    """
    def __init__(self, log_file):
        """
        Arguments:
            log_file -- str( path and filename bGeigie measurement log files )
        """
        # file contents -> DataFrame

        names = ['hdr', 'id', 'date', '1min', '5sec', 'counts', 'flag', 'lat', 
                 'hemi', 'long', 'e_w', 'alt', 'gps', 'sats', 'hdop', 'chks']

        extracted = []
        with open(log_file) as data:
            for line in data.readlines():
                if line[0] != '#':
                    extracted.append(line.replace('*',',').split(','))

        self.logfile = log_file

        self.dataset = pd.DataFrame(extracted, columns=names)
        self.dataset['date'] = pd.to_datetime(self.dataset['date'])
        self.dataset['1min'] = self.dataset['1min'].astype(int)
        self.dataset['5sec'] = self.dataset['5sec'].astype(int)
        self.dataset['counts'] = self.dataset['counts'].astype(int)

    def __len__(self):
        """ total number of data points """
        return len(self.dataset)

    def elapsed(self):
        """ elapsed time in the dataset """
        dates =  list(self.dataset['date'])
        return (dates[-1] - dates[0])

    def otsu(self):
        """ www.labbookpages.co.uk/software/imgProc/otsuThreshold.html """

        data = list(self.dataset['1min'])
        hist = np.histogram(data, bins = range(max(data)+1))[0]
        total = len(self)

        current_max_variance, threshold = 0, 0
        sum_total, sum_foreground, sum_background = 0, 0, 0

        for i in range(0, len(hist)):
            sum_total += i * hist[i]

        weight_background, weight_foreground = 0, 0
        variance_between, mean_background, mean_foreground = 0, 0, 0

        for i in range(len(hist)):        
            weight_background += hist[i]
            weight_foreground = total - weight_background        
            if weight_foreground == 0: 
                break
            sum_background += i*hist[i]
            sum_foreground = sum_total - sum_background
            mean_background = sum_background / weight_background
            mean_foreground = sum_foreground / weight_foreground
            variance_between = weight_background * weight_foreground
            variance_between *= (mean_background-mean_foreground)**2
            if variance_between > current_max_variance:
                current_max_variance = variance_between
                threshold = i 
        return threshold

    def binarize(self):
        """  """
        threshold = self.otsu()
        data = list(self.dataset['1min'])
        binarized = []
        for item in data:
            if item >= threshold:
                binarized.append(1)
            else:
                binarized.append(0)
        print binarized
        fig,axes = plt.subplots(1, 1, sharex=True, sharey=False)
        axes.plot(binarized, 'k-', linewidth=0.1)
        plt.show()



    plot_label_style = {'fontsize':12, 
                        'fontstyle':'italic', 
                        'weight':'ultralight'}

    def draw_plot(self):
        """ draw plots of raw count data vs time """

        fig,axes = plt.subplots(3, 1, sharex=True, sharey=False)
        fig.suptitle('bGeigie Log File Raw Data', **self.plot_label_style)

        axes[0].plot(self.dataset['1min'], 'k-', linewidth=0.1)
        axes[0].set_ylabel('1min', **self.plot_label_style)

        axes[1].plot(self.dataset['5sec'], 'k-', linewidth=0.05)
        axes[1].set_ylabel('5sec', **self.plot_label_style)
        
        axes[2].plot(self.dataset['counts'], 'k-', linewidth=1)
        axes[2].set_ylabel('accum', **self.plot_label_style)

        plt.xlabel('Filename: ' + self.logfile.split()[-1].split('.')[0], 
                   **self.plot_label_style)

        plt.subplots_adjust(hspace=0.1)
        plt.show()
        
    def draw_histogram(self):
        """ draw plot of frequency vs count rate """
        data = list(self.dataset['1min'])
        plt.hist(self.dataset['1min'], bins = math.ceil(math.sqrt(len(data))))
        plt.xlim(0, max(data))
        plt.show()

def test():
    """ Self test for GeigieLog class """
    log_dir = os.path.join(os.getcwd(),'logs')
    log_file = os.path.join(log_dir, '2016_08_29 - Tomo_Under_Console.LOG')
    bgeigie = GeigieLog(log_file)

    print len(bgeigie), 'measurements'
    print (bgeigie.elapsed()), 'h:m:s'
    print bgeigie.elapsed().seconds/5 - len(bgeigie), 'second discrepancy'

    print "threshold is:", bgeigie.otsu()

    #bgeigie.draw_plot()
    #bgeigie.draw_histogram()
    bgeigie.binarize()


if __name__ == "__main__":
    test()
    print '--------------------'
