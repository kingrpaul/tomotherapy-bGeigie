# -*- coding: utf-8 -*-
"""
Created on Thu Feb 09 14:39:41 2017
This is free software, under the terms of the GNU General Public License
Version 3 (www.gnu.org/licenses) without any implied warranty of
merchantability or fitness for a particular purpose.
@author: king.r.paul@gmail.com
"""

import os
#import numpy as np
import pandas as pd

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

        names = ['hdr', 'id', 'date', '1min', '5sec', 'counts', 
                 'flag', 'lat', 'hemi', 'long', 'e_w', 'alt', 
                 'gps', 'sats', 'hdop', 'chks']

        extracted = []
        with open(log_file) as data:
            for line in data.readlines():
                if line[0] != '#':
                    extracted.append(line.replace('*',',').split(','))

        self.dataset = pd.DataFrame(extracted, columns=names)
        self.dataset['date'] = pd.to_datetime(self.dataset['date'])
        self.dataset['1min'] = self.dataset['1min'].astype(int)
        self.dataset['5sec'] = self.dataset['5sec'].astype(int)
        self.dataset['counts'] = self.dataset['counts'].astype(int)

    def otsu(self):
        """ Find beam-on / beam-off times. """
        pass #stub

    def make_graph(self):
        """ Make graphical output """
        import matplotlib.pyplot as plt 
        plt.plot(self.dataset['date'], self.dataset['1min'])
        plt.show()

def test():
    """ Self test for GeigieLog class """
    log_dir = os.path.join(os.getcwd(),'logs')
    log_file = os.path.join(log_dir, '2016_08_25 - Tomo_In_Long_Room.LOG')
    GeigieLog(log_file).make_graph()
    
if __name__ == "__main__":
    test()
    print '--------------------'
