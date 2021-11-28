# Copyright (c) 2021 Pedro B. Castro
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
from scipy.constants import golden

class Mplibcfg(object):
    def __init__(self, fs=33, mscale=2.5,lenscale=2.7,lscale=1,lgndscale=1,ncol=1, extras={}, **kwargs):
        """

        :param fs: fontsize for label
        :param mscale: marker scale
        :param lenscale: lenght scale of the ticks.
        :param lscale: label scale
        :param lgndscale: lenght scale
        :param ncol: number of cols
        :param extras: extra params
        :param kwargs:
        """
        label_args = {
        'fontsize': fs,
        'fontname': 'Arial',
        }

        plt_args = {
        'markersize': label_args['fontsize']/pow(golden,mscale),
        'mec': 'black'
        }
        tick_args = {
        'axis': 'both',
        'direction': 'in',
        'right': True,
        'left': True,
        'top': True,
        'bottom': True,
        'length': label_args['fontsize']/pow(golden, lenscale),
        'labelsize': label_args['fontsize']/pow(golden,lscale)
        }
        legend_args = {
        'fontsize': label_args['fontsize']/pow(golden,lgndscale),
        'ncol': ncol,
        'frameon': False,
        'fancybox': False
        }
        self.fontsize=fs
        self.mscale=mscale
        self.lenscale=lenscale
        self.labelscale=lscale
        self.legendscale=lgndscale
        self.ncol=ncol
        self.largs=label_args
        self.pargs=plt_args
        self.targs=tick_args
        self.lgdnargs=legend_args
        self._args = {"fs": fs, 'mscale': mscale, 'lenscale': lenscale,'lscale': lscale,'lgndscale': lgndscale,'ncol': ncol}
        self.extras = extras

    def update(self, **kwargs):
        for key, value in zip(kwargs.keys(), kwargs.values()):
            print(key, value)
            self._args[key] = value
        self.__init__(extras=self.extras,**self._args)
    def get_labelargs(self):
        return self.largs
    def get_plotargs(self):
        return self.pargs
    def get_tickargs(self):
        return self.targs
    def get_legendargs(self):
        return self.lgdnargs
    def add_extra_params(self, param_dic, name):
        self.extras.update({name: param_dic})
    def get_extra(self, name):
        try:
            return self.extras[name]
        except:
            return "No param with" + name + "found in the extra params"
