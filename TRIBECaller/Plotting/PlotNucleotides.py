# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue

from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D



def make_nucleotides_elements(text, x=0.0, y=0.0, width=1.0, height=1.0, color='blue',
                       font = FontProperties(family='Arial',weight=900)):
    tp = TextPath((0.0, 0.0), text, size=1, prop=font)
    bbox = tp.get_extents()
    bwidth = bbox.x1 - bbox.x0
    bheight = bbox.y1 - bbox.y0
    trafo = Affine2D()
    trafo.translate(-bbox.x0, -bbox.y0)
    trafo.scale(1 / bwidth * width, 1 / bheight * height)
    trafo.translate(x,y)
    tp = tp.transformed(trafo)
    return PathPatch(tp, facecolor=color, edgecolor=(0.0, 0.0, 0.0, 0))