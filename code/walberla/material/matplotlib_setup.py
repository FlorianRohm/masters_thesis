import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from IPython import get_ipython
ipython = get_ipython()


# Matplotlib Setup
ipython.magic("matplotlib inline")                      # Show plots as images embedded in iPython notebook
#plt.style.use('ggplot')                                 # choose nice style
matplotlib.rcParams['figure.figsize'] = (10.0, 7.0 )    # increase size of plot images


# Code for embedding matplotlib animations as videos in iPython notebook

from tempfile import NamedTemporaryFile
import base64

VIDEO_TAG = """<video controls>
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""

def anim_to_html(anim, fps):
    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            anim.save(f.name, fps=fps, extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p', '-profile:v', 'baseline', '-level','3.0' ] )
            video = open(f.name, "rb").read()
        anim._encoded_video = base64.b64encode(video).decode('ascii')           
    
    return VIDEO_TAG.format(anim._encoded_video)


from IPython.display import HTML

def display_animation(anim, fps=30, show=True):
    plt.close(anim._fig)
    res = anim_to_html(anim,fps)
    if show:
        return HTML(res)
    else:
        return HTML("")

def makeImshowAnimation( grid, gridUpdateFunction, frames=90, **kwargs ):
    from functools import partial
    fig = plt.figure()
    im = plt.imshow( grid, interpolation='none' )
    def updatefig(*args, **kwargs):
        image = kwargs['image']
        image = gridUpdateFunction( image )
        im.set_array( image )
        return im,
    return animation.FuncAnimation(fig, partial(updatefig,image=grid), frames=frames )
