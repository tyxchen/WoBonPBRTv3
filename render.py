import argparse
import subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path


class Render:
    def __init__(self, args):
        self.args = args

    def render_result(self):
        scene_file = Path(self.args.datafile)
        scene_file = scene_file.with_stem(scene_file.stem + '-scene')

        print("### Rendering scene ###")
        subprocess.run([self.args.pbrt, scene_file])

        print("### Rendering domain ###")
        subprocess.run([self.args.pbrt, self.args.datafile])

    def composite(self):
        print("### Compositing ###")
        data_file = Path(self.args.datafile).with_suffix('.png')
        output_file = data_file.with_stem(data_file.stem + '-plt')

        data = pd.read_csv(str(data_file) + '.txt', header=None).to_numpy()
        mask = pd.read_csv(str(data_file) + '.mask', header=None).to_numpy()

        # Remove outermost "ring" of interior points, these are miscalculated
        mask = np.pad(mask, (1, 1), 'constant', constant_values=0)
        conv_filter = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        s = tuple(np.subtract(mask.shape, conv_filter.shape) + 1) + conv_filter.shape
        mask_sub = np.lib.stride_tricks.as_strided(mask, shape=s, strides=mask.strides * 2)
        mask_filtered = np.sum(mask_sub * conv_filter, axis=(2, 3))
        mask_filtered = mask_filtered >= 4

        # Mask data to include only the interior
        data_ma = np.ma.array(data, mask=~mask_filtered)

        # Read in rendered scene file
        bg = plt.imread(data_file.with_stem(data_file.stem + '-scene'))
        bg_mask = np.repeat(np.reshape(mask_filtered, (*bg.shape[:-1], 1)), bg.shape[-1], axis=2)
        bg_ma = np.ma.array(bg, mask=bg_mask).filled(0)

        # Clamp data and composite
        lims = np.percentile(data_ma.compressed(), [0.5, 99.5])
        data_ma -= lims[0]
        data_ma /= lims[1] - lims[0]
        img = bg_ma + mpl.colormaps['jet'](data_ma.clip(0, 1), alpha=1)[:, :, :3]
        plt.imsave(output_file, img)

    def run(self):
        if not self.args.composite_only:
            self.render_result()
        self.composite()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('datafile', help="Path to PBRT data file")
    parser.add_argument('-c', '--composite-only', action='store_true',
                        help="Only perform compositing; don't run PBRT")
    parser.add_argument('--pbrt', help="Path to PBRT executable", default='pbrt')

    renderer = Render(parser.parse_args())
    renderer.run()
