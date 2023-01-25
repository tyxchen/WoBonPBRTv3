import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path


def main():
    data_file = Path(sys.argv[1]).with_suffix('.png')
    output_file = data_file.with_stem(data_file.stem + '-plt')

    data = pd.read_csv(str(data_file) + '.txt', header=None).to_numpy()
    mask = pd.read_csv(str(data_file) + '.mask', header=None).to_numpy()

    # remove outermost "ring" of interior points, these are miscalculated
    mask = np.pad(mask, (1, 1), 'constant', constant_values=0)
    conv_filter = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    s = tuple(np.subtract(mask.shape, conv_filter.shape) + 1) + conv_filter.shape
    mask_sub = np.lib.stride_tricks.as_strided(mask, shape=s, strides=mask.strides * 2)
    mask_filtered = np.sum(mask_sub * conv_filter, axis=(2, 3))
    mask_filtered = mask_filtered >= 4

    data_ma = np.ma.array(data, mask=~mask_filtered)

    bg = plt.imread(data_file.with_stem(data_file.stem + '-raw'))
    bg_mask = np.repeat(np.reshape(mask_filtered, (*bg.shape[:-1], 1)), bg.shape[-1], axis=2)
    bg_ma = np.ma.array(bg, mask=bg_mask).filled(0)

    lims = [-.5, .5]  #np.percentile(data_ma.compressed(), [.5, 99.5])
    data_ma -= lims[0]
    data_ma /= lims[1] - lims[0]
    img = bg_ma + mpl.colormaps['jet'](data_ma.clip(0, 1), alpha=1)[:, :, :3]
    plt.imsave(output_file, img)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: render.py data-file")
        exit(1)
    main()
