import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path


def main():
    data_file = Path(sys.argv[1])
    output_file = Path(sys.argv[2])

    data = pd.read_csv(str(data_file) + '.txt', header=None).to_numpy()
    mask = pd.read_csv(str(data_file) + '.mask', dtype=np.float, header=None).to_numpy()

    # remove outermost "ring" of interior points, these are miscalculated
    mask = np.pad(mask, (1, 1), 'constant', constant_values=0)
    conv_filter = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    s = tuple(np.subtract(mask.shape, conv_filter.shape) + 1) + conv_filter.shape
    mask_sub = np.lib.stride_tricks.as_strided(mask, shape=s, strides=mask.strides * 2)
    mask_filtered = np.sum(mask_sub * conv_filter, axis=(2, 3))
    mask_filtered = mask_filtered >= 4

    data_ma = np.ma.array(data, mask=~mask_filtered)

    lims = np.percentile(data_ma, [0, 100])

    plt.imsave(output_file, data_ma, cmap='jet', vmin=lims[0], vmax=lims[1])


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: render.py data-file output-file")
        exit(1)
    main()
