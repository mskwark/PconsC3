import sys
import h5py

import numpy as np
import joblib

treedir = sys.argv[1]

for i in range(6):
    tlst = open(treedir + '/tlayer' + str(i) + '/tree.list').read().strip().split('\n')
    trees = []

    with h5py.File('tlayer' + str(i) + '.hdf5', "w") as f:

        for tree in tlst:
            t = joblib.load(treedir + '/tlayer' + str(i) + '/' + tree.split('/')[-1])[:5]
            grp = f.create_group(tree)
            trunks = np.vstack((t[0], t[2:4]))
            compares = t[1]
            leaf = t[4][:, 0, :]
            leafs = leaf[:, 1] / leaf.sum(axis=1)

            grp.create_dataset('trunks', data=trunks)
            grp.create_dataset('compares', data=compares)
            grp.create_dataset('leafs', data=leafs)


for i in range(6):
    for j in range(5,10):
        tlst = open(treedir + '/tlayer' + str(i) + '-' + str(j) + '/tree.list').read().strip().split('\n')
        trees = []

        with h5py.File('tlayer' + str(i) + '-' + str(j) + '.hdf5', "w") as f:

            for tree in tlst:
                t = joblib.load(treedir + '/tlayer' + str(i) + '-' + str(j) + '/' + tree.split('/')[-1])[:5]
                grp = f.create_group(tree)
                trunks = np.vstack((t[0], t[2:4]))
                compares = t[1]
                leaf = t[4][:, 0, :]
                leafs = leaf[:, 1] / leaf.sum(axis=1)

                grp.create_dataset('trunks', data=trunks)
                grp.create_dataset('compares', data=compares)
                grp.create_dataset('leafs', data=leafs)
