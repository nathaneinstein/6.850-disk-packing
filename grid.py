import numpy as np
from itertools import combinations
from matplotlib import pyplot as plt


class Grid(object):
    
    def __init__(self, k, width, height):
        """ Defines a mesh grid subdivided into ceil(width/k) x ceil(height/k) buckets
            for storing unit-radius disks.
        
            k:          bucket (grid cell) width
            width:      grid width
            height:     grid height
            vert_divs:  x coordinates of k-width vertical dividers between buckets
            horiz_divs: y coordinates of k-width horizontal dividers between buckets
        """
        self.k = k
        self.width = width
        self.height = height
        # shift bin dividers by random amount, keeping k-width between dividers
        self.vert_divs = np.arange(0, width, k)[1:] #- np.random.uniform(0,k)
        self.horiz_divs = np.arange(0, height, k)[1:] #- np.random.uniform(0,k)
        self.grid = {}
    
    
    def add_unit_disks(self, disk_geoms):
        """ Adds a list of 1+ disks (plt.Point objects) to the correct grid bucket. 
        """

        def _get_bucket_coords(centroid):
            """ Finds the row and column indices of bucket (grid cell) containing 
                the given disk. """
            x, y = centroid
            ncols = len(self.vert_divs) + 1
            nrows = len(self.horiz_divs) + 1
            col_idx = ncols - sum(x < np.asarray(self.vert_divs)) - 1
            row_idx = nrows - sum(y < np.asarray(self.horiz_divs)) - 1
            return col_idx, row_idx


        for geom in disk_geoms:
            if geom.center == (None, None):
                # clicked outside of grid
                continue

            bucket_coords = _get_bucket_coords(geom.center)

            if bucket_coords in self.grid:
                self.grid[bucket_coords].append(geom.center)
            else:
                self.grid[bucket_coords] = [geom.center]
            
  
    def _get_maximal_disjoint_subset(self, disks):
        """ Performs brute-force search for the maximal subset of non-overlapping 
            disks. In particular, if there are n disks, it sequentially checks 
            each of the nCr possible subsets of size r, for sequentially smaller 
            values of r. 
        """

        def _any_overlapping(disks):
            """ Iterates through all nC2 pairs in passed list of disks checking for 
                pair-wise overlap. """
            for i in range(len(disks)):
                di = disks[i]
                for j in range(i+1, len(disks)):
                    dj = disks[j]
                    # disks di and dj overlap iff ||centroid1 - centroid2|| < r1 + r2
                    if np.sqrt((di[0]-dj[0])**2 + (di[1]-dj[1])**2) < 2: 
                        return True
            return False

        
        r = len(disks)
        while r > 0:
            for subset in combinations(disks, r):
                if not _any_overlapping(subset):
                    return subset
            r = r-1
            
    
    def _overlaps_bucket(self, disk):
        """ Checks whether the given disk overlaps multiple buckets or the edge of the grid. 
        """
        x, y = disk
        x_min = x-1; x_max = x+1; y_min = y-1; y_max = y+1

        # check for overlap with outer edges of grid
        if (x_min < 0) or (y_min < 0) or (x_max > self.width) or (y_max > self.height):
            return True

        # check for overlap with bucket divisions
        for div in self.vert_divs:
            if x_min < div < x_max:
                return True
        for div in self.horiz_divs:
            if y_min < div < y_max:
                return True

        return False
    
    
    def get_approx_maximal_packing(self):
        """ Returns an approximation of the maximal packing of (non-overlapping) unit 
            disks currently in the grid. 
            
            Procedure:
                For each bucket (grid cell):
                1. Remove any disks that overlap the bucket edges
                2. Perform brute-force search for maximal disjoint subset of remaining disks
        """
        max_packing = []
        all_candidates = []
        all_discards = []
        
        for bucket_idx in self.grid.keys():
            candidates = self.grid[bucket_idx].copy()
            discards = []
            # 1. remove any disks that fall in multiple buckets
            for disk in candidates:
                if self._overlaps_bucket(disk):
                    discards.append(disk)
            [candidates.remove(d) for d in discards]
            all_discards = all_discards + discards
            
            # 2. find the maximal subset of non-overlapping disks in this bucket
            if len(candidates) > 0:
                all_candidates = all_candidates + candidates
                max_packing_in_bucket = list(self._get_maximal_disjoint_subset(candidates))
                max_packing = max_packing + max_packing_in_bucket
        
        return max_packing, all_candidates, all_discards
            
    
    def demo_alg(self, fig=None, ax=None, figsize=(5,5), calc_n_optimal=True):
        
        def refresh_plot(fig, pause_len=1):
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(pause_len)
            
        if fig is None:
            fig, ax = plt.subplots(figsize=figsize)
        ax.set(xlim=(0, self.width), ylim = (0, self.height))

        orig_disk_set = [disk for bucket in list(self.grid.values()) for disk in bucket]
        max_packing, all_candidates, all_discards = self.get_approx_maximal_packing()

        disk_geoms = []
        for d in orig_disk_set:
            disk_geom = plt.Circle(d, radius=1, alpha=0.3, ec='b')
            ax.add_artist(disk_geom)
            disk_geoms.append(disk_geom)
        
        ax.vlines(self.vert_divs, ymin=0, ymax=self.height)
        ax.hlines(self.horiz_divs, xmin=0, xmax=self.width)
        refresh_plot(fig)
        
        for d in disk_geoms:
            if d.center in all_discards:
                d.set_edgecolor('goldenrod')
                d.set_facecolor('none')
                #d.set_linestyle(":")
                d.set_linewidth(1)
                d.set_alpha(0.4)
            elif d.center in max_packing:
                d.set_facecolor('tomato')
                d.set_edgecolor('none')
                d.set_alpha(0.8)
            elif d.center in all_candidates:
                d.set_facecolor('none')
                d.set_edgecolor('cornflowerblue')
                d.set_alpha(0.5)
                d.set_linewidth(2)
            else:
                raise Exception
        refresh_plot(fig)

        if calc_n_optimal:
            n_optimal = len(self._get_maximal_disjoint_subset(orig_disk_set))
        else:
            n_optimal = None
        
        plt.show()

        return len(orig_disk_set), len(max_packing), len(all_discards), n_optimal


