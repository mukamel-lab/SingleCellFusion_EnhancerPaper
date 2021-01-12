nohup: ignoring input
2 5 10 20 50 100 200 500 1000
2, 2
11/09/2020 12:42:59 PM  config_mop_smarter_cells_snmcseq_gene_ka2_knn2_201108.py
11/09/2020 12:42:59 PM * Begin integration
11/09/2020 12:42:59 PM Metadata snmcseq_gene (9366, 32)
sys:1: DtypeWarning: Columns (64) have mixed types.Specify dtype option on import or set low_memory=False.
11/09/2020 12:43:00 PM Metadata smarter_cells (6244, 129)
11/09/2020 12:43:25 PM Feature matrix snmcseq_gene (4754, 9366)
11/09/2020 12:43:26 PM Feature matrix smarter_cells (5743, 6244)
11/09/2020 12:43:26 PM Done reading data
11/09/2020 12:43:26 PM Subsampling 1/3
11/09/2020 12:43:26 PM snmcseq_gene (7446, 32) (4754, 7446)
11/09/2020 12:43:27 PM smarter_cells (4995, 129) (5743, 4995)
11/09/2020 12:43:27 PM Smoothing within modalities...
./01.run_mc.sh: line 25: 72980 Illegal instruction     (core dumped) /cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling.py -c ${scf_config} -s ${subsample_frac} -sn ${subsample_times}
11/09/2020 12:43:31 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka2_knn2_201108
11/09/2020 12:43:31 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka2_knn2_201108
11/09/2020 12:43:31 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka2_knn2_201108
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 12:43:31 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 12:43:32 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 12:43:32 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 12:43:32 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 12:43:34 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 12:43:34 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 12:43:34 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 12:43:45 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 12:43:45 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 12:43:45 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 12:44:08 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 12:44:08 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 12:44:08 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 12:44:12 PM (5030, 7453)_(5030, 5030)_(32324, 5030)_(233514, 7453)_(233514, 7453)
11/09/2020 12:44:12 PM (5029, 7451)_(5029, 5029)_(32324, 5029)_(233514, 7451)_(233514, 7451)
11/09/2020 12:44:17 PM (4995, 7444)_(4995, 4995)_(32324, 4995)_(233514, 7444)_(233514, 7444)
11/09/2020 12:44:44 PM (32324, 4995)_(233514, 4995)_(233514, 4995)
11/09/2020 12:44:47 PM (32324, 5029)_(233514, 5029)_(233514, 5029)
11/09/2020 12:44:47 PM (32324, 5030)_(233514, 5030)_(233514, 5030)
11/09/2020 12:44:50 PM Note: NumExpr detected 32 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
11/09/2020 12:44:50 PM NumExpr defaulting to 8 threads.
11/09/2020 12:44:55 PM Note: NumExpr detected 32 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
11/09/2020 12:44:55 PM NumExpr defaulting to 8 threads.
11/09/2020 12:44:56 PM Note: NumExpr detected 32 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
11/09/2020 12:44:56 PM NumExpr defaulting to 8 threads.
11/09/2020 12:46:14 PM Imputing data... (No effect if sufficient_coverage_fraction=1)
11/09/2020 12:46:15 PM (7, 4995)
/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/stats/stats.py:2500: RuntimeWarning: invalid value encountered in true_divide
  return (a - mns) / sstd
11/09/2020 12:46:19 PM 0
11/09/2020 12:46:20 PM 0
11/09/2020 12:46:20 PM 0
11/09/2020 12:46:20 PM (164,)_(164,)_(164,)
11/09/2020 12:46:20 PM 2182_994112
11/09/2020 12:46:27 PM Imputing data... (No effect if sufficient_coverage_fraction=1)
11/09/2020 12:46:28 PM (11, 5030)
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
11/09/2020 12:46:28 PM Imputing data... (No effect if sufficient_coverage_fraction=1)
11/09/2020 12:46:29 PM (13, 5029)
/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/stats/stats.py:2500: RuntimeWarning: invalid value encountered in true_divide
  return (a - mns) / sstd
11/09/2020 12:46:33 PM 0
11/09/2020 12:46:33 PM 0
11/09/2020 12:46:33 PM 0
11/09/2020 12:46:34 PM (268,)_(268,)_(268,)
11/09/2020 12:46:34 PM 2182_996043
/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/stats/stats.py:2500: RuntimeWarning: invalid value encountered in true_divide
  return (a - mns) / sstd
11/09/2020 12:46:34 PM 0
11/09/2020 12:46:35 PM 0
11/09/2020 12:46:35 PM 0
11/09/2020 12:46:35 PM (279,)_(279,)_(279,)
11/09/2020 12:46:35 PM 2182_996043
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
../correlation_analysis_celllevel_mc.py:245: RuntimeWarning: divide by zero encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:245: RuntimeWarning: invalid value encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:246: RuntimeWarning: divide by zero encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:246: RuntimeWarning: invalid value encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:287: RuntimeWarning: invalid value encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: divide by zero encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: invalid value encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:245: RuntimeWarning: divide by zero encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:246: RuntimeWarning: divide by zero encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:246: RuntimeWarning: invalid value encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:287: RuntimeWarning: invalid value encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: divide by zero encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: invalid value encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:245: RuntimeWarning: divide by zero encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:245: RuntimeWarning: invalid value encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:246: RuntimeWarning: invalid value encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:246: RuntimeWarning: divide by zero encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:287: RuntimeWarning: divide by zero encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:287: RuntimeWarning: invalid value encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: invalid value encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: divide by zero encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
5, 5
11/09/2020 01:23:10 PM  config_mop_smarter_cells_snmcseq_gene_ka5_knn5_201108.py
11/09/2020 01:23:10 PM * Begin integration
11/09/2020 01:23:10 PM Metadata snmcseq_gene (9366, 32)
sys:1: DtypeWarning: Columns (64) have mixed types.Specify dtype option on import or set low_memory=False.
11/09/2020 01:23:11 PM Metadata smarter_cells (6244, 129)
11/09/2020 01:23:33 PM Feature matrix snmcseq_gene (4754, 9366)
11/09/2020 01:23:34 PM Feature matrix smarter_cells (5743, 6244)
11/09/2020 01:23:34 PM Done reading data
11/09/2020 01:23:34 PM Subsampling 1/3
11/09/2020 01:23:34 PM snmcseq_gene (7543, 32) (4754, 7543)
11/09/2020 01:23:34 PM smarter_cells (5010, 129) (5743, 5010)
11/09/2020 01:23:34 PM Smoothing within modalities...
./01.run_mc.sh: line 25: 73593 Illegal instruction     (core dumped) /cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling.py -c ${scf_config} -s ${subsample_frac} -sn ${subsample_times}
11/09/2020 01:23:37 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka5_knn5_201108
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 01:23:37 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka5_knn5_201108
11/09/2020 01:23:37 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka5_knn5_201108
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 01:23:37 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 01:23:37 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 01:23:37 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 01:23:37 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 01:23:38 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 01:23:38 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 01:23:38 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 01:23:39 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 01:23:39 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 01:23:39 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 01:23:41 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 01:23:41 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 01:23:41 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 01:23:45 PM (4977, 7512)_(4977, 4977)_(32324, 4977)_(233514, 7512)_(233514, 7512)
11/09/2020 01:23:45 PM (4961, 7579)_(4961, 4961)_(32324, 4961)_(233514, 7579)_(233514, 7579)
Traceback (most recent call last):
  File "../correlation_analysis_celllevel_mc.py", line 132, in <module>
    knn_xy = knn_xy.tocsr()[cell_idx_xaxis,:].tocsc()[:,cell_idx_yaxis] # x-by-y
  File "/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/sparse/_index.py", line 33, in __getitem__
    row, col = self._validate_indices(key)
  File "/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/sparse/_index.py", line 137, in _validate_indices
    row = self._asindices(row, M)
  File "/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/sparse/_index.py", line 169, in _asindices
    raise IndexError('index (%d) out of range' % max_indx)
IndexError: index (5009) out of range
11/09/2020 01:24:32 PM (32324, 4977)_(233514, 4977)_(233514, 4977)
11/09/2020 01:24:32 PM (32324, 4961)_(233514, 4961)_(233514, 4961)
11/09/2020 01:24:37 PM Note: NumExpr detected 32 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
11/09/2020 01:24:37 PM NumExpr defaulting to 8 threads.
11/09/2020 01:24:37 PM Note: NumExpr detected 32 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
11/09/2020 01:24:37 PM NumExpr defaulting to 8 threads.
11/09/2020 01:25:52 PM Imputing data... (No effect if sufficient_coverage_fraction=1)
11/09/2020 01:25:53 PM (1814, 4961)
11/09/2020 01:25:54 PM Imputing data... (No effect if sufficient_coverage_fraction=1)
11/09/2020 01:25:55 PM (1784, 4977)
/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/stats/stats.py:2500: RuntimeWarning: invalid value encountered in true_divide
  return (a - mns) / sstd
11/09/2020 01:25:57 PM 0
/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/stats/stats.py:2500: RuntimeWarning: invalid value encountered in true_divide
  return (a - mns) / sstd
11/09/2020 01:25:58 PM 0
11/09/2020 01:26:04 PM 0
11/09/2020 01:26:05 PM 0
11/09/2020 01:26:11 PM 0
11/09/2020 01:26:11 PM 0
11/09/2020 01:26:17 PM (37679,)_(37679,)_(37679,)
11/09/2020 01:26:17 PM 2182_999930
11/09/2020 01:26:18 PM (37052,)_(37052,)_(37052,)
11/09/2020 01:26:18 PM 2182_999930
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
No XVisualInfo for format QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Falling back to using screens root_visual.
../correlation_analysis_celllevel_mc.py:287: RuntimeWarning: invalid value encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: divide by zero encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: invalid value encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:287: RuntimeWarning: invalid value encountered in true_divide
  fdr = cdf_shuff/cdf
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: divide by zero encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
../correlation_analysis_celllevel_mc.py:288: RuntimeWarning: invalid value encountered in true_divide
  tracks_hist_ratios[label] = hist/hist_shuff
10, 10
11/09/2020 01:26:40 PM  config_mop_smarter_cells_snmcseq_gene_ka10_knn10_201108.py
11/09/2020 01:26:40 PM * Begin integration
11/09/2020 01:26:40 PM Metadata snmcseq_gene (9366, 32)
sys:1: DtypeWarning: Columns (64) have mixed types.Specify dtype option on import or set low_memory=False.
11/09/2020 01:26:40 PM Metadata smarter_cells (6244, 129)
11/09/2020 01:27:03 PM Feature matrix snmcseq_gene (4754, 9366)
11/09/2020 01:27:03 PM Feature matrix smarter_cells (5743, 6244)
11/09/2020 01:27:03 PM Done reading data
11/09/2020 01:27:03 PM Subsampling 1/3
11/09/2020 01:27:03 PM snmcseq_gene (7523, 32) (4754, 7523)
11/09/2020 01:27:04 PM smarter_cells (4968, 129) (5743, 4968)
11/09/2020 01:27:04 PM Smoothing within modalities...
./01.run_mc.sh: line 25: 73768 Illegal instruction     (core dumped) /cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling.py -c ${scf_config} -s ${subsample_frac} -sn ${subsample_times}
11/09/2020 01:27:06 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka10_knn10_201108
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 01:27:07 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka10_knn10_201108
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 01:27:07 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 01:27:07 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 01:27:07 PM corr_analysis_smarter_cells_snmcseq_gene_mop_smarter_cells_snmcseq_gene_ka10_knn10_201108
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(6244,)_cell_smarter_cells.txt
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(9364,)_cell_snmcseq_gene.txt
11/09/2020 01:27:07 PM <class 'numpy.ndarray'>_(32324,)_gene_smarter_cells.txt
11/09/2020 01:27:07 PM <class 'pandas.core.frame.DataFrame'>_(233514, 3)_enh_snmcseq_gene.tsv
11/09/2020 01:27:07 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 01:27:07 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 01:27:07 PM <class 'scipy.sparse.csc.csc_matrix'>_(32324, 6244)_mat_smarter_cells.npz
11/09/2020 01:27:08 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 01:27:08 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 01:27:08 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_mcg_snmcseq_gene.npz
11/09/2020 01:27:10 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 01:27:10 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 01:27:10 PM <class 'scipy.sparse.csc.csc_matrix'>_(233514, 9364)_mat_cg_snmcseq_gene.npz
11/09/2020 01:27:14 PM (4982, 7536)_(4982, 4982)_(32324, 4982)_(233514, 7536)_(233514, 7536)
11/09/2020 01:27:14 PM (4968, 7475)_(4968, 4968)_(32324, 4968)_(233514, 7475)_(233514, 7475)
Traceback (most recent call last):
  File "../correlation_analysis_celllevel_mc.py", line 132, in <module>
    knn_xy = knn_xy.tocsr()[cell_idx_xaxis,:].tocsc()[:,cell_idx_yaxis] # x-by-y
  File "/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/sparse/_index.py", line 33, in __getitem__
    row, col = self._validate_indices(key)
  File "/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/sparse/_index.py", line 146, in _validate_indices
    col = self._asindices(col, N)
  File "/cndd2/fangming/venvs/routine/lib/python3.8/site-packages/scipy/sparse/_index.py", line 169, in _asindices
    raise IndexError('index (%d) out of range' % max_indx)
IndexError: index (7522) out of range
11/09/2020 01:28:22 PM (32324, 4982)_(233514, 4982)_(233514, 4982)
11/09/2020 01:28:23 PM (32324, 4968)_(233514, 4968)_(233514, 4968)
11/09/2020 01:28:28 PM Note: NumExpr detected 32 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
11/09/2020 01:28:28 PM NumExpr defaulting to 8 threads.
11/09/2020 01:28:30 PM Note: NumExpr detected 32 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
11/09/2020 01:28:30 PM NumExpr defaulting to 8 threads.
