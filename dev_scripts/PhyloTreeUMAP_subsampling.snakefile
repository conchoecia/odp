"""
The point of this script is to subsample the availble rbh files based on taxonomy.

TODO: - We should also include the sentinel value as a wildcard to enable users to try other sentinel values.
"""

from PhyloTreeUMAP import (algcomboix_file_to_dict,
                           construct_coo_matrix_from_sampledf,
                           filter_sample_df_by_clades,
                           rbh_to_samplename,
                           rbh_directory_to_distance_matrix,
                           rbh_to_distance_gbgz,
                           ALGrbh_to_algcomboix,
                           plot_umap_from_files,
                           plot_umap_pdf,
                           plot_umap_phylogeny_pdf,
                           taxids_to_analyses,
                           taxids_of_interest_to_analyses,
                           topoumap_genmatrix,
                           mgt_mlt_umap,
                           mgt_mlt_plot_HTML,
                           sampleToRbhFileDict_to_sample_matrix,
                           odog_pairwise_distance_matrix,
                           odog_iter_pairwise_distance_matrix, # these are new functions with the pairwise option
                           plot_precomputed_umap,
                           plot_umap_from_files_just_df)

from phylotreeumap_subsample import (generate_subsample_priorities,
                                     subsample_phylogenetically,
                                     make_subsampling_summary_table,
                                     make_subsampling_report_breadcrumbs,
                                     make_subsampling_report_tree)

# get the path of this script, so we know where to look for the plotdfs file
# This block imports fasta-parser as fasta
import os
import sys
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))

import pandas as pd
from scipy.sparse import coo_matrix, lil_matrix, save_npz, load_npz

# Ignore all ResourceWarnings
import warnings
warnings.filterwarnings("ignore", category=ResourceWarning)

configfile: "config.yaml"

# CHECKS
# First verify whether the ALGname has been defined.
if not "ALGname" in config:
    print("ALGname not in config file, assuming \"BCnSSimakov2022\"")
    config["ALGname"] = "BCnSSimakov2022"
# check that rbh_directory exists in the config file
if not "rbh_directory" in config:
    raise ValueError("rbh_directory not in config file")
# make sure that the user has provided the rbh file
if not "ALG_rbh_file" in config:
    raise ValueError("ALG_rbh_file not in config file")

# check that the ALG rbh file exists
if not os.path.exists(config["rbh_directory"]):
    raise ValueError(f"rbh_directory {config['rbh_directory']} does not exist")

# check that there are some NCBI taxids to plot in the config file
if "taxids" in config:
    config["taxids"] = taxids_to_analyses(config["taxids"])
else:
    config["taxids"] = taxids_of_interest_to_analyses()

config["rbh_directory"] = os.path.abspath(config["rbh_directory"])

if "umap" not in config:
    config["umap"] = {}

# It is possible that this could be called many times. Like once per rule run.
# For this reason this needs to be removed s.t. a list of files is created once.
config["rbh_files"] = list(sorted([os.path.join(config["rbh_directory"], f)
                           for f in os.listdir(config["rbh_directory"])
                           if f.endswith('.rbh')], reverse = True))
print("ACCESSING THE RBH FILES")
config["sample_to_rbh_file"] = {rbh_to_samplename(x, config["ALGname"]): x
                                for x in config["rbh_files"]}
print("We identified {} RBH files.".format(len(config["sample_to_rbh_file"])))

# Results_base_directory is the prefix to which everything will be saved
results_base_directory = "GTUMAP"
if results_base_directory.endswith("/"):
    results_base_directory = results_base_directory[:-1]

# Use these parameters for the full space exploration
# These are the parameters in the 2024 biorxiv supplementary figure S7
odog_n    = [20, 35, 50, 75, 100, 150, 250]
odog_m    = [0.0, 0.1, 0.2, 0.5, 0.75, 0.9, 1.0]
odog_n    = [50, 150]
odog_m    = [0.75, 1.0]

# MGT_sentinel_size is now a list of integers representing sentinel values for missing data
# Common values: 0 (small), 999999999999 (large), or any other integer

if "MGT_sentinel_size" not in config:
    config["MGT_sentinel_size"] = [999_999_999_999]
# check that config["MGT_sentinel_size"] is a list of integers
if not isinstance(config["MGT_sentinel_size"], list):
    # raise an error
    raise ValueError("MGT_sentinel_size must be a list of integers")
# check that everything in the list is an integer
for x in config["MGT_sentinel_size"]:
    if not isinstance(x, int):
        raise ValueError("MGT_sentinel_size must be a list of integers")

if "MGT_sentinel_size" in config:
    MGT_sentinel_size = config["MGT_sentinel_size"]
if "odog_n" in config:
    odog_n = config["odog_n"]
if "odog_m" in config:
    odog_m = config["odog_m"]

# get the subsampling dictionary to figure out what ranks to do this for
if "subsample_options" not in config:
    config["subsample_options"] = {}

if "smallest_level" not in config["subsample_options"]:
    config["subsample_options"]["smallest_level"] = "genus"
if "largest_level" not in config["subsample_options"]:
    config["subsample_options"]["largest_level"] = "phylum"
if "max_samples_per_rank" not in config["subsample_options"]:
    config["subsample_options"]["max_samples_per_rank"] = 10

# set the subsampling options
subsample_dict = {x[0]: x for x in generate_subsample_priorities(
      smallest_level = config["subsample_options"]["smallest_level"],
      largest_level  = config["subsample_options"]["largest_level"]
    )}

wildcard_constraints:
    rank="[A-Za-z]+",
    sizeNaN= r"[+-]?\d+",   # allows -42, +7, 0, 123
    n="[0-9]+",

rule all:
    input:
        #    ┓
        # ┏┓┏┫┏┓┏┓ - One-Dot-One-Genome plots - SUBSAMPLING
        # ┗┛┗┻┗┛┗┫   Each dot represents a single genome, and the data vector is the distance pairs
        #        ┛
        # odog old inefficient implementation that uses all locus pairs
        expand(results_base_directory + "/coo_matrices/subsample_{rank}.coo.npz",
                rank = list(subsample_dict.keys())),
        expand(results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
                n = odog_n,
                m = odog_m,
                rank = list(subsample_dict.keys()),
                sizeNaN = MGT_sentinel_size),
        expand(results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html",
                n = odog_n,
                m = odog_m,
                rank = list(subsample_dict.keys()),
                sizeNaN = MGT_sentinel_size),
        expand(results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.pdf",
                n = odog_n,
                m = odog_m,
                rank = list(subsample_dict.keys()),
                sizeNaN = MGT_sentinel_size),
        expand(results_base_directory + "/subsample_umaps/subsample_{rank}.missing_{sizeNaN}.paramsweep.pdf",
                rank = list(subsample_dict.keys()),
                sizeNaN = MGT_sentinel_size),
        expand(results_base_directory + "/subsample_umaps/allranks.neighbors_{n}.mind_{m}.missing_{sizeNaN}.phyloresample.pdf",
                n = odog_n,
                m = odog_m,
                sizeNaN = MGT_sentinel_size),
        # newer implementation
        #results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
        #results_base_directory + "/reduced/pcs.tsv", # this is from feature selection
        #results_base_directory + "/reduced/pca.explained_variance.tsv", # this is also from feature selection
        #expand(results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{metric}.df",
        #        n = odog_n,
        #        m = odog_m,
        #        metric = "euclidean"),
        #expand(results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{metric}.bokeh.html",
        #        n = odog_n,
        #        m = odog_m,
        #        metric = "euclidean"),
        #expand(results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{metric}.pdf",
        #        n = odog_n,
        #        m = odog_m,
        #        metric = "euclidean"),
        #expand(results_base_directory + "/reduced/umap/allsamples.missing_{metric}.paramsweep.pdf",
        #        n = odog_n,
        #        m = odog_m,
        #        metric = "euclidean"),
        # odog new implementation
        #results_base_directory + "/allsamples.sentinel.csr.npz",
        #results_base_directory + "/allsamples.rows.tsv",
        #results_base_directory + "/allsamples.cols.tsv",
        #expand(results_base_directory + "/all_samples_umaps/allsamples.umap.neighbors_{n}.mind_{m}.euclidean.df",
        #        n = odog_n,
        #        m = odog_m,
        #        metric = "euclidean"),
        #expand(results_base_directory + "/allsamples/allsamples.missing_{sizeNaN}.paramsweep.pdf",
        #        sizeNaN = MGT_sentinel_size)
        #    ┓     ┏┓┓ ┏┓┳┓┏┓┏┓
        # ┏┓┏┫┏┓┏┓ ┃ ┃ ┣┫┃┃┣ ┗┓ - One-Dot-One-Genome plots FOR SPECIFIC CLADES
        # ┗┛┗┻┗┛┗┫ ┗┛┗┛┛┗┻┛┗┛┗┛   Each dot represents a single genome, and the data vector is the distance pairs
        #        ┛
        #expand(results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
        #       n       = codog_n,
        #       m       = codog_m,
        #       sizeNaN = codog_size,
        #       taxanalysis = config["taxids"]),
        #expand(results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html",
        #       n       = codog_n,
        #       m       = codog_m,
        #       sizeNaN = codog_size,
        #       taxanalysis = config["taxids"]),
        #expand(results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.phylogeny.pdf",
        #       n       = codog_n,
        #       m       = codog_m,
        #       sizeNaN = codog_size,
        #       taxanalysis = config["taxids"]),
        #expand(results_base_directory + "/ODOG/clades/{taxanalysis}.missing_{sizeNaN}.paramsweep.pdf",
        #       sizeNaN = codog_size,
        #       taxanalysis = config["taxids"]),
        #expand(results_base_directory + "/ODOG/clades/{taxanalysis}.pairwise_distance.missing_{sizeNaN}.tsv",
        #       sizeNaN = codog_size,
        #       taxanalysis = config["taxids"])

def generic_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 2000,
                   2: 4000,
                   3: 8000,
                   4: 16000}
    return attemptdict[attempt]

def pdf_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 4000,
                   2: 8000,
                   3: 16000,
                   4: 32000,
                  }
    return attemptdict.get(attempt, attemptdict[max(attemptdict.keys())])

def odogPairwise_get_mem_mb(wildcards, attempt):
    """Estimate RAM for pairwise matrix generation."""
    attemptdict = {1: 16000,
                   2: 32000,
                   3: 64000,
                   4: 128000}
    return attemptdict[attempt]

def odogPairwise_get_runtime(wildcards, attempt):
    """Estimate runtime for pairwise matrix generation."""
    attemptdict = {1: 60,
                   1: 120,
                   2: 240,
                   4: 480}
    return attemptdict[attempt]

rule odog_paircache_manifest:
    input:
        sampletsv = results_base_directory + "/sampledf.tsv"
    output:
        manifest = results_base_directory + "/pairwise_cache/manifest.json"
    threads: 8
    resources:
        mem_mb  = odogPairwise_get_mem_mb,
        runtime = odogPairwise_get_runtime,
        highio  = 20
    run:
        import os, re, json, tempfile
        from pathlib import Path
        import pandas as pd
        import numpy as np

        cache_dir = os.path.abspath(os.path.dirname(output.manifest))  # <- absolute
        Path(cache_dir).mkdir(parents=True, exist_ok=True)

        cdf = pd.read_csv(input.sampletsv, sep="\t", index_col=0)
        col = "dis_filepath_abs" if "dis_filepath_abs" in cdf.columns else "dis_filepath"
        if col not in cdf.columns:
            raise ValueError("sampledf.tsv must have 'dis_filepath' (or 'dis_filepath_abs').")

        if "sample" in cdf.columns:
            samples = list(cdf["sample"])
            paths   = cdf.set_index("sample")[col].to_dict()
        else:
            samples = list(cdf.index)
            paths   = cdf[col].to_dict()

        def _safe_name(s: str) -> str:
            return re.sub(r"[^A-Za-z0-9._+-]", "_", s)

        distance_dtype = str(config.get("pairwise", {}).get("distance_dtype", "float32"))

        def _pairs_to_uint64_hash(rbh1: pd.Series, rbh2: pd.Series):
            h1 = pd.util.hash_pandas_object(rbh1, index=False).values
            h2 = pd.util.hash_pandas_object(rbh2, index=False).values
            return (h1 ^ h2).astype(np.uint64)

        npz_paths = {}
        for s in samples:
            src = paths[s]
            if not os.path.exists(src):
                raise IOError(f"Missing input file: {src}")

            dst = os.path.abspath(os.path.join(cache_dir, f"{_safe_name(s)}.npz"))  # <- absolute

            # salvage old stray file if present
            stray = dst.replace(".npz", ".npz.tmp.npz")
            if (not os.path.exists(dst)) and os.path.exists(stray):
                os.replace(stray, dst)

            need_build = (not os.path.exists(dst)) or (os.path.getmtime(src) > os.path.getmtime(dst))
            if need_build:
                df = pd.read_csv(src, sep="\t", compression="gzip",
                                 usecols=["rbh1","rbh2","distance"])
                pair_id = _pairs_to_uint64_hash(df["rbh1"], df["rbh2"])
                dist    = df["distance"].astype(distance_dtype).to_numpy()
                ordix   = np.argsort(pair_id, kind="mergesort")
                pair_id = pair_id[ordix]
                dist    = dist[ordix]

                # write to a real .npz temp file (avoid .npz.tmp.npz)
                with tempfile.NamedTemporaryFile(dir=cache_dir,
                                                 prefix=_safe_name(s) + ".",
                                                 suffix=".npz",
                                                 delete=False) as tf:
                    tmp_path = tf.name
                    np.savez_compressed(tf, pair_id=pair_id, dist=dist)
                    tf.flush()
                    try:
                        os.fsync(tf.fileno())
                    except Exception:
                        pass
                os.replace(tmp_path, dst)

            npz_paths[s] = dst

        mani = {
            "samples": samples,
            "npz_paths": npz_paths,            # now absolute values
            "metric": config.get("pairwise", {}).get("metric", "mad"),
            "root_dir": os.path.abspath(os.getcwd())  # helpful for any residual relative paths
        }

        with open(output.manifest, "w") as fh:
            json.dump(mani, fh, indent=2)

def odogFeature_get_runtime(wildcards, attempt):
    """Estimate runtime for pairwise matrix generation."""
    attemptdict = {1: 120,
                   2: 240,
                   3: 480,
                   4: 960}
    return attemptdict[attempt]

def odogFeature_get_mem_mb(wildcards, attempt):
    """Estimate RAM for pairwise matrix generation."""
    attemptdict = {1: 24000,
                   2: 80000,
                   3: 160000,
                   4: 320000}
    return attemptdict[attempt]

rule odog_feature_candidates:
    """
    This picks N candidate features, then picks the top K from those candidates via SVD.

    2.5M candidates takes 1965GB VIRT 800M RAM RES 15m to find candidates
       2446 VIRT 1338 RES for the filtering
    """
    input:
        manifest  = results_base_directory + "/pairwise_cache/manifest.json",
        sampletsv = results_base_directory + "/sampledf.tsv"
    output:
        candidates = results_base_directory + "/reduced/candidates.top.npy",
        stats      = results_base_directory + "/reduced/candidates.stats.tsv"
    params:
        K       = config["feature_select"].get("K", 5000),
        mult    = config["feature_select"].get("candidate_multiplier", 8),
        subsamp = config["feature_select"].get("subsample", 0),
        log1p   = config["feature_select"].get("log1p", True),
        rnd     = config["feature_select"].get("random_state", 42),
        clade_level = config["feature_select"].get("clade_level", 5),
        alpha   = config["feature_select"].get("diversity_alpha", 0.5),
        assume_unique_ids = config["feature_select"].get("assume_unique_ids", True)
    threads: 8
    resources:
        mem_mb  = odogFeature_get_mem_mb,
        runtime = odogFeature_get_runtime,
    run:
        import os, json, math, heapq, time
        import numpy as np
        import pandas as pd
        from collections import Counter

        os.makedirs(os.path.dirname(output.candidates), exist_ok=True)

        # --------------------------
        # Load manifest + sampledf
        # --------------------------
        with open(input.manifest) as fh:
            mani   = json.load(fh)
        samples  = list(mani["samples"])
        npzmap   = mani["npz_paths"]

        sampledf = pd.read_csv(input.sampletsv, sep="\t")
        print(sampledf["sample"].head())

        # quick sanity check
        missing = set(samples) - set(sampledf["sample"])
        if missing:
            raise ValueError(f"Samples in manifest but not in sampledf: {list(missing)[:5]} (total {len(missing)})")

        # --- get clade label safely ---
        def get_clade_label(row, level=params.clade_level):
            lineage = str(row["taxid_list_str"]).split(";")
            if len(lineage) < level:
                return lineage[-1]
            idx = level - 1
            label = lineage[idx]
            # If Eumetazoa, drop one deeper
            if label == "6072" and idx + 1 < len(lineage):
                idx += 1
                label = lineage[idx]
            # If Bilateria, drop two deeper
            if label == "33213" and idx + 2 < len(lineage):
                idx += 2
                label = lineage[idx]
            return label

        sampledf["clade_label"] = sampledf.apply(get_clade_label, axis=1)

        # Map sample -> clade & weight
        clade_sizes = sampledf.groupby("clade_label").size().to_dict()
        print(clade_sizes)
        sample2clade = dict(zip(sampledf["sample"], sampledf["clade_label"]))
        sample2weight = {
            s: 1.0 / clade_sizes[c]
            for s, c in zip(sampledf["sample"], sampledf["clade_label"])
        }
        # print the weights for the first 10 samples
        print("Sample weights (first 10):")
        for s in samples[:10]:
            print(f"  {s}: clade={sample2clade[s]}, weight={sample2weight[s]:.4f}")

        # optional subsample of samples
        if params.subsamp and params.subsamp > 0 and params.subsamp < len(samples):
            rnd = np.random.default_rng(params.rnd)
            idx = rnd.choice(len(samples), size=params.subsamp, replace=False)
            subsamples = [samples[i] for i in idx]
        else:
            subsamples = samples

        def _maybe_log(stage, done, total, start_t, last_t, extras=""):
            now = time.time()
            if (done == 1) or (done == total) or (done % 100 == 0) or (now - last_t >= 10.0):
                elapsed = now - start_t
                rate = (done / elapsed) if elapsed > 0 else 0.0
                eta = ((total - done) / rate) if rate > 0 else float("inf")
                pct = 100.0 * done / max(1, total)
                print(f"[{stage}] {done:,}/{total:,} ({pct:5.1f}%) | "
                      f"elapsed {elapsed:7.1f}s | rate {rate:6.1f}/s | ETA {eta:7.1f}s{extras}",
                      flush=True)
                return now
            return last_t

        # --------------------------
        # Stage A: weighted counts (parallelized)
        # --------------------------
        from concurrent.futures import ThreadPoolExecutor, as_completed
        import heapq

        target = params.K * params.mult
        counts = {}
        heap   = []   # (weighted_count, pid) min-heap, bounded implicitly by counts size
        startA = time.time()

        print(f"[CAND-A] Starting Stage A on {len(subsamples):,} samples | "
              f"K={params.K:,}, multiplier={params.mult}, target pool≈{target:,}",
              flush=True)

        def load_sample(sample):
            """Worker: read npz for one sample and return prepped data for the aggregator."""
            w = sample2weight[sample]
            z = np.load(npzmap[sample], mmap_mode="r")
            ids = z["pair_id"]  # sorted np.uint64
            if params.assume_unique_ids:
                # Fast uniqueness check since ids are sorted
                if ids.size > 1 and np.any(ids[1:] == ids[:-1]):
                    raise ValueError(
                        f"Non-unique pair_ids detected in sample {sample} with assume_unique_ids=True"
                    )
                return sample, w, ids, None  # aggregator will treat each pid once
            else:
                uniq, cnt = np.unique(ids, return_counts=True)
                return sample, w, uniq, cnt   # aggregator consumes (pid, count)

        def _agg_update(sample, w, ids, cnt):
            """Main-thread aggregator: update global counts+heap exactly like your original logic."""
            if cnt is None:
                # assume_unique_ids: each pid contributes w once
                for pid in ids:
                    v = w
                    if pid in counts:
                        counts[pid] += v
                    elif len(counts) < max(10*params.K, target):
                        counts[pid] = v
                        heapq.heappush(heap, (v, int(pid)))
                    else:
                        minc, minpid = heap[0]
                        if v > minc:
                            counts.pop(minpid, None)
                            heapq.heapreplace(heap, (v, int(pid)))
                            counts[int(pid)] = v
            else:
                # safe path with counts per pid
                for pid, c in zip(ids, cnt):
                    v = c * w
                    if pid in counts:
                        counts[pid] += v
                    elif len(counts) < max(10*params.K, target):
                        counts[pid] = v
                        heapq.heappush(heap, (v, int(pid)))
                    else:
                        minc, minpid = heap[0]
                        if v > minc:
                            counts.pop(minpid, None)
                            heapq.heapreplace(heap, (v, int(pid)))
                            counts[int(pid)] = v

        # Launch I/O workers, aggregate as each finishes (no large intermediate storage)
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futures = [ex.submit(load_sample, s) for s in subsamples]
            for i, f in enumerate(as_completed(futures), 1):
                sample, w, ids, cnt = f.result()
                _agg_update(sample, w, ids, cnt)

                # Progress
                if i == 1 or i == len(subsamples) or (i % 100 == 0):
                    elapsed = time.time() - startA
                    rate = i / elapsed if elapsed > 0 else 0.0
                    eta  = (len(subsamples) - i) / rate if rate > 0 else float("inf")
                    print(f"[CAND-A] {i:,}/{len(subsamples):,} ({100*i/len(subsamples):5.1f}%) | "
                          f"elapsed {elapsed:7.1f}s | rate {rate:6.1f}/s | ETA {eta:7.1f}s "
                          f"| cand_map~{len(counts):,}",
                          flush=True)

        print(f"[CAND-A] Done. Building top list from {len(counts):,} tracked pairs...", flush=True)

        # Same finishing step as before (exact semantics preserved)
        top_items = heapq.nlargest(max(target, params.K), [(v, k) for k, v in counts.items()])
        cand_ids  = np.array([int(k) for (v, k) in top_items], dtype=np.uint64)

        print(f"[CAND] Candidate pool size: {cand_ids.size:,}", flush=True)

        # --------------------------
        # Stage B: variance + diversity
        # --------------------------
        id2col = {pid:i for i, pid in enumerate(cand_ids)}
        n   = np.zeros(cand_ids.size, dtype=np.int32)
        mean= np.zeros(cand_ids.size, dtype=np.float64)
        m2  = np.zeros(cand_ids.size, dtype=np.float64)
        clade_counts = {pid: Counter() for pid in cand_ids}

        totalB = len(subsamples)
        startB = lastB  = time.time()

        for i, s in enumerate(subsamples, 1):
            z    = np.load(npzmap[s], mmap_mode="r")
            ids  = z["pair_id"]
            vals = z["dist"].astype(np.float64)
            if params.log1p:
                vals = np.log1p(vals)
            mask = np.in1d(ids, cand_ids, assume_unique=False)
            if mask.any():
                sel_ids  = ids[mask]
                sel_vals = vals[mask]
                clade = sample2clade[s]
                for pid, x in zip(sel_ids, sel_vals):
                    j = id2col[int(pid)]
                    n[j] += 1
                    delta = x - mean[j]
                    mean[j] += delta / n[j]
                    m2[j]   += delta * (x - mean[j])
                    clade_counts[int(pid)][clade] += 1

            extras = f" | seen_features~{int((n>0).sum()):,}" if (i==1 or i==totalB or (i%100==0)) else ""
            lastB = _maybe_log("CAND-B", i, totalB, startB, lastB, extras)

        var = np.zeros_like(mean)
        nz  = n > 1
        var[nz] = m2[nz] / (n[nz] - 1)
        cov = n.astype(np.float64) / max(1, len(subsamples))

        def shannon_entropy(clade_counter):
            total = sum(clade_counter.values())
            return -sum((v/total)*math.log(v/total) for v in clade_counter.values() if v>0)

        entropy = np.array([shannon_entropy(clade_counts[pid]) for pid in cand_ids])
        n_clades = np.array([len(clade_counts[pid]) for pid in cand_ids])

        score = var * cov * (1.0 + params.alpha * entropy)

        # --------------------------
        # Final DataFrame + filtering
        # --------------------------
        df = pd.DataFrame({
            "pair_id": cand_ids.astype(np.uint64),
            "count": n.astype(int),
            "coverage": cov,
            "mean": mean,
            "var": var,
            "entropy": entropy,
            "n_clades": n_clades,
            "score": score
        })

        # Require ≥2 clades
        df = df[df["n_clades"] >= 2]

        keep = max(params.K * params.mult, params.K)
        print(f"[CAND] Selecting top {keep:,} candidates (after >=2 clade filter)...", flush=True)

        df_top = df.sort_values("score", ascending=False).head(keep).reset_index(drop=True)
        df_top.to_csv(output.stats, sep="\t", index=False)
        np.save(output.candidates, df_top["pair_id"].to_numpy(dtype=np.uint64))
        print(f"[CAND] Wrote: {output.stats} and {output.candidates}", flush=True)

rule odog_build_candidate_matrix:
    input:
        manifest   = results_base_directory + "/pairwise_cache/manifest.json",
        candidates = results_base_directory + "/reduced/candidates.top.npy"
    output:
        X_csr = results_base_directory + "/reduced/X_candidates.csr.npz",
        cols  = results_base_directory + "/reduced/candidates.columns.tsv",
        rows  = results_base_directory + "/reduced/samples.tsv"
    params:
        log1p = config["feature_select"].get("log1p", True)
    threads: 8
    resources:
        mem_mb  = odogPairwise_get_mem_mb,
        runtime = odogPairwise_get_runtime,
        highio  = 40
    run:
        import os, json
        import numpy as np
        import pandas as pd
        from scipy.sparse import csr_matrix, save_npz

        os.makedirs(os.path.dirname(output.X_csr), exist_ok=True)

        with open(input.manifest) as fh:
            mani   = json.load(fh)
        samples  = list(mani["samples"])
        npzmap   = mani["npz_paths"]

        cand_ids = np.load(input.candidates).astype(np.uint64)
        id2col   = {int(pid): i for i, pid in enumerate(cand_ids)}

        data = []
        rows = []
        cols = []

        for i, s in enumerate(samples):
            z    = np.load(npzmap[s], mmap_mode="r")
            ids  = z["pair_id"]
            vals = z["dist"].astype(np.float32)
            if params.log1p:
                # compress dynamic range; missing remain 0
                vals = np.log1p(vals, dtype=np.float32)

            # filter to candidate ids
            mask = np.in1d(ids, cand_ids, assume_unique=False)
            if not mask.any():
                continue
            sel_ids  = ids[mask]
            sel_vals = vals[mask]

            # map ids -> cols
            # pair_id unique per sample; safe to iterate
            for pid, x in zip(sel_ids, sel_vals):
                j = id2col[int(pid)]
                rows.append(i)
                cols.append(j)
                data.append(float(x))

        N = len(samples)
        C = cand_ids.size
        X = csr_matrix((np.array(data, dtype=np.float32),
                        (np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32))),
                       shape=(N, C))
        save_npz(output.X_csr, X)
        pd.Series(cand_ids, name="pair_id").to_csv(output.cols, sep="\t", index=False)
        pd.Series(samples,  name="sample").to_csv(output.rows, sep="\t", index=False)


rule odog_select_topK_via_svd:
    input:
        X_csr = results_base_directory + "/reduced/X_candidates.csr.npz",
        cols  = results_base_directory + "/reduced/candidates.columns.tsv"
    output:
        sel_ids = results_base_directory + "/reduced/selected.topK.npy",
        svdinfo = results_base_directory + "/reduced/candidate_svd.info.tsv"
    params:
        K        = config["feature_select"].get("K", 5000),
        ncomp    = config["feature_select"].get("svd_components", 50),
        niter    = config["feature_select"].get("svd_iter", 7),
        rnd      = config["feature_select"].get("random_state", 42),
        blas_thr = config["feature_select"].get("threads_blas", 8)
    threads: 8
    resources:
        mem_mb  = odogPairwise_get_mem_mb,
        runtime = odogPairwise_get_runtime,
        highio  = 20
    run:
        import os
        for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
                os.environ[v] = str(threads)

        import numpy as np
        import pandas as pd
        from scipy.sparse import load_npz
        from sklearn.decomposition import TruncatedSVD

        X = load_npz(input.X_csr)   # (N x C), sparse
        cols = pd.read_csv(input.cols, sep="\t")["pair_id"].to_numpy(dtype=np.uint64)

        # SVD on candidate matrix (no centering). Works like LSA on sparse data.
        svd = TruncatedSVD(n_components=min(params.ncomp, min(X.shape)-1),
                           n_iter=params.niter, random_state=params.rnd)
        Z = svd.fit_transform(X)   # not used here; we need component loadings

        # Feature importance via sum of squared loadings weighted by explained variance
        # components_: (n_components x n_features), each row normalized to unit norm roughly
        comp = svd.components_               # shape (p, C)
        evr  = svd.explained_variance_ratio_ # shape (p,)
        weights = (comp**2) * evr[:, None]   # (p x C)
        importance = weights.sum(axis=0)     # (C,)

        order = np.argsort(importance)[::-1]
        keep  = order[:params.K]
        sel_ids = cols[keep]

        np.save(output.sel_ids, sel_ids)
        pd.DataFrame({
            "pair_id": cols,
            "importance": importance
        }).sort_values("importance", ascending=False).to_csv(output.svdinfo, sep="\t", index=False)

# --- change this one to CONSUME samples.tsv instead of producing it ---
rule odog_build_topK_matrix:
    input:
        manifest = results_base_directory + "/pairwise_cache/manifest.json",
        sel_ids  = results_base_directory + "/reduced/selected.topK.npy",
        rows     = results_base_directory + "/reduced/samples.tsv"   # <— now an input
    output:
        X    = results_base_directory + "/reduced/X_topK.csr.npz",
        cols = results_base_directory + "/reduced/topK.columns.tsv"
    params:
        log1p = config["feature_select"].get("log1p", True)
    threads: 8
    resources:
        mem_mb  = odogPairwise_get_mem_mb,
        runtime = odogPairwise_get_runtime,
        highio  = 40
    run:
        import os, json
        import numpy as np
        import pandas as pd
        from scipy.sparse import csr_matrix, save_npz

        with open(input.manifest) as fh:
            mani   = json.load(fh)
        samples  = pd.read_csv(input.rows, sep="\t")["sample"].tolist()
        npzmap   = mani["npz_paths"]

        sel = np.load(input.sel_ids).astype(np.uint64)
        id2col = {int(pid): i for i, pid in enumerate(sel)}

        data, rows_idx, cols_idx = [], [], []
        for i, s in enumerate(samples):
            z    = np.load(npzmap[s], mmap_mode="r")
            ids  = z["pair_id"]
            vals = z["dist"].astype(np.float32)
            if params.log1p:
                vals = np.log1p(vals, dtype=np.float32)
            mask = np.in1d(ids, sel, assume_unique=False)
            if not mask.any():
                continue
            sel_ids  = ids[mask]
            sel_vals = vals[mask]
            for pid, x in zip(sel_ids, sel_vals):
                j = id2col[int(pid)]
                rows_idx.append(i); cols_idx.append(j); data.append(float(x))

        N, K = len(samples), sel.size
        X = csr_matrix((np.array(data, dtype=np.float32),
                        (np.array(rows_idx, dtype=np.int32), np.array(cols_idx, dtype=np.int32))),
                       shape=(N, K))
        os.makedirs(os.path.dirname(output.X), exist_ok=True)
        save_npz(output.X, X)
        pd.Series(sel, name="pair_id").to_csv(output.cols, sep="\t", index=False)

rule odog_pca_topK:
    input:
        X    = results_base_directory + "/reduced/X_topK.csr.npz",
        rows = results_base_directory + "/reduced/samples.tsv",
        cols = results_base_directory + "/reduced/topK.columns.tsv"
    output:
        pcs  = results_base_directory + "/reduced/pcs.tsv",
        evr  = results_base_directory + "/reduced/pca.explained_variance.tsv"
    params:
        ncomp    = config["feature_select"].get("svd_components", 50),
        niter    = config["feature_select"].get("svd_iter", 7),
        rnd      = config["feature_select"].get("random_state", 42),
        blas_thr = config["feature_select"].get("threads_blas", 8)
    threads: 8
    resources:
        mem_mb  = odogPairwise_get_mem_mb,
        runtime = odogPairwise_get_runtime,
        highio  = 10
    run:
        import os
        for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
                os.environ[v] = str(threads)

        import numpy as np
        import pandas as pd
        from scipy.sparse import load_npz
        from sklearn.decomposition import TruncatedSVD

        X = load_npz(input.X)  # (N x K), sparse
        samples = pd.read_csv(input.rows, sep="\t")["sample"].tolist()

        svd = TruncatedSVD(n_components=min(params.ncomp, min(X.shape)-1),
                           n_iter=params.niter, random_state=params.rnd)
        Z = svd.fit_transform(X)  # (N x ncomp)

        # Write PCs with sample names
        pc_cols = [f"PC{i+1}" for i in range(Z.shape[1])]
        df = pd.DataFrame(Z, columns=pc_cols, index=samples)
        df.index.name = "sample"
        df.to_csv(output.pcs, sep="\t")

        # Write explained variance
        pd.DataFrame({
            "component": pc_cols,
            "explained_variance_ratio": svd.explained_variance_ratio_,
            "singular_values": svd.explained_variance_ ** 0.5
        }).to_csv(output.evr, sep="\t", index=False)

rule odog_umap_from_pcs:
    input:
        pcs      = results_base_directory + "/reduced/pcs.tsv",
        sampletsv= results_base_directory + "/sampledf.tsv"
    output:
        df   = results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{metric}.df"
    params:
        metric       = config["umap"].get("metric", "euclidean"),
        densmap      = config["umap"].get("densmap", False),
        random_state = config["umap"].get("random_state", 42),
        init         = config["umap"].get("init", "random"),
        low_memory   = config["umap"].get("low_memory", False)
    threads: 4
    resources:
        mem_mb  = odogPairwise_get_mem_mb,
        runtime = odogPairwise_get_runtime
    run:
        import os, time, re
        import numpy as np
        import pandas as pd
        import umap

        # --- load PCs robustly ---
        DF = pd.read_csv(input.pcs, sep="\t", index_col=0, engine="python")
        # trim stray whitespace in column names
        DF.columns = [str(c).strip() for c in DF.columns]
        # keep only columns named like PC<number>, sorted by number
        pc_cols = [c for c in DF.columns if re.fullmatch(r"PC\d+", str(c))]
        if not pc_cols:
            # fallback: any column starting with 'PC'
            pc_cols = [c for c in DF.columns if str(c).startswith("PC")]
        # sort by numeric suffix if possible
        def _pc_key(c):
            s = str(c)[2:]
            return int(s) if s.isdigit() else 10**9
        pc_cols = sorted(pc_cols, key=_pc_key)
        X = DF[pc_cols].to_numpy(dtype=np.float32, copy=False)
        samples = DF.index.tolist()

        n_neighbors = int(wildcards.n)
        min_dist    = float(wildcards.m)
        metric      = str(wildcards.metric)

        print(DF)

        if X.shape[0] <= n_neighbors:
            print(f"[UMAP-PC] n_samples={X.shape[0]} <= n_neighbors={n_neighbors}. Writing empty outputs.")
            open(output.df, "w").close()
            open(output.html, "w").close()
            raise SystemExit(0)

        print(f"[UMAP-PC] X shape={X.shape}, PCs used={len(pc_cols)}, "
              f"n_neighbors={n_neighbors}, min_dist={min_dist}, metric={metric}, densmap={params.densmap}", flush=True)

        reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            densmap=bool(params.densmap),
            init=params.init,
            random_state=int(params.random_state),
            low_memory=bool(params.low_memory),
        )

        # fit like your existing function (so mapper has graph for bokeh helper)
        t0 = time.time()
        mapper = reducer.fit(X)
        t1 = time.time()
        print(f"[UMAP-PC] fit time: {t1 - t0:.1f}s", flush=True)

        # write df
        from PhyloTreeUMAP import umap_mapper_to_df, umap_mapper_to_bokeh
        umap_df = umap_mapper_to_df(mapper, pd.read_csv(input.sampletsv, sep="\t", index_col=0))
        os.makedirs(os.path.dirname(output.df), exist_ok=True)
        umap_df.to_csv(output.df, sep="\t", index=True)

rule odog_umap_plot_html:
    input:
        df = results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{metric}.df"
    output:
        html = results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{metric}.bokeh.html",
    threads: 1
    retries: 4
    resources:
        mem_mb  = pdf_get_mem_mb,
        runtime = 10
    run:
        # write bokeh
        n_neighbors = int(wildcards.n)
        min_dist    = float(wildcards.m)
        metric      = str(wildcards.metric)
        mgt_mlt_plot_HTML(input.df, output.html,
                          plot_title=f"UMAP (PCs) n={n_neighbors}, min_dist={min_dist}, metric={metric}",
                          analysis_type = "MGT")

rule odog_umap_plot_pdf:
    input:
        df = results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{metric}.df"
    output:
        pdf  = results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{metric}.pdf"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        runtime = 10
    run:
        title = f"UMAP (PCs) n={wildcards.n}, min_dist={wildcards.m}, metric={wildcards.metric}"
        plot_umap_pdf(input.df, output.pdf,
                      title, color_by_clade = True)

rule odogSweep:
    """
    Makes a pdf of the parameter sweep of the all-species UMAP plots.
    """
    input:
        dfs =  expand(results_base_directory + "/reduced/umap/pcs.umap.neighbors_{n}.mind_{m}.{{metric}}.df",
                     n = odog_n,
                     m = odog_m),
        plotdfs = os.path.join(snakefile_path, "PhyloTreeUMAP_plotdfs.py")
    output:
        pdf    = results_base_directory + "/reduced/umap/allsamples.missing_{metric}.paramsweep.pdf",
        html   = results_base_directory + "/reduced/umap/allsamples.missing_{metric}.paramsweep.html",
    params:
        prefix = results_base_directory + "/reduced/umap/allsamples.missing_{metric}.paramsweep"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        highio = 1,
        runtime = 5
    shell:
        """
        python {input.plotdfs} -f "{input.dfs}" -p {params.prefix} --pdf --html
        """

#  ------------ Feature selection stops here --------------

# ┏┓    ┓        ┓         ┓
# ┗┓┏┓┏┓┃┏┏┓┏┳┓┏┓┃┏┏┓  ┏┓┓┏┃┏┓┏
# ┗┛┛┗┗┻┛┗┗ ┛┗┗┗┻┛┗┗   ┛ ┗┻┗┗ ┛
#

def is_valid_gzip(file_path):
    """
    Just check if the gzip file is valid or complete.
    """
    try:
        with gzip.open(file_path, 'rb') as f:
            while f.read(1024 * 1024):  # Read in chunks of 1MB
                pass
        return True
    except (OSError, gzip.BadGzipFile):
        return False

rule samples_and_gzipped:
    """
    This version works as originally planned. It runs once per file that we need to make.
    If you run too many at once it will freeze the IO on the disk system.

    The samples need to be converted to a different format before stiching them together into
      distance matrices in numpy. In an earlier version of this pipeline, I had this rule run once
      per file. This caused a high IO burden, so now I am running this in a for loop. This actually
      runs much faster, and each operation takes less than one second typically.
    """
    input:
        rbh_file = lambda wildcards: config["sample_to_rbh_file"][wildcards.sample]
    output:
        gbgz = results_base_directory + "/distance_matrices/{sample}.gb.gz"
    threads: 1
    params:
        ALGname   = config["ALGname"]
    retries: 3
    resources:
        mem_mb  = generic_get_mem_mb,
        runtime = 10,
        #highio  = 1,
    run:
        rbh_to_distance_gbgz(input.rbh_file, output.gbgz, params.ALGname)

rule sample_df_all:
    input:
        rbh_files = [config["sample_to_rbh_file"][x]
                     for x in config["sample_to_rbh_file"]]
    output:
        sampletsv = results_base_directory + "/sampledf.tsv"
    retries: 6
    threads: 1
    params:
        gbgz_directory = results_base_directory + "/distance_matrices"
    resources:
        mem_mb = 1000,
        highio = 1,
        runtime = 10
    run:
        sampleToRbhFileDict_to_sample_matrix(config["sample_to_rbh_file"],
                                             config["ALGname"],
                                             params.gbgz_directory,
                                             output.sampletsv)

rule combo_to_index:
    input:
        rbhfile = config["ALG_rbh_file"]
    output:
        outfile = results_base_directory + "/combo_to_index.txt"
    threads: 1
    retries: 3
    resources:
        mem_mb = generic_get_mem_mb,
        highio = 1,
        runtime = 5
    run:
        # We will need to calculate the rbh combo to index dictionary.
        # DO NOT bother reading in the existing file. It takes 5x longer
        #  to read in the file than it does to generate it.
        alg_combo_to_ix = ALGrbh_to_algcomboix(input.rbhfile)
        # save the dictionary pairs
        with open(output.outfile, "w") as f:
            for key, value in alg_combo_to_ix.items():
                f.write(f"{key}\t{value}\n")

def coo_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 10000,
                   2: 25000,
                   3: 50000,
                   4: 100000,
                   5: 200000,
                   6: 300000,
                   7: 400000}
    return attemptdict[attempt]

def coo_get_runtime(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 60,
                   2: 120,
                   3: 150,
                   4: 180,
                   5: 210,
                   6: 240,
                   7: 270}
    return attemptdict[attempt]

#    ┓
# ┏┓┏┫┏┓┏┓ - One-Dot-One-Genome plots
# ┗┛┗┻┗┛┗┫   Each dot represents a single genome, and the data vector is the distance pairs
#        ┛

def allsample_coo_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the number of input genomes.
    The script needs about 0.05 GB per genome.
    So 5831 genomes needs around 300 GB of RAM.
    """
    import os

    sampledf_path = results_base_directory + f"/subsample_info/subsample_{wildcards.rank}_sampledf.tsv"

    # Count genomes (lines - header), ignore blank/comment lines
    n = 0
    try:
        with open(sampledf_path, "r", encoding="utf-8", errors="ignore") as f:
            next(f, None)  # header
            for line in f:
                s = line.strip()
                if s and not s.startswith("#"):
                    n += 1
    except FileNotFoundError:
        # Safe fallback if file isn't available at DAG eval time
        n = int(config.get("umap_default_n_genomes", 1000))

    n = max(0, n)

    # Per-attempt tiers (MiB). 0.05 GiB/genome ⇒ 0.05*1024 MiB/genome.
    tiers = {
        1: max(30000, int(0.051 * n * 1024)),  # ≥ 30 GiB
        2: max(35000, int(0.061 * n * 1024)),  # ≥ 35 GiB
        3: max(40000, int(0.071 * n * 1024)),  # ≥ 40 GiB
        4: max(50000, int(0.081 * n * 1024)),  # ≥ 50 GiB
    }
    return tiers.get(attempt, tiers[max(tiers)])  # use last tier for attempt>=5

rule subsampling_df:
    """
    This time, use the subsampling function to make lists of which samples per subsampling to use.
    These lists will be used later to make the smaller coo files.
    """
    input:
        sampletsv = results_base_directory + "/sampledf.tsv"
    output:
        breadcrumbs = results_base_directory + "/subsample_info/subsample_{rank}_breadcrumbs.txt",
        selectedsam = results_base_directory + "/subsample_info/subsample_{rank}_selected_samples.txt",
        summarytabl = results_base_directory + "/subsample_info/subsample_{rank}_summary_table.tsv",
        sampletree  = results_base_directory + "/subsample_info/subsample_{rank}_sample_tree.txt"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 10
    run:
        from pathlib import Path
        df = pd.read_csv(input.sampletsv, sep="\t", index_col=0)
        list_of_ranks = subsample_dict[wildcards.rank]
        if wildcards.rank == "allsamples":
            selected_buckets, flat = subsample_phylogenetically(df, select_all = True)
        else:
            selected_buckets, flat = subsample_phylogenetically(df,
                                                                bucket_priority=list_of_ranks,
                                                                max_per_bucket=10,
                                                                priority=True)
        summary_txt     = make_subsampling_summary_table(selected_buckets)
        breadcrumbs_txt = make_subsampling_report_breadcrumbs(selected_buckets)
        tree_txt        = make_subsampling_report_tree(selected_buckets)

        basedir = os.path.dirname(output.breadcrumbs)
        outdir = Path(basedir)
        outdir.mkdir(parents=True, exist_ok=True)
        Path(output.breadcrumbs).write_text(breadcrumbs_txt, encoding="utf-8")
        Path(output.selectedsam).write_text("\n".join(sorted(flat)) + "\n", encoding="utf-8")
        Path(output.summarytabl).write_text(summary_txt, encoding="utf-8")
        Path(output.sampletree).write_text(tree_txt, encoding="utf-8")

rule sampletsv_of_subsampling:
    """
    Make a sample table just for the subsampled genomes.
    Read in the full sample table, then use the list of selected samples to filter it down.
    """
    input:
        sampletsv = results_base_directory + "/sampledf.tsv",
        selectedsam = results_base_directory + "/subsample_info/subsample_{rank}_selected_samples.txt",
    output:
        outtsv = results_base_directory + "/subsample_info/subsample_{rank}_sampledf.tsv"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 5
    run:
        import pandas as pd
        selected = set()
        with open(input.selectedsam, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if s:
                    selected.add(s)
        df = pd.read_csv(input.sampletsv, sep="\t", index_col=0)
        # select on the "sample" column
        df_sub = df[df["sample"].isin(selected)].copy()
        df_sub = df_sub.reset_index(drop=True)
        df_sub.to_csv(output.outtsv, sep="\t", index=True)

rule odogCooGen_subsample:
    """
    This generates a coo matrix of the distance matrices for all genomes.
      Each row is a genome, and each column is the distance between every locus pair.
      There is no pre-processing performed on the distance matrices. These are the raw
      distances between every locus pair.

    This is for the old, inefficient version of the UMAP plotting that uses all of the distances.

    This takes up a lot of RAM and time. For 3600 genomes, it took:
        CPU Efficiency: 97.32% of 01:57:32 core-walltime
        Job Wall-clock time: 01:57:32
        Memory Utilized: 216.61 GB
        Memory Efficiency: 110.91% of 195.31 GB
    """
    input:
        gbgzfiles = expand(results_base_directory + "/distance_matrices/{sample}.gb.gz",
                           sample = config["sample_to_rbh_file"].keys()),
        combotoindex = results_base_directory + "/combo_to_index.txt",
        sampletsv = results_base_directory + "/subsample_info/subsample_{rank}_sampledf.tsv"
    output:
        coo    = results_base_directory + "/coo_matrices/subsample_{rank}.coo.npz"
    threads: 1
    retries: 3
    shadow: "copy-minimal"
    resources:
        slurm_extra = "--no-requeue",
        tmpdir = "/tmp",
        mem_mb = allsample_coo_get_mem_mb,
        #highio = 10, # I don't think we need this now since we stage to local disk
        runtime = 300
    run:
        import os, shutil, tempfile, time
        from pathlib import Path
        import pandas as pd
        from scipy.sparse import save_npz

        work = None
        try:
            # -------- scratch directory selection ----------
            scratch_base = (
                os.environ.get("SLURM_TMPDIR")
                or os.environ.get("TMPDIR")
                or getattr(resources, "tmpdir", None)
                or "/tmp"
            )
            work = Path(tempfile.mkdtemp(prefix="odogcoo_", dir=scratch_base))
            print(f"[DOG-COO] Staging to {work}", flush=True)

            # -------- copy inputs with progress ----------
            t0 = time.time()
            combo_local    = work / Path(input.combotoindex).name
            sampledf_local = work / Path(input.sampletsv).name
            shutil.copy2(input.combotoindex, combo_local)
            shutil.copy2(input.sampletsv,    sampledf_local)

            # get the gb.gz list from sample TSV
            sampledf  = pd.read_csv(sampledf_local, sep="\t", index_col=0)
            pathcol   = "dis_filepath_abs"
            samplecol = "sample"
            for checking_this_col in (pathcol, samplecol):
                if checking_this_col not in sampledf.columns:
                    raise ValueError(f"Expected column '{checking_this_col}' in sample TSV")

            # check for duplicate sample names - this would break everything
            duplicate_sample_names = sampledf[samplecol][sampledf[samplecol].duplicated()].unique()
            if len(duplicate_sample_names) > 0:
                raise ValueError(f"Duplicate sample names found in sample TSV: {duplicate_sample_names}")

            # ---------- Stage only those files -----------------
            local_gbgz = {}
            N = len(sampledf)
            tcopy = time.time()

            for i, (_, row) in enumerate(sampledf.iterrows(), 1):
                samp = str(row[samplecol])
                src  = Path(str(row[pathcol]))
                if not src.name.endswith(".gb.gz"):
                    raise ValueError(f"Sample '{samp}' path does not end with .gb.gz: {src}")
                if not src.exists():
                    raise IOError(f"Input .gb.gz file does not exist: {src}")

                dst = work / src.name
                shutil.copy2(src, dst)
                local_gbgz[samp] = str(dst)

                if i % max(1, N // 10) == 0 or i == N:
                    elapsed = time.time() - tcopy
                    print(f"[DOG-COO] Copied {i}/{N} files (elapsed {elapsed:.1f}s)", flush=True)

            # -------- run computation using local copies ----------
            print("[DOG-COO] Generating COO...", flush=True)
            alg_combo_to_ix = algcomboix_file_to_dict(str(combo_local))

            coo = construct_coo_matrix_from_sampledf(
                sampledf,
                alg_combo_to_ix,
                print_prefix="[DOG-COO] ",
                gbgz_paths=local_gbgz,          # force local paths
                path_column="dis_filepath_abs", # ignored since list is passed
            )
            Path(output.coo).parent.mkdir(parents=True, exist_ok=True)
            save_npz(output.coo, coo)
            print(f"[DOG-COO] Done. Total stage+compute time {time.time()-t0:.1f}s", flush=True)

        finally:
            # Always remove scratch directory; ignore errors so cleanup never crashes the job
            try:
                if 'work' in locals() and work is not None:
                    if hasattr(work, "exists") and work.exists():
                        shutil.rmtree(work, ignore_errors=True)
            except Exception as e:
                print(f"[DOG-COO] Warning: failed to remove scratch dir {work}: {e}", flush=True)

def odogPlotUMAP_old_inefficient_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the number of input genomes.
    The script needs about 0.05 GB per genome.
    So 5831 genomes needs around 300 GB of RAM.
    """
    import os

    sampledf_path = results_base_directory + f"/subsample_info/subsample_{wildcards.rank}_sampledf.tsv"

    # Count genomes (lines - header), ignore blank/comment lines
    n = 0
    try:
        with open(sampledf_path, "r", encoding="utf-8", errors="ignore") as f:
            next(f, None)  # header
            for line in f:
                s = line.strip()
                if s and not s.startswith("#"):
                    n += 1
    except FileNotFoundError:
        # Safe fallback if file isn't available at DAG eval time
        n = int(config.get("umap_default_n_genomes", 1000))

    n = max(0, n)

    # Per-attempt tiers (MiB). 0.05 GiB/genome ⇒ 0.05*1024 MiB/genome.
    tiers = {
        1: max(30000, int(0.07 * n * 1024)),  # ≥ 30 GiB
        2: max(35000, int(0.08 * n * 1024)),  # ≥ 35 GiB
        3: max(40000, int(0.09 * n * 1024)),  # ≥ 40 GiB
        4: max(50000, int(0.10 * n * 1024)),  # ≥ 50 GiB
    }
    return tiers.get(attempt, tiers[max(tiers)])  # use last tier for attempt>=5

def odogPlotUMAP_old_inefficient_get_runtime(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    sampledf_path = results_base_directory + f"/subsample_info/subsample_{wildcards.rank}_sampledf.tsv"

    # Count genomes (lines - header), ignore blank/comment lines
    n = 0
    try:
        with open(sampledf_path, "r", encoding="utf-8", errors="ignore") as f:
            next(f, None)  # header
            for line in f:
                s = line.strip()
                if s and not s.startswith("#"):
                    n += 1
    except FileNotFoundError:
        # Safe fallback if file isn't available at DAG eval time
        n = int(config.get("umap_default_n_genomes", 1000))

    n = max(10, n)

    # this is the number of minutes per sample. Scale it by threads
    # was 1 min per sample for single-threaded. It took about 9 hours for 5831 samples.
    # Now we have 32 threads, so it should be much faster.
    # Tuning it to 1/4 of a minute per sample with 32 threads.
    attemptdict = {1: int(1.0 * n * 0.33),
                   2: int(1.2 * n * 0.33),
                   3: int(1.4 * n * 0.33),
                   4: int(1.6 * n * 0.33)}
    return attemptdict[attempt]

rule odogPlotUMAP_old_inefficient:
    input:
        sampletsv    = results_base_directory + "/subsample_info/subsample_{rank}_sampledf.tsv",
        coo          = results_base_directory + "/coo_matrices/subsample_{rank}.coo.npz",
        combotoindex = results_base_directory + "/combo_to_index.txt"
    output:
        df   = results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df"
    threads: 32
    retries: 3
    resources:
        mem_mb  = odogPlotUMAP_old_inefficient_get_mem_mb,
        runtime = odogPlotUMAP_old_inefficient_get_runtime,
        bigUMAPSlots = 1,
    benchmark:
        results_base_directory + "/benchmarks/odogPlotUMAP_old_inefficient/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.tsv"
    run:
        import os, shutil, tempfile, time
        from pathlib import Path
        import pandas as pd

        os.makedirs(os.path.dirname(output.df), exist_ok=True)

        # check if the input sampletsv has enough samples, else write the header of sampleTSV plus the columns "UMAP1" and "UMAP2"
        df = pd.read_csv(input.sampletsv, sep="\t", index_col=0)
        if df.shape[0] <= int(wildcards.n):
            print(f"[ODOG-UMAP] n_samples={df.shape[0]} <= n_neighbors={wildcards.n}. Writing empty outputs.")
            # write header only
            with open(output.df, "w", encoding="utf-8") as fh:
                column_string = "\t".join(df.columns)
                print("Column string is: ", column_string)
                fh.write(f"{column_string}\tUMAP1\tUMAP2\n")
            return # exit cleanly. SystemExit(0) would be a failure in Snakemake

        # There are enough samples, so the UMAP should always produce something.
        work = None
        try:
            # -------- scratch directory selection ----------
            scratch_base = (
                os.environ.get("SLURM_TMPDIR")
                or os.environ.get("TMPDIR")
                or getattr(resources, "tmpdir", None)
                or "/tmp"
            )
            work = Path(tempfile.mkdtemp(prefix="odogumap_", dir=scratch_base))
            print(f"[ODOG-UMAP] Staging inputs to {work}", flush=True)

            # -------- copy inputs ----------
            sample_local    = work / "sampledf.tsv"
            combo_local     = work / "combo_to_index.txt"
            coo_local       = work / "allsamples.coo.npz"

            shutil.copy2(input.sampletsv,    sample_local)
            shutil.copy2(input.combotoindex, combo_local)
            shutil.copy2(input.coo,          coo_local)

            # -------- run DF-only UMAP on local copies ----------
            print(f"[ODOG-UMAP] Starting the UMAP calculation.", flush=True)
            # Convert sizeNaN wildcard to integer sentinel value
            sentinel_value = int(wildcards.sizeNaN)
            plot_umap_from_files_just_df(
                sampledffile     = str(sample_local),
                ALGcomboixfile   = str(combo_local),
                coofile          = str(coo_local),
                sample           = "allsamples",
                smalllargeNaN    = sentinel_value,  # Now an integer sentinel value
                n_neighbors      = int(wildcards.n),
                min_dist         = float(wildcards.m),
                dfoutfilepath    = output.df,
                missing_value_as = sentinel_value,  # Use the same sentinel value
                threads          = threads,
                print_prefix     = "[ODOG-UMAP] ",
            )
            print(f"[ODOG-UMAP] Wrote {output.df}", flush=True)

        finally:
            # Always remove scratch directory; ignore errors so cleanup never crashes the job
            try:
                if work is not None and work.exists():
                    print(f"[ODOG-UMAP] Removing files from the scratch directory {work}", flush=True)
                    shutil.rmtree(work, ignore_errors=True)
            except Exception as e:
                print(f"[ODOG-UMAP] Warning: failed to remove scratch dir {work}: {e}", flush=True)

rule odog_umap_oldinefficient_plot_html:
    input:
        df   = results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
    output:
        html = results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html",
    threads: 1
    retries: 4
    resources:
        mem_mb  = pdf_get_mem_mb,
        runtime = 10
    run:
        import os
        import pandas as pd
        from pandas.errors import EmptyDataError
        os.makedirs(os.path.dirname(output.html), exist_ok=True)

        n_neighbors = int(wildcards.n)
        min_dist    = float(wildcards.m)
        sizeNaN     = str(wildcards.sizeNaN)
        plot_title = f"UMAP (PCs) n={n_neighbors}, min_dist={min_dist}, missing size={sizeNaN}"

        # Detect "empty df" robustly
        empty = False
        try:
            if not os.path.exists(input.df) or os.path.getsize(input.df) == 0:
                empty = True
            else:
                _df = pd.read_csv(input.df, sep="\t", index_col=0, nrows=1)
                empty = (_df.shape[0] == 0)
        except EmptyDataError:
            empty = True

        if empty:
            # Minimal HTML with your message
            msg = (f"The input DF with n_neighbors={n_neighbors}, "
                   f"min_dist={min_dist}, missing size={sizeNaN} had an empty df")
            with open(output.html, "w", encoding="utf-8") as fh:
                fh.write(f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
                         f"<title>Empty UMAP</title></head><body>"
                         f"<h3>{msg}</h3></body></html>")
        else:
            # Normal plotting path
            mgt_mlt_plot_HTML(input.df, output.html, plot_title=plot_title, analysis_type="MGT", plot_sizing_mode = "scale_both")

rule odog_umap_oldinefficient_plot_pdf:
    input:
        df   = results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
    output:
        pdf  = results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.pdf"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        runtime = 10
    run:
        title = f"UMAP (PCs) n={wildcards.n}, min_dist={wildcards.m}, missing size={wildcards.sizeNaN}"
        plot_umap_pdf(input.df, output.pdf,
                      title, color_by_clade = True)

rule odogSweep_oldinefficient:
    """
    Makes a pdf of the parameter sweep of the all-species UMAP plots.
    """
    input:
        dfs = expand(results_base_directory + "/subsample_umaps/{{rank}}/subsample_{{rank}}.neighbors_{n}.mind_{m}.missing_{{sizeNaN}}.df",
                     n = odog_n,
                     m = odog_m),
        plotdfs = os.path.join(snakefile_path, "PhyloTreeUMAP_plotdfs.py")
    output:
        pdf    = results_base_directory + "/subsample_umaps/subsample_{rank}.missing_{sizeNaN}.paramsweep.pdf",
        html   = results_base_directory + "/subsample_umaps/subsample_{rank}.missing_{sizeNaN}.paramsweep.html",
    params:
        prefix = results_base_directory + "/subsample_umaps/subsample_{rank}.missing_{sizeNaN}.paramsweep"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        highio = 1,
        runtime = 5
    shell:
        """
        python {input.plotdfs} -f "{input.dfs}" -p {params.prefix} --pdf --html
        """

rule odogSweep_subsampling:
    input:
        dfs = expand(results_base_directory + "/subsample_umaps/{rank}/subsample_{rank}.neighbors_{{n}}.mind_{{m}}.missing_{{sizeNaN}}.df",
                        rank = subsample_dict.keys()),
        plotdfs = os.path.join(snakefile_path, "PhyloTreeUMAP_plotdfs.py")
    output:
        pdf    = results_base_directory + "/subsample_umaps/allranks.neighbors_{n}.mind_{m}.missing_{sizeNaN}.phyloresample.pdf",
    params:
        prefix = results_base_directory + "/subsample_umaps/allranks.neighbors_{n}.mind_{m}.missing_{sizeNaN}"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        highio = 1,
        runtime = 5
    shell:
        """
        python {input.plotdfs} --phylolist "{input.dfs}" -p {params.prefix} --pdf
        """

def odogPlotUMAP_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 25000,
                   2: 30000,
                   3: 35000,
                   4: 40000}
    return attemptdict.get(attempt, attemptdict[max(attemptdict.keys())])

def odogPlotUMAP_get_runtime(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {#1: 300,
                   1: 360,
                   2: 420,
                   3: 480,
                   4: 540}
    return attemptdict[attempt]

rule odog_full_sentinel_csr:
    input:
        manifest  = results_base_directory + "/pairwise_cache/manifest.json",
        sampletsv = results_base_directory + "/sampledf.tsv"
    output:
        X = results_base_directory + "/allsamples.sentinel.csr.npz",
        rows = results_base_directory + "/allsamples.rows.tsv",
        cols = results_base_directory + "/allsamples.cols.tsv"
    params:
        log1p    = config["feature_select"].get("log1p", False),  # your choice for the *full* run
        M_raw    = int(config.get("sentinel", 999_999_999_999)),  # 999 billion by default
        use_sqlite_map = True
    threads: 8
    shadow: "copy-minimal"   # <<< ensures inputs copied to fast local disk
    resources:
        mem_mb = coo_get_mem_mb,
        runtime = coo_get_runtime
    run:
        import os, json, sqlite3, math, time, gc, shutil
        import numpy as np
        import pandas as pd
        from scipy.sparse import csr_matrix, save_npz
        from concurrent.futures import ThreadPoolExecutor

        # ---------- load manifest ----------
        with open(input.manifest) as fh:
            mani = json.load(fh)
        samples = list(mani["samples"])
        npzmap  = dict(mani["npz_paths"])  # sample -> ABS path
        root_dir = mani.get("root_dir", None)

        # ---------- stage npz locally (hardlink if possible) ----------
        local_root   = os.environ.get("TMPDIR", os.getcwd())
        local_npzdir = os.path.join(local_root, "npz_local")
        os.makedirs(local_npzdir, exist_ok=True)

        def _resolve_src(src):
            if os.path.isabs(src):
                return os.path.realpath(src)
            if root_dir:
                cand = os.path.realpath(os.path.join(root_dir, src))
                if os.path.exists(cand):
                    return cand
            raise FileNotFoundError(
                f"Relative npz path '{src}' with no 'root_dir' in manifest. "
                "Rebuild the manifest with absolute paths."
            )

        # preflight trace
        some_src_in = next(iter(npzmap.values()))
        print(f"[STAGE] example in manifest: {some_src_in}", flush=True)
        print(f"[STAGE] resolves to       : {_resolve_src(some_src_in)}", flush=True)

        def _stage_one(item):
            s, src_in = item
            src = _resolve_src(src_in)
            if not os.path.exists(src):
                raise FileNotFoundError(f"Staging failed: {src_in} → {src}")
            dst = os.path.join(local_npzdir, os.path.basename(src))
            try:
                if os.path.exists(dst) and os.path.getmtime(dst) >= os.path.getmtime(src):
                    return s, dst
            except Exception:
                pass
            try:
                if os.path.exists(dst):
                    os.remove(dst)
                os.link(src, dst)   # fast if same FS
            except OSError:
                shutil.copy2(src, dst)
            return s, dst

        with ThreadPoolExecutor(max_workers=threads) as ex:
            npzmap = dict(ex.map(_stage_one, npzmap.items()))

        # ---------- choose sentinel baseline in the working domain ----------
        if params.log1p:
            M = np.float32(math.log1p(params.M_raw))
            def to_residual(xf32):
                return np.log1p(xf32).astype(np.float32, copy=False) - M
        else:
            M   = float(params.M_raw)
            inv = np.float32(1.0 / M)
            def to_residual(xf32):
                return (xf32 * inv) - np.float32(1.0)  # sentinel -> 0; observed <= 0

        # ---------- on-disk pid→col map (pid as TEXT to avoid uint64 overflow) ----------
        db_path = os.path.join(os.getcwd(), "pidmap.sqlite")
        try:
            os.remove(db_path)
        except FileNotFoundError:
            pass
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur.execute("PRAGMA journal_mode=OFF;")
        cur.execute("PRAGMA synchronous=OFF;")
        cur.execute("CREATE TABLE map (col INTEGER PRIMARY KEY AUTOINCREMENT, pid TEXT UNIQUE);")
        cur.execute("CREATE INDEX IF NOT EXISTS idx_map_pid ON map(pid);")
        con.commit()

        def _u64_to_str_list(arr):
            # arr: np.ndarray of np.uint64
            return [str(int(x)) for x in arr]

        # ---------- Pass 1: count nnz and register unique pids ----------
        t0 = time.time()
        nnz_total = 0
        _last = t0

        # (optional) wrap inserts into an explicit transaction for speed
        cur.execute("BEGIN")

        for i, s in enumerate(samples, 1):
            z = np.load(npzmap[s], mmap_mode="r")
            pids = z["pair_id"]                # np.uint64 (sorted)
            nnz_total += pids.size

            up = np.unique(pids)
            CH = 10000
            for start in range(0, up.size, CH):
                chunk = up[start:start+CH]
                pid_strs = [str(int(x)) for x in chunk]
                cur.executemany("INSERT OR IGNORE INTO map(pid) VALUES (?)",
                                [(ps,) for ps in pid_strs])

            now = time.time()
            if (i == 1) or (i == len(samples)) or (i % 25 == 0) or (now - _last > 5):
                print(f"[CSR-P1] {i:,}/{len(samples):,} | nnz~{nnz_total:,} "
                      f"| elapsed {now - t0:,.1f}s", flush=True)
                _last = now
        con.commit()
        print(f"[CSR-P1] done in {time.time()-t0:.1f}s | nnz={nnz_total:,}", flush=True)

        # columns = number of unique pair_ids
        cur.execute("SELECT COUNT(*) FROM map;")
        n_cols = int(cur.fetchone()[0])
        N = len(samples)
        print(f"[CSR] rows={N:,} cols={n_cols:,} nnz={nnz_total:,}", flush=True)

        # ---------- Preallocate CSR ----------
        data    = np.empty(nnz_total, dtype=np.float32)
        indices = np.empty(nnz_total, dtype=np.int32)
        indptr  = np.empty(N + 1, dtype=np.int64)
        indptr[0] = 0

        # helper: bulk pid→col with caching + IN() chunks (pid TEXT)
        def map_pids_to_cols(pids, cache):
            out = np.empty(pids.size, dtype=np.int32)
            CH = 10000
            pos = 0
            for start in range(0, pids.size, CH):
                chunk = pids[start:start+CH]
                pid_strs = [str(int(x)) for x in chunk]
                need = [ps for ps in pid_strs if ps not in cache]
                if need:
                    qmarks = ",".join("?" for _ in need)
                    rows = cur.execute(f"SELECT pid,col FROM map WHERE pid IN ({qmarks})", need).fetchall()
                    cache.update({pid: col for pid, col in rows})
                out[pos:pos+len(pid_strs)] = [cache[ps] for ps in pid_strs]
                pos += len(pid_strs)
            return out

        # ---------- Pass 2: fill CSR arrays ----------
        t1 = time.time()
        offset = 0
        pid_cache = {}
        for i, s in enumerate(samples, 1):
            z    = np.load(npzmap[s], mmap_mode="r")
            ids  = z["pair_id"]
            vals = z["dist"].astype(np.float32, copy=False)

            cols_i = map_pids_to_cols(ids, pid_cache)
            res    = to_residual(vals)

            k = ids.size
            data[offset:offset+k]    = res
            indices[offset:offset+k] = cols_i
            indptr[i] = indptr[i-1] + k
            offset += k

            if (i == 1) or (i % 250 == 0) or (i == len(samples)):
                print(f"[CSR-P2] {i:,}/{len(samples):,} | fill {indptr[i]:,}/{nnz_total:,}", flush=True)

        X = csr_matrix((data, indices, indptr), shape=(N, n_cols))
        del data, indices, indptr; gc.collect()

        # ---------- Save outputs ----------
        os.makedirs(os.path.dirname(output.X), exist_ok=True)
        save_npz(output.X, X, compressed=True)
        pd.Series(samples, name="sample").to_csv(output.rows, sep="\t", index=False)

        # save column pid (as TEXT) in column order
        pid_in_order = [pid for (pid,) in cur.execute("SELECT pid FROM map ORDER BY col")]
        pd.Series(pid_in_order, name="pair_id").to_csv(output.cols, sep="\t", index=False)

        con.close()
        print(f"[CSR] wrote {output.X}", flush=True)

rule odog_umap_full_sentinel:
    """
    This gets up to 15GB RES, 23 VIRT so far while doing this on 5800 genomes.
    """
    input:
        X         = results_base_directory + "/allsamples.sentinel.csr.npz",
        rows      = results_base_directory + "/allsamples.rows.tsv",
        cols      = results_base_directory + "/allsamples.cols.tsv",   # <- add this
        sampletsv = results_base_directory + "/sampledf.tsv"
    output:
        df   = results_base_directory + "/all_samples_umaps/allsamples.umap.neighbors_{n}.mind_{m}.euclidean.df"
    params:
        random_state = config["umap"].get("random_state", 42),
        densmap      = config["umap"].get("densmap", False),
        low_memory   = True
    threads: 8
    retries: 4
    shadow: "copy-minimal"
    resources:
        mem_mb  = odogPlotUMAP_get_mem_mb,
        runtime = odogPlotUMAP_get_runtime
    run:
        # --- THREAD LIMITS: set *before* importing numpy/scipy/umap ---
        import os
        for v in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "BLIS_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "NUMBA_NUM_THREADS",
        ):
            os.environ[v] = str(max(1, threads - 1))
        os.environ.setdefault("NUMBA_THREADING_LAYER", "omp")  # or "tbb" if that's your stack

        # optional: avoid OpenMP fork issues on some clusters
        os.environ.setdefault("KMP_INIT_AT_FORK", "FALSE")
        os.environ.setdefault("KMP_AFFINITY", "disabled")

        import time, math
        import numpy as np
        import pandas as pd
        from scipy.sparse import load_npz, issparse

        # try to clamp BLAS inside fit() too (safe if unavailable)
        try:
            from threadpoolctl import threadpool_limits, threadpool_info
        except Exception:
            threadpool_limits = None
            threadpool_info = None

        # also ask numba to honor our thread count (in addition to env)
        try:
            import numba
            numba.set_num_threads(threads)
        except Exception:
            pass

        # Defer umap import until after env + numba setup
        import umap
        from PhyloTreeUMAP import umap_mapper_to_df

        # --- small helper for memory debug ---
        def _rss_mb():
            try:
                with open("/proc/self/statm") as f:
                    pages = int(f.read().split()[1])
                return (pages * os.sysconf("SC_PAGE_SIZE")) / (1024**2)
            except Exception:
                return float("nan")

        print(f"[UMAP-FULL] threads reservation={threads}", flush=True)
        try:
            import numba
            print(f"[UMAP-FULL] numba layer={getattr(numba, 'threading_layer', lambda: 'n/a')()} "
                  f"num_threads={getattr(numba, 'get_num_threads', lambda: lambda: 'n/a')()}", flush=True)
        except Exception:
            pass
        if threadpool_info:
            try:
                infos = [(lib.get('internal_api'), lib.get('num_threads')) for lib in threadpool_info()]
                print(f"[UMAP-FULL] BLAS libs: {infos}", flush=True)
            except Exception:
                pass

        # --- load matrix and sample order ---
        t_load0 = time.time()
        X = load_npz(input.X)
        t_load1 = time.time()
        print(f"[UMAP-FULL] loaded CSR in {t_load1 - t_load0:.1f}s | RSS≈{_rss_mb():.0f} MB", flush=True)
        if not issparse(X):
            raise ValueError("X is not sparse; expected CSR.")

        nnz = X.nnz
        N, C = X.shape
        density = (nnz / (N * C)) if N*C else float("nan")
        print(f"[UMAP-FULL] X shape={X.shape} | nnz={nnz:,} | density={density:.6f} | dtype={X.dtype}", flush=True)

        samples = pd.read_csv(input.rows, sep="\t")["sample"].tolist()

        # --- quick file consistency checks ---
        def _count_lines(path):
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                return sum(1 for _ in f)

        n_cols_file = _count_lines(input.cols) - 1
        n_rows_file = _count_lines(input.rows) - 1
        if X.shape[1] != n_cols_file:
            raise ValueError(f"Column mismatch: X has {X.shape[1]} cols but cols.tsv has {n_cols_file} rows.")
        if X.shape[0] != n_rows_file:
            raise ValueError(f"Row mismatch: X has {X.shape[0]} rows but rows.tsv has {n_rows_file} rows.")

        # --- UMAP params (we'll inject precomputed kNN to avoid RP-forest) ---
        n_neighbors = int(wildcards.n)
        min_dist    = float(wildcards.m)
        if n_neighbors >= X.shape[0]:
            raise ValueError(f"n_neighbors ({n_neighbors}) must be < n_samples ({X.shape[0]}).")

        # Build exact kNN with sklearn on the CSR (algorithm='brute' keeps it sparse-friendly)
        print("[UMAP-FULL] computing exact kNN with sklearn (brute, euclidean) ...", flush=True)
        t_knn0 = time.time()
        try:
            from sklearn.neighbors import NearestNeighbors
            # n_jobs uses OpenMP under the hood for pairwise distances; be conservative
            nn = NearestNeighbors(
                n_neighbors=n_neighbors,
                algorithm="brute",
                metric="euclidean",
                n_jobs=max(1, threads - 1),
            )
            nn.fit(X)
            dists, indices = nn.kneighbors(X, return_distance=True)
            dists   = dists.astype(np.float32, copy=False)
            indices = indices.astype(np.int32,  copy=False)
            print(f"[UMAP-FULL] kNN done in {time.time()-t_knn0:.1f}s", flush=True)
            print(f"[UMAP-FULL] kNN shapes: indices={indices.shape}, dists={dists.shape} "
                  f"| min/max dist=({float(dists.min()):.4g}, {float(dists.max()):.4g})",
                  flush=True)

        except Exception as e:
            # Fallback: PyNNDescent but *without* RP-tree init
            print(f"[UMAP-FULL] sklearn brute kNN failed: {e} ; falling back to NNDescent(rp_tree_init=False)", flush=True)
            from pynndescent import NNDescent
            t_knn0 = time.time()
            nn = NNDescent(
                X,
                n_neighbors=n_neighbors,
                metric="euclidean",
                rp_tree_init=False,   # <-- critical: avoid RP forest
                n_trees=0,
                leaf_size=15,
                low_memory=True,
                #random_state=int(params.random_state), commenting this out to stop numba from using many cores
            )
            indices, dists = nn.neighbor_graph
            dists   = dists.astype(np.float32, copy=False)
            indices = indices.astype(np.int32,  copy=False)
            print(f"[UMAP-FULL] NNDescent kNN done in {time.time()-t_knn0:.1f}s", flush=True)

        # ------ checks HERE (covers both branches) ------
        if not np.isfinite(dists).all():
            # optional: show a few offending positions
            bad = np.argwhere(~np.isfinite(dists))
            raise ValueError(f"Non-finite distances detected in kNN results (example indices: {bad[:5].tolist()}).")

        if indices.shape != dists.shape or indices.shape[1] != n_neighbors or indices.shape[0] != X.shape[0]:
            raise ValueError(f"kNN shape mismatch: indices {indices.shape}, dists {dists.shape}, "
                             f"expected ({X.shape[0]}, {n_neighbors}).")

        print(f"[UMAP-FULL] kNN shapes: indices={indices.shape}, dists={dists.shape} | "
              f"min/max/mean dist=({float(dists.min()):.4g}, {float(dists.max()):.4g}, {float(dists.mean()):.4g})",
              flush=True)
        # ----------------------------------------------

        # Construct reducer; keep sparse-friendly settings, but we won't let it find neighbors
        kwargs = dict(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric="euclidean",
            densmap=bool(params.densmap),
            init="random",
            random_state=int(params.random_state),
            low_memory=bool(params.low_memory),
            verbose=True,
        )
        reducer = umap.UMAP(**kwargs)

        # Pass precomputed neighbors directly if this UMAP version supports it
        fit_accepts_knn = "knn" in reducer.fit.__code__.co_varnames
        print(f"[UMAP-FULL] reducer.fit accepts knn={fit_accepts_knn}", flush=True)

        t0 = time.time()
        if fit_accepts_knn:
            print(f"[UMAP-FULL] starting fit with precomputed kNN | RSS≈{_rss_mb():.0f} MB", flush=True)
            mapper = reducer.fit(X, knn=(indices, dists))
        else:
            # Robust fallback: build the fuzzy graph ourselves and hand it to UMAP internals
            print("[UMAP-FULL] UMAP version lacks fit(..., knn=...). Using manual fuzzy graph.", flush=True)
            from umap.umap_ import fuzzy_simplicial_set
            rng = np.random.RandomState(int(params.random_state))
            # Note: we give the raw data X and our kNN so UMAP won’t try to re-run NN search
            graph, sigmas, rhos = fuzzy_simplicial_set(
                X,
                n_neighbors=n_neighbors,
                random_state=rng,
                metric="euclidean",
                knn_indices=indices,
                knn_dists=dists,
                set_op_mix_ratio=1.0,
                local_connectivity=1.0,
            )
            # Prime reducer state and embed
            reducer._raw_data = X
            reducer._knn_indices = indices
            reducer._knn_dists   = dists
            reducer.graph_       = graph
            # use the public API to finish the embedding
            mapper = reducer.fit(X)
        print(f"[UMAP-FULL] fit done in {time.time()-t0:.1f}s | RSS≈{_rss_mb():.0f} MB", flush=True)

        # --- write output ---
        t_df0 = time.time()
        cdf = pd.read_csv(input.sampletsv, sep="\t", index_col=0)
        umap_df = umap_mapper_to_df(mapper, cdf)
        os.makedirs(os.path.dirname(output.df), exist_ok=True)
        umap_df.to_csv(output.df, sep="\t", index=True)
        t_df1 = time.time()
        print(f"[UMAP-FULL] wrote {output.df} in {t_df1 - t_df0:.1f}s "
              f"| rows={len(umap_df):,} cols={umap_df.shape[1]:,}", flush=True)


##    ┓     ┏┓┓ ┏┓┳┓┏┓┏┓
## ┏┓┏┫┏┓┏┓ ┃ ┃ ┣┫┃┃┣ ┗┓ - One-Dot-One-Genome plots FOR SPECIFIC CLADES
## ┗┛┗┻┗┛┗┫ ┗┛┗┛┛┗┻┛┗┛┗┛   Each dot represents a single genome, and the data vector is the distance pairs
##        ┛
#rule odogClCooGen:
#    """
#    This generates a coo matrix of the distance matrices for the genomes in this clade.
#
#    **This is a modification of the rule called odogCooGen.**
#      Each row is a genome, and each column is the distance between every locus pair.
#      There is no pre-processing performed on the distance matrices. These are the raw
#      distances between every locus pair.
#
#    This takes up a lot of RAM and time. For 3600 genomes, it took:
#        CPU Efficiency: 97.32% of 01:57:32 core-walltime
#        Job Wall-clock time: 01:57:32
#        Memory Utilized: 216.61 GB
#        Memory Efficiency: 110.91% of 195.31 GB
#    """
#    input:
#        gbgzfiles = expand(results_base_directory + "/distance_matrices/{sample}.gb.gz",
#                           sample = config["sample_to_rbh_file"].keys()),
#        combotoindex = results_base_directory + "/combo_to_index.txt",
#        sampletsv    = results_base_directory + "/sampledf.tsv"
#    output:
#        coo       = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.coo.npz",
#        sampletsv = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.sampledf.tsv"
#    threads: 1
#    retries: 6
#    resources:
#        mem_mb = coo_get_mem_mb,
#        highio = 1,
#        bigUMAPSlots = 1,
#        runtime = 300
#    params:
#        taxids_to_keep    = lambda wildcards: config["taxids"][wildcards.taxanalysis][0],
#        taxids_to_remove  = lambda wildcards: config["taxids"][wildcards.taxanalysis][1]
#    run:
#        # filter the sample df
#        sampledf = pd.read_csv(input.sampletsv, sep = "\t", index_col = 0)
#        sampledf = filter_sample_df_by_clades(sampledf, params.taxids_to_keep, params.taxids_to_remove)
#        # Reset the index. We need to reset the indices so that we can associate the annotated sampledf with the
#        #   metadata-less coo matrix.
#        sampledf = sampledf.reset_index(inplace = False)
#        # save the sampledf. Keep the column.
#        sampledf.to_csv(output.sampletsv, sep = "\t")
#        print(sampledf)
#
#        # read in the combo_to_index file as a df, then convert to a dict
#        alg_combo_to_ix = algcomboix_file_to_dict(input.combotoindex)
#        coo = construct_coo_matrix_from_sampledf(sampledf, alg_combo_to_ix)
#        save_npz(output.coo, coo)
#
#
#def odogPlotCladeUMAP_get_mem_mb(wildcards, attempt):
#    """
#    The amount of RAM needed for the script depends on the size of the input genome.
#    """
#    attemptdict = {1: 20000,
#                   2: 40000,
#                   3: 80000,
#                   4: 160000}
#    return attemptdict[attempt]
#
#def odogPlotCladeUMAP_get_runtime(wildcards, attempt):
#    """
#    The amount of RAM needed for the script depends on the size of the input genome.
#    """
#    attemptdict = {1: 5,
#                   2: 10,
#                   3: 20,
#                   4: 40}
#    return attemptdict[attempt]
#
#rule odogCl_calculate_umap:
#    input:
#        sampletsv    = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.sampledf.tsv",
#        combotoindex = results_base_directory + "/combo_to_index.txt",
#        coo          = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.coo.npz",
#    output:
#        df   = results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df"
#    threads: 1
#    #retries: 3
#    resources:
#        mem_mb  = odogPlotCladeUMAP_get_mem_mb,
#        runtime = odogPlotCladeUMAP_get_runtime,
#        bigUMAPSlots = 1
#    run:
#        # note 20250612 - updated this method to just generate the df. All plotting is handled elsewhere.
#        mgt_mlt_umap(input.sampletsv, input.combotoindex, input.coo,
#                     wildcards.sizeNaN, int(wildcards.n), float(wildcards.m),
#                     output.df, missing_value_as = 9999999999)
#
#rule odogCl_pairwise_distance:
#    input:
#        sampletsv    = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.sampledf.tsv",
#        combotoindex = results_base_directory + "/combo_to_index.txt",
#        coo          = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.coo.npz",
#    output:
#        matrix = results_base_directory + "/ODOG/clades/{taxanalysis}.pairwise_distance.missing_{sizeNaN}.tsv"
#    threads: 1
#    resources:
#        mem_mb  = odogPlotCladeUMAP_get_mem_mb,
#        runtime = odogPlotCladeUMAP_get_runtime,
#        bigUMAPSlots = 1
#    run:
#        odog_pairwise_distance_matrix(input.sampletsv, input.combotoindex, input.coo,
#                                      wildcards.sizeNaN, output.matrix,
#                                      missing_value_as = 9999999999)
#
#def odol_plot_get_mem_mb(wildcards, attempt):
#    """
#    The amount of RAM needed for the script depends on the size of the input genome.
#    """
#    attemptdict = {1: 2000,
#                   2: 30000,
#                   3: 50000,
#                   4: 100000,
#                   5: 200000,
#                   6: 300000,
#                   7: 400000}
#    return attemptdict[attempt]
#
#def odol_plot_get_runtime(wildcards, attempt):
#    """
#    The amount of RAM needed for the script depends on the size of the input genome.
#    """
#    attemptdict = {1: 5,
#                   2: 120,
#                   3: 150,
#                   4: 180,
#                   5: 210,
#                   6: 240,
#                   7: 270}
#    return attemptdict[attempt]
#
#rule odogCl_GenHTML:
#    """
#    This makes an interactive .html file of the clade-specific UMAP plot.
#    """
#    input:
#        df   = results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df"
#    output:
#        html = results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html"
#    params:
#        outdir = results_base_directory + "/subchrom",
#    retries: 6
#    threads: 1
#    resources:
#        mem_mb = odol_plot_get_mem_mb,
#        highio = 1,
#        runtime = odol_plot_get_runtime
#    run:
#        plottitle = f"MGT plot of {wildcards.taxanalysis}, n = {wildcards.n} m={wildcards.m} NaN Size={wildcards.sizeNaN}"
#        mgt_mlt_plot_HTML(input.df, output.html, plot_title=plottitle, analysis_type = "MGT")
#
#rule odogClPDF:
#    input:
#        df  = results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
#    output:
#        pdf = results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.phylogeny.pdf"
#    threads: 1
#    retries: 3
#    resources:
#        mem_mb = pdf_get_mem_mb,
#        highio = 1,
#        runtime = 10
#    run:
#        plot_umap_phylogeny_pdf(input.df, output.pdf, wildcards.taxanalysis, wildcards.sizeNaN, wildcards.n, wildcards.m)
#
#rule odogClSweep:
#    """
#    Makes a pdf of the parameter sweep of the all-species UMAP plots.
#    """
#    input:
#        dfs = expand(results_base_directory + "/ODOG/clades/{{taxanalysis}}.neighbors_{n}.mind_{m}.missing_{{sizeNaN}}.df",
#                n = codog_n,
#                m = codog_m),
#        plotdfs = os.path.join(snakefile_path, "PhyloTreeUMAP_plotdfs.py")
#    output:
#        pdf    = results_base_directory + "/ODOG/clades/{taxanalysis}.missing_{sizeNaN}.paramsweep.pdf",
#        html   = results_base_directory + "/ODOG/clades/{taxanalysis}.missing_{sizeNaN}.paramsweep.html"
#    params:
#        prefix = results_base_directory + "/ODOG/clades/{taxanalysis}.missing_{sizeNaN}.paramsweep"
#    threads: 1
#    retries: 4
#    resources:
#        mem_mb = pdf_get_mem_mb,
#        highio = 1,
#        runtime = 5
#    shell:
#        """
#        python {input.plotdfs} -f "{input.dfs}" -p {params.prefix} --pdf --html
#        """
#
