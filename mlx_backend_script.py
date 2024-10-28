import mlx.core as mx
import numpy as np
import pandas as pd
from datetime import datetime
import sys

def madhyper_process(prefix, min_wells=4):
    print("start load:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    bigmas = mx.array(np.loadtxt(prefix+'_bigmas.tsv', delimiter='\t', dtype=np.float32))
    bigmbs = mx.array(np.loadtxt(prefix+'_bigmbs.tsv', delimiter='\t', dtype=np.float32))
    mdh = mx.array(np.loadtxt(prefix+'_mdh.tsv', delimiter='\t', dtype=np.int32))
    rowinds_bigmas = mx.arange(bigmas.shape[0])
    rowinds_bigmbs = mx.arange(bigmbs.shape[0])

    # Apply the min_wells filter
    non_zero_counts_bigmas = mx.sum(bigmas > 0, axis=1)
    non_zero_counts_bigmbs = mx.sum(bigmbs > 0, axis=1)
    valid_rows_bigmas = mx.nonzero(non_zero_counts_bigmas > min_wells)[0]
    valid_rows_bigmbs = mx.nonzero(non_zero_counts_bigmbs > min_wells)[0]

    bigmas = bigmas[valid_rows_bigmas]
    bigmbs = bigmbs[valid_rows_bigmbs]
    rowinds_bigmas = rowinds_bigmas[valid_rows_bigmas]
    rowinds_bigmbs = rowinds_bigmbs[valid_rows_bigmbs]

    results = []
    chunk_size = 500
    total_chunks = bigmas.shape[0] // chunk_size
    b_total = mx.sum(bigmbs > 0, axis=1, keepdims=True)
    bigmbs = (bigmbs > 0).T.astype(mx.float32)
    print("start time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    for ch in range(0, bigmas.shape[0], chunk_size):
        percent_complete = int((ch // chunk_size + 1) / total_chunks * 100)
        if percent_complete % 10 == 0 and percent_complete > 0:
            print(f'Progress: {ch} ({percent_complete}%)')

        chunk_end = min(ch + chunk_size, bigmas.shape[0])
        row_range = slice(ch, chunk_end)
        a_total = mx.sum(bigmas[row_range] > 0, axis=1, keepdims=True)
        overlaps = mx.matmul((bigmas[row_range] > 0).astype(mx.float32), bigmbs)
        mask_condition = -(overlaps.T - b_total).T < mdh[overlaps.astype(mx.int16), -(overlaps - a_total).astype(mx.int16)]
        pairs = mx.argwhere(mask_condition)

        result = {
            'alpha_nuc': 1 + rowinds_bigmas[row_range][pairs[:, 0]],
            'beta_nuc': 1 + rowinds_bigmbs[pairs[:, 1]],
            'wij': overlaps[pairs[:, 0], pairs[:, 1]],
            'wa': a_total[:, 0][pairs[:, 0]],
            'wb': b_total[:, 0][pairs[:, 1]]
        }
        results.append(result)

    print("end time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    results_df = pd.concat([pd.DataFrame(result) for result in results])
    results_df.to_csv(prefix+'_madhyperesults.csv', index=False)
    print(f"Number of pairs: {results_df.shape[0]}")

def correlation_process(prefix, min_wells=4):
    print("start load:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    bigmas = mx.array(np.loadtxt(prefix+'_bigmas.tsv', delimiter='\t', dtype=np.float32))
    bigmbs = mx.array(np.loadtxt(prefix+'_bigmbs.tsv', delimiter='\t', dtype=np.float32))
    rowinds_bigmas = mx.arange(bigmas.shape[0])
    rowinds_bigmbs = mx.arange(bigmbs.shape[0])
    print('file read done')

    # Apply the min_wells filter
    non_zero_counts_bigmas = mx.sum(bigmas > 0, axis=1)
    non_zero_counts_bigmbs = mx.sum(bigmbs > 0, axis=1)
    valid_rows_bigmas = mx.nonzero(non_zero_counts_bigmas > min_wells)[0]
    valid_rows_bigmbs = mx.nonzero(non_zero_counts_bigmbs > min_wells)[0]

    bigmas = bigmas[valid_rows_bigmas]
    bigmbs = bigmbs[valid_rows_bigmbs]
    rowinds_bigmas = rowinds_bigmas[valid_rows_bigmas]
    rowinds_bigmbs = rowinds_bigmbs[valid_rows_bigmbs]

    results = []
    chunk_size = 500
    total_chunks = bigmas.shape[0] // chunk_size
    bigmb_w1_scaled = bigmbs - mx.mean(bigmbs, axis=1, keepdims=True)
    bigmb_w1_scaled = (bigmb_w1_scaled / mx.linalg.norm(bigmb_w1_scaled, ord=2, axis=1, keepdims=True)).T
    b_total = mx.sum(bigmbs > 0, axis=1, keepdims=True)
    bigmbs = (bigmbs > 0).T.astype(np.float32)
    print("start processing time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    for ch in range(0, bigmas.shape[0], chunk_size):
        percent_complete = int((ch // chunk_size + 1) / total_chunks * 100)
        if percent_complete % 10 == 0 and percent_complete > 0:
            print(f'Progress: {ch} ({percent_complete}%)')

        chunk_end = min(ch + chunk_size, bigmas.shape[0])
        row_range = slice(ch, chunk_end)
        bigma_w1_scaled = bigmas[row_range] - mx.mean(bigmas[row_range], axis=1, keepdims=True)
        bigma_w1_scaled = bigma_w1_scaled / mx.linalg.norm(bigma_w1_scaled, ord=2, axis=1, keepdims=True)
        a_total = mx.sum(bigmas[row_range] > 0, axis=1, keepdims=True)
        pairwise_cors_method2 = mx.matmul(bigma_w1_scaled, bigmb_w1_scaled)
        overlaps = mx.matmul((bigmas[row_range] > 0).astype(np.float32), bigmbs)
        mask_condition = (-pairwise_cors_method2 <= mx.partition(pairwise_cors_method2 * -1, 3, axis=1)[:, 2:3])

        pairs = mx.argwhere(mask_condition)

        result = {
            'alpha_nuc': 1 + rowinds_bigmas[row_range][pairs[:, 0]],
            'beta_nuc': 1 + rowinds_bigmbs[pairs[:, 1]],
            'r': pairwise_cors_method2[pairs[:, 0], pairs[:, 1]],
            'wij': overlaps[pairs[:, 0], pairs[:, 1]],
            'wa': a_total[:, 0][pairs[:, 0]],
            'wb': b_total[:, 0][pairs[:, 1]]
        }
        results.append(result)

    print("end time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    results_df = pd.concat([pd.DataFrame(result) for result in results])
    results_df.to_csv(prefix + '_corresults.csv', index=False)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <prefix>")
        sys.exit(1)
        
    prefix = sys.argv[1]    
    madhyper_process(prefix)
    correlation_process(prefix)