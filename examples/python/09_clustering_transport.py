"""Clustering, fingerprints ensembles, transport."""
from sci_form import (
    embed,
    embed_batch,
    butina_cluster,
    rmsd_matrix,
    filter_diverse,
    pack_conformers,
    split_worker_tasks,
    estimate_workers,
    ml_descriptors,
    charges,
)
from common import require_ok


def run():
    print("=== 09 clustering/transport ===")
    results = embed_batch(["CCCC"] * 10, seed=42)
    coords = [r.coords for r in results if r.is_ok()]
    print("coords count", len(coords))
    clusters = butina_cluster(coords, rmsd_cutoff=0.5)
    print("clusters", clusters.n_clusters, "sizes", clusters.cluster_sizes, "assignments", clusters.assignments[:5])

    mx = rmsd_matrix(coords[:4])
    print("rmsd matrix shape", len(mx), len(mx[0]))

    diverse_idx = filter_diverse(coords, rmsd_cutoff=1.0)
    print("diverse idx", diverse_idx)

    batch = pack_conformers([r for r in results if r.is_ok()])
    print("record batch rows", batch.num_rows, "columns", batch.column_names)
    print("float data keys", list(batch.float_data.keys()))

    tasks = split_worker_tasks(["CCO", "CCCC", "c1ccccc1"] * 3, n_workers=2, seed=42)
    print("tasks", len(tasks), tasks[0] if tasks else None)

    print("estimate_workers", estimate_workers(100, max_workers=8))

    # molecular descriptors
    from sci_form import ml_descriptors
    md = ml_descriptors("CCO")
    print("mol descriptors", md)


if __name__ == "__main__":
    run()
