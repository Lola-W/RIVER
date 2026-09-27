"""Combine dataset support within methods and classify each candidate origin."""
from collections import defaultdict

from .models import Allele


def summarize(origins, samples, read_rows):
    """Return variant and summary records, consuming read rows only once."""
    supported = defaultdict(set)
    for row in read_rows:
        if row["qualifying_support"]:
            allele = Allele(row["chrom"], int(row["pos"]), row["ref"], row["alt"])
            supported[allele].add(row["sample_id"])
    cell_methods = sorted({sample.method for sample in samples.values() if sample.type == "single_cell" and sample.bam is not None})
    variants, groups = [], {}
    seen_origins = set()
    for origin in origins:
        if origin in seen_origins:
            continue
        seen_origins.add(origin)
        allele = origin.allele
        source = samples[origin.source_sample_id]
        datasets = sorted(supported.get(allele, set()))
        methods = sorted({samples[sid].method for sid in datasets})
        cells = {method: [sid for sid in datasets if samples[sid].method == method] for method in cell_methods}
        cell_counts = {method: len(values) for method, values in cells.items()}
        cell_classes = {method: "0" if count == 0 else "1" if count == 1 else ">=2" for method, count in cell_counts.items()}
        other = [method for method in methods if method != source.method]
        category = "0" if not other else "1" if len(other) == 1 else ">=2"
        variants.append({"chrom": allele.chrom, "pos": allele.pos, "ref": allele.ref, "alt": allele.alt,
                         "source_sample_id": source.sample_id, "source_method": source.method, "callset": origin.callset,
                         "supporting_samples": datasets, "supporting_methods": methods,
                         "supporting_cells_by_method": cells, "supporting_cell_counts_by_method": cell_counts,
                         "supporting_cell_classes_by_method": cell_classes,
                         "other_supporting_methods": other, "n_other_supporting_methods": len(other),
                         "other_method_support_class": category})
        key = (source.sample_id, source.method, origin.callset)
        group = groups.setdefault(key, {"source_sample_id": source.sample_id, "source_method": source.method,
                                      "callset": origin.callset, "n_variants": 0,
                                      "n_support_0": 0, "n_support_1": 0, "n_support_ge2": 0})
        group["n_variants"] += 1
        field = "n_support_0" if not other else "n_support_1" if len(other) == 1 else "n_support_ge2"
        group[field] += 1
    summaries = []
    for group in groups.values():
        for suffix in ("0", "1", "ge2"):
            group[f"proportion_support_{suffix}"] = group[f"n_support_{suffix}"] / group["n_variants"] if group["n_variants"] else None
        summaries.append(group)
    return variants, summaries
