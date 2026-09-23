#!/opt/bin/python3

"""
When a user uploads "three-tab" file  bundle it contains the following files:

  - expression.tab
  - genes.tab
  - observations.tab

It's expected that the values in the 1st column of the observations tab file
match positionally with the values in the row headers of the expression.tab
file.  If not, you just get this error from scanpy:

  Exception: Observation IDs from 'expressions' and 'observations' files are not the same.

Doesn't help you find where the issue is.  This script reports where they first disagree.
"""


import argparse
import difflib
import sys

def main():
    parser = argparse.ArgumentParser( description="Checks that observation IDs match in the expression tab file")
    parser.add_argument('-o', '--observations_file', type=str, required=True, help='Path to an input observations tab file' )
    parser.add_argument('-e', '--expression_file', type=str, required=True, help='Path to an input expression file' )
    parser.add_argument('-w', '--context_window', type=int, default=3, help='Rows of surrounding context to show around each mismatch (default: 3)' )
    parser.add_argument('-m', '--max_mismatch_reports', type=int, default=3, help='Maximum number of mismatch regions to report (default: 3)' )
    args = parser.parse_args()

    if args.context_window < 0:
        print("--context_window must be >= 0", file=sys.stderr)
        sys.exit(2)

    if args.max_mismatch_reports < 1:
        print("--max_mismatch_reports must be >= 1", file=sys.stderr)
        sys.exit(2)

    obs_ids = get_observation_ids(args.observations_file)

    obs_columns = get_expression_observation_columns(args.expression_file)

    print("Got {0} obs IDs from {1}".format(len(obs_ids), args.observations_file))
    print("Got {0} obs IDs from {1}".format(len(obs_columns), args.expression_file))

    mismatch_indices = get_mismatch_indices(obs_ids, obs_columns)

    if not mismatch_indices:
        print("All observation values matched exactly")
        return

    print("Found {0} mismatched position(s).".format(len(mismatch_indices)), file=sys.stderr)
    print_similarity_report(obs_ids, obs_columns)
    print_shift_diagnostic(obs_ids, obs_columns)

    context_window = args.context_window
    max_context_reports = args.max_mismatch_reports
    print("\nContext around first {0} mismatch(es):".format(min(max_context_reports, len(mismatch_indices))), file=sys.stderr)

    for idx in mismatch_indices[:max_context_reports]:
        print_mismatch_context(obs_ids, obs_columns, idx, context_window)

    sys.exit(1)


def get_expression_observation_columns(infile):
    with open(infile) as fh:
        for line in fh:
            # only the first line matters
            return line.rstrip("\n\r").split("\t")[1:]

    return list()


def get_mismatch_indices(obs_ids, obs_columns):
    mismatch_indices = list()
    max_len = max(len(obs_ids), len(obs_columns))

    for i in range(max_len):
        of_id = obs_ids[i] if i < len(obs_ids) else None
        ef_id = obs_columns[i] if i < len(obs_columns) else None

        if of_id != ef_id:
            mismatch_indices.append(i)

    return mismatch_indices


def print_similarity_report(obs_ids, obs_columns):
    max_len = max(len(obs_ids), len(obs_columns))
    positional_matches = 0

    if max_len > 0:
        for i in range(max_len):
            of_id = obs_ids[i] if i < len(obs_ids) else None
            ef_id = obs_columns[i] if i < len(obs_columns) else None

            if of_id == ef_id:
                positional_matches += 1

    positional_ratio = (positional_matches / max_len) if max_len else 1.0
    set_union = set(obs_ids) | set(obs_columns)
    set_intersection = set(obs_ids) & set(obs_columns)
    set_overlap_ratio = (len(set_intersection) / len(set_union)) if set_union else 1.0
    sequence_ratio = difflib.SequenceMatcher(a=obs_ids, b=obs_columns).ratio()

    if positional_ratio >= 0.95 and len(obs_ids) == len(obs_columns):
        classification = "nearly identical"
    elif positional_ratio <= 0.20 and set_overlap_ratio <= 0.20:
        classification = "completely off"
    else:
        classification = "partially matching"

    print("\nSimilarity summary:", file=sys.stderr)
    print("  Classification: {0}".format(classification), file=sys.stderr)
    print("  Positional match: {0:.2%} ({1}/{2})".format(positional_ratio, positional_matches, max_len), file=sys.stderr)
    print("  Set overlap: {0:.2%} ({1}/{2} unique IDs)".format(set_overlap_ratio, len(set_intersection), len(set_union)), file=sys.stderr)
    print("  Sequence similarity: {0:.2%}".format(sequence_ratio), file=sys.stderr)

    if set(obs_ids) == set(obs_columns) and obs_ids != obs_columns:
        print("  Note: Same IDs appear in both files, but order differs.", file=sys.stderr)


def print_shift_diagnostic(obs_ids, obs_columns):
    matcher = difflib.SequenceMatcher(a=obs_ids, b=obs_columns)
    opcodes = matcher.get_opcodes()
    non_equal_ops = [op for op in opcodes if op[0] != "equal"]

    if len(non_equal_ops) != 1:
        return

    tag, i1, i2, j1, j2 = non_equal_ops[0]
    equal_total = sum(i2_ - i1_ for (tag_, i1_, i2_, j1_, j2_) in opcodes if tag_ == "equal")
    min_len = min(len(obs_ids), len(obs_columns))
    equal_ratio = (equal_total / min_len) if min_len else 0.0

    # Only emit this hint when the two lists are mostly the same except for one insertion/deletion.
    if equal_ratio < 0.90:
        return

    if tag == "insert":
        extra_ids = obs_columns[j1:j2]
        print("\nPossible global shift detected:", file=sys.stderr)
        print("  Expression IDs contain {0} extra entr{1} near position {2}.".format(len(extra_ids), "y" if len(extra_ids) == 1 else "ies", j1 + 1), file=sys.stderr)
        print("  This can shift all following columns even if the rest are aligned.", file=sys.stderr)
        print("  Extra expression IDs (first 5): {0}".format(", ".join(extra_ids[:5]) if extra_ids else "<none>"), file=sys.stderr)
    elif tag == "delete":
        extra_ids = obs_ids[i1:i2]
        print("\nPossible global shift detected:", file=sys.stderr)
        print("  Observation IDs contain {0} extra entr{1} near position {2}.".format(len(extra_ids), "y" if len(extra_ids) == 1 else "ies", i1 + 1), file=sys.stderr)
        print("  This can shift all following columns even if the rest are aligned.", file=sys.stderr)
        print("  Extra observation IDs (first 5): {0}".format(", ".join(extra_ids[:5]) if extra_ids else "<none>"), file=sys.stderr)


def print_mismatch_context(obs_ids, obs_columns, mismatch_idx, window):
    max_len = max(len(obs_ids), len(obs_columns))
    start = max(0, mismatch_idx - window)
    end = min(max_len, mismatch_idx + window + 1)

    print("\n- Mismatch at position {0}".format(mismatch_idx + 1), file=sys.stderr)

    for i in range(start, end):
        of_id = obs_ids[i] if i < len(obs_ids) else "<missing>"
        ef_id = obs_columns[i] if i < len(obs_columns) else "<missing>"
        marker = ">>" if i == mismatch_idx else "  "
        status = "MATCH" if of_id == ef_id else "DIFF "
        print("{0} {1:4d} [{2}] expression=({3}) observations=({4})".format(marker, i + 1, status, ef_id, of_id), file=sys.stderr)

def get_observation_ids(infile):
    ids = list()
    line_count = 0

    with open(infile) as fh:
        for line in fh:
            line_count += 1

            if line_count == 1:
                continue

            ids.append(line.rstrip("\n\r").split('\t')[0])

    return ids

if __name__ == '__main__':
    main()
