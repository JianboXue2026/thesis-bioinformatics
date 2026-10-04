"""
fix_gff3_parent_refs.py
Fix invalid `Parent` attribute references in a GFF3 file.

The raw T2T-YAO annotation GFF3 contains records whose `Parent` points to IDs
that do not exist, which makes gff3ToGenePred fail with errors like:
    Can't find annotation record "ENSG00000282458.1_1" referenced by
    "ENST00000632292.1_1" Parent attribute

Strategy:
  1. collect all valid IDs
  2. find records whose Parent references a missing ID
  3. drop the orphaned parents from the attribute list (and drop records that
     only exist to reference them), then write a "_modified.gff3" file

Note: this script repairs the `Parent`-reference class of errors only — it does
not fix transcript-structure problems such as txStart >= txEnd or overlapping
exons (for those, use an AGAT-cleaned GFF3).
"""

# --- User-configurable section ----------------------------------------------
gff3_file_path = "/path/to/project/VEP_ref/T2T-Yao-hp.v1.1.gff3"
# -----------------------------------------------------------------------------


def parse_gff3(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    header = [line for line in lines if line.startswith("#")]
    entries = [line for line in lines if not line.startswith("#")]
    return header, entries


def extract_attributes(attribute_string):
    return {key: value for key, value in
            [attr.split("=") for attr in attribute_string.split(";")]}


def validate_gff3(entries):
    valid_ids = set()
    for line in entries:
        parts = line.strip().split("\t")
        attributes = parts[8]
        attr_dict = extract_attributes(attributes)
        if "ID" in attr_dict:
            valid_ids.add(attr_dict["ID"])

    errors = []
    for line in entries:
        parts = line.strip().split("\t")
        attributes = parts[8]
        attr_dict = extract_attributes(attributes)
        if "Parent" in attr_dict:
            parent_ids = attr_dict["Parent"].split(",")
            for pid in parent_ids:
                if pid not in valid_ids:
                    errors.append((parts, attr_dict, pid))
    return errors


def fix_gff3(entries, errors):
    fixed_entries = []
    error_ids = {error[2] for error in errors}
    for line in entries:
        parts = line.strip().split("\t")
        attributes = parts[8]
        attr_dict = extract_attributes(attributes)
        # Drop records whose ID is an orphaned parent
        if "ID" in attr_dict and attr_dict["ID"] in error_ids:
            continue
        if "Parent" in attr_dict:
            parent_ids = attr_dict["Parent"].split(",")
            valid_parent_ids = [pid for pid in parent_ids if pid not in error_ids]
            if not valid_parent_ids:
                del attr_dict["Parent"]
            else:
                attr_dict["Parent"] = ",".join(valid_parent_ids)
        attributes = ";".join([f"{key}={value}" for key, value in attr_dict.items()])
        parts[8] = attributes
        fixed_entries.append("\t".join(parts))
    return fixed_entries


def write_gff3(file_path, header, entries):
    with open(file_path, "w") as f:
        for line in header:
            f.write(line)
        for entry in entries:
            f.write(f"{entry}\n")


# --- Main --------------------------------------------------------------------
header, entries = parse_gff3(gff3_file_path)
errors = validate_gff3(entries)

if errors:
    print(f"Found {len(errors)} errors. Fixing them...")
    fixed_entries = fix_gff3(entries, errors)
    fixed_gff3_file_path = gff3_file_path.replace(".gff3", "_modified.gff3")
    write_gff3(fixed_gff3_file_path, header, fixed_entries)
    print(f"Fixed GFF3 file written to {fixed_gff3_file_path}")
else:
    print("No errors found in GFF3 file.")
