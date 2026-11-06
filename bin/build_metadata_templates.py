#!/usr/bin/env python3

"""
Builds the metadata templates offered by the dataset uploader from the field definitions in
lib/gear/metadata_fields.py. Re-run this after changing those definitions.

Writes, to www/user_templates/ by default:
    metadata_template.xlsx          dataset-level metadata, for users to fill in
    metadata_example.xlsx           the same, filled in with the example values
    sample_metadata_template.xlsx   sample-level (observation) metadata columns, with requirements and descriptions
    metadata_fields.json            field rules read by the uploader form (www/js/upload_dataset.js)
"""

import argparse
import json
import os
import sys

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.worksheet.datavalidation import DataValidation

lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'lib')
sys.path.append(lib_path)

from gear import metadata_fields as mf

DEFAULT_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'www', 'user_templates')

# Placeholder contact for the example file, instead of the template's blank examples
EXAMPLE_OVERRIDES = {
    "contact_name": "Jane Doe",
    "contact_email": "demo@email.com",
    "contact_institute": "Institute for Genome Sciences",
}


def requirement_label(f: mf.MetadataField) -> str:
    if f.required:
        return "Y"
    if f.required_if:
        return "Y if " + ", ".join(f"{other} is '{value}'" for other, value in f.required_if.items()) + ", N if not"
    return "N"


def write_dataset_workbook(path: str, filled: bool) -> None:
    """Sheet 'metadata' with columns field | value | is_required | instructions | example (read by gear.metadata.Metadata)."""
    wb = Workbook()
    ws = wb.active
    ws.title = "metadata"
    ws.append(["field", "value", "is_required", "instructions", "example"])

    for f in mf.FIELDS:
        value = EXAMPLE_OVERRIDES.get(f.name, f.example) if filled else None
        ws.append([f.name, value, requirement_label(f), f.instructions, f.example])

        row = ws.max_row
        if f.choices:
            # Dataset types are shown by their user-facing labels (gosling as epigenome)
            choices = mf.TEMPLATE_DATASET_TYPES if f.name == "dataset_type" else f.choices
            dv = DataValidation(type="list", formula1='"{0}"'.format(",".join(choices)), allow_blank=not f.required)
            ws.add_data_validation(dv)
            dv.add(f"B{row}")
        elif f.max_length is not None and f.max_length < 65535:
            dv = DataValidation(type="textLength", operator="lessThanOrEqual", formula1=str(f.max_length),
                                error=f"Must be {f.max_length} characters or fewer.", showErrorMessage=True)
            ws.add_data_validation(dv)
            dv.add(f"B{row}")

    format_sheet(ws, widths={"A": 34, "B": 40, "C": 14, "D": 90, "E": 30})
    wb.save(path)


def write_sample_workbook(path: str) -> None:
    """
    Sheet 'sample_metadata' with one column per field. Under the column names are a row saying
    whether each column is required, a row describing it and a row with an example value.
    """
    columns = [(f.name, f.requirement, f.instructions, f.example) for f in mf.SAMPLE_FIELDS]
    columns.append(("(other columns)", "Optional", mf.SAMPLE_FIELD_NOTE, ""))

    wb = Workbook()
    ws = wb.active
    ws.title = "sample_metadata"
    for row in zip(*columns):
        ws.append(list(row))

    # The required and description rows are guidance, so set them apart from the example row
    guidance_fill = PatternFill(fill_type="solid", start_color="F2F2F2", end_color="F2F2F2")
    for cell in ws[1]:
        cell.font = Font(bold=True)
    for cell in ws[2]:
        cell.fill = guidance_fill
        cell.font = Font(bold=cell.value == "Required")
    for cell in ws[3]:
        cell.fill = guidance_fill
        cell.font = Font(italic=True)
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    for col_cells in ws.columns:
        ws.column_dimensions[col_cells[0].column_letter].width = 32
    ws.freeze_panes = "A4"
    wb.save(path)


def format_sheet(ws, widths: dict) -> None:
    for cell in ws[1]:
        cell.font = Font(bold=True)
    for letter, width in widths.items():
        ws.column_dimensions[letter].width = width
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "A2"


def write_fields_json(path: str) -> None:
    data = {
        "version": mf.METADATA_VERSION,
        "fields": [f.to_json() for f in mf.FIELDS],
        "max_keywords": mf.MAX_KEYWORDS,
        "max_keyword_length": mf.MAX_KEYWORD_LENGTH,
        # Widths of the dataset columns that hold only part of a field's value
        "annotation_release_max_length": mf.ANNOTATION_RELEASE_MAX_LENGTH,
        "geo_id_max_length": mf.GEO_ID_MAX_LENGTH,
    }
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")


def main():
    parser = argparse.ArgumentParser(description="Builds the uploader's metadata templates from lib/gear/metadata_fields.py")
    parser.add_argument('-o', '--output_dir', type=str, default=DEFAULT_OUTPUT_DIR, help='Directory to write the templates to')
    args = parser.parse_args()

    write_dataset_workbook(os.path.join(args.output_dir, "metadata_template.xlsx"), filled=False)
    write_dataset_workbook(os.path.join(args.output_dir, "metadata_example.xlsx"), filled=True)
    write_sample_workbook(os.path.join(args.output_dir, "sample_metadata_template.xlsx"))
    write_fields_json(os.path.join(args.output_dir, "metadata_fields.json"))


if __name__ == '__main__':
    main()
