import numpy as np
import os
import sys


def read_coeff_file(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()

    if not lines:
        raise RuntimeError(f"Empty file: {file_path}")

    if lines[0].strip().lower().startswith("sep="):
        header_line = lines[1].rstrip("\n")
        data_lines = lines[2:]
    else:
        header_line = lines[0].rstrip("\n")
        data_lines = lines[1:]

    header_fields = [field.strip() for field in header_line.split(";")]
    ncols = len(header_fields)

    rows = []

    for line in data_lines:
        if not line.strip():
            continue

        fields = [field.strip() for field in line.rstrip("\n").split(";")]

        if len(fields) != ncols:
            raise RuntimeError(
                f"Column mismatch in {file_path}: "
                f"expected {ncols}, got {len(fields)}.\n"
                f"Line was:\n{line}"
            )

        rows.append([float(field) for field in fields])

    data = np.asarray(rows, dtype=float)

    return header_fields, data


def format_number(value):
    if not np.isfinite(value):
        return ""

    # Decimal dot, no comma, no surrounding spaces.
    return f"{float(value):.15g}"


def post_process(foldername):
    batch_dir = os.path.join(foldername, "Batch_Results")

    files = []
    data_list = []
    headers = []

    for filename in os.listdir(batch_dir):
        if filename.startswith("aerodynamic_coefficients_") and filename.endswith(".csv"):
            file_path = os.path.join(batch_dir, filename)

            print("Opening", filename)

            header_fields, data = read_coeff_file(file_path)

            print(data.shape)

            if data.shape[0] == 0:
                continue

            files.append(filename)
            data_list.append(data)
            headers.append(header_fields)

    if not data_list:
        raise RuntimeError("No valid aerodynamic_coefficients_*.csv files found.")

    max_rows = max(data.shape[0] for data in data_list)

    selected_files = []
    selected_data = []
    selected_headers = []

    for filename, data, header_fields in zip(files, data_list, headers):
        if data.shape[0] == max_rows:
            selected_files.append(filename)
            selected_data.append(data)
            selected_headers.append(header_fields)

    reference_shape = selected_data[0].shape
    reference_header = selected_headers[0]

    for filename, data, header_fields in zip(selected_files, selected_data, selected_headers):
        if data.shape != reference_shape:
            raise RuntimeError(
                f"Shape mismatch in selected file {filename}: "
                f"expected {reference_shape}, got {data.shape}"
            )

        if header_fields != reference_header:
            raise RuntimeError(f"Header mismatch in selected file {filename}.")

    print("Maximum number of rows:", max_rows)
    print("Averaging files:")

    for filename in selected_files:
        print("  " + filename)

    data_averaged = np.mean(np.stack(selected_data, axis=0), axis=0)

    output_path = os.path.join(batch_dir, "aerodynamic_coefficients.csv")

    with open(output_path, "w", newline="") as f:
        f.write("sep=;\n")
        f.write(";".join(reference_header) + "\n")

        for row in data_averaged:
            f.write(";".join(format_number(value) for value in row) + "\n")

    print("Saved:", output_path)


if __name__ == "__main__":
    post_process(sys.argv[1])