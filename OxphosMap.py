import argparse
import pandas as pd
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex
from pymol import cmd

# ----- Function to check for missing columns --------------------

def check_required_columns(df, required, name):
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing column(s) in {name}: {', '.join(sorted(missing))}")

# ----- Function to convert hex color to RGB decimals --------------------

def hex_to_rgb(hex_color):
    """
    Convert a hex color string to an "RGB tuple" with floats from 0 to 1.
    """

    hex_color = hex_color.lstrip('#')
    
    if len(hex_color) != 6:
        raise ValueError(f"Invalid hex color format: '{hex_color}'")

    try:
        return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    
    except ValueError:
        raise ValueError(f"Invalid hex digits in color: '{hex_color}'")

# ----- Parser --------------------

def parse_arguments():
    
    parser = argparse.ArgumentParser(
        description="Process command-line arguments."
    )

    parser.add_argument(
        "-l",
        "--logfc_limit", 
        type=float, 
        required=True, 
        help="log2fc limit")
    
    parser.add_argument(
        "-p",
        "--file_path", 
        type=str, 
        required=True, 
        help="Path to dataset")

    parser.add_argument(
        "-s",
        "--structure_info_path", 
        type=str, 
        required=True, 
        help="Path to structure metadata")

    parser.add_argument(
        "-c",
        "--complex", 
        type=str, 
        required=True, 
        help="Which complex to visualize")

    return parser.parse_args()

# ----- Load data --------------------

def load_data(file_path: str, structure_info_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:

    """
    Loads data and checks for the required columns. 
    Then trims the data frames to only include those columns. 
    """

    fc_data = pd.read_csv(file_path)

    structure_info = pd.read_csv(structure_info_path)

    check_required_columns(fc_data, {"Log2FC", "GeneID"}, "LogFC data")    

    check_required_columns(structure_info, {"GeneID", "pdb", "chain", "complex"}, "Structure metadata")

    # Trim the data to avoid extra columns
    fc_subset = fc_data[["Log2FC", "GeneID"]].copy()

    # Ensure GeneID is string
    fc_subset["GeneID"] = fc_subset["GeneID"].astype("string")

    structure_info["GeneID"] = structure_info["GeneID"].astype("string")

    structure_metadata_subset = structure_info[["GeneID", "pdb", "chain", "complex"]].copy()

    return fc_subset, structure_metadata_subset

# ----- Filter and clip logFC --------------------

def filter_and_limit_log2fc(fc_data, structure_info, logfc_limit):
    """
    Filters fc_data by GeneID, ensuring only OXPHOS genes in the dataset.
    Also clips the logFC based on the input limit.
    """
    
    # Filter by gene IDs in structure_info to have only the relevant genes in the dataset
    mask = fc_data["GeneID"].isin(structure_info["GeneID"])

    # Apply filter and reset index
    fc_data_filtered = fc_data[mask].copy().reset_index(drop=True)

    # Clip log2FC column to avoid extreme values. This prevents outliers from messing up the color scales. 
    fc_data_filtered["Log2FC"] = fc_data_filtered["Log2FC"].clip(lower= -logfc_limit, upper=logfc_limit)

    return fc_data_filtered
    
# ----- Get the hex colors from logFC values --------------------

def convert_logfc_to_hex(fc_data_filtered):
    """
    Defines the colors of the structure based on the log2FCs in the dataset.
    Adds a new 'color' column with the corresponding hex values.
    """

    color_negative = "#003070"
    color_positive = "#700020"
    color_middle = "#FFFFFF"

    # Get the max absolute log2FC in the dataset to define the edges of the color scale below
    max_absolute = fc_data_filtered["Log2FC"].abs().max()

    # Create color map
    colormap = LinearSegmentedColormap.from_list("custom_gradient", [color_negative, color_middle, color_positive])

    # Create normalization function for coloring proteins
    norm = Normalize(vmin= -max_absolute, vmax = max_absolute)

    # Apply color map to the log2fc values
    fc_data_filtered["color"] = fc_data_filtered["Log2FC"].apply(lambda x: to_hex(colormap(norm(x))))

    return fc_data_filtered

# ----- merge structure metadata with logFC data containing the hex colors --------------------

def merge_and_color_nas(structure_info, fc_data_colors, na_color="#3E3E3E"):
    """
    Merges the structure metadata and the logFC data, to have a combined file to use for PyMOL visualization
    Colors any missing color values according to the input na_color
    """

    # Merge
    merged = structure_info.merge(fc_data_colors, on='GeneID', how='left')

    # Change all NA or empty color values to na_color
    mask = merged["Log2FC"].isna() | (merged["Log2FC"] == "")
    
    merged.loc[mask, "color"] = na_color

    return merged

# ----- Write the merged data to a csv file --------------------

def write_processed_data_to_csv(structure_merged, file_path):
    """
    Ads a '_processed' suffix to the logFC file path and saves a csv of the merged 
    logFC data and metadata. 
    """

    # Edit the input file name and add 'processed' for exporting the processed data 
    processed_path = file_path.with_name(file_path.stem + '_processed' + file_path.suffix)

    print(processed_path)

    structure_merged.to_csv(processed_path, index=False)

    print(f"\n\nProcessed data save as:\n{processed_path}\n\n")

# ----- Prepare data for PyMOL visalization --------------------

def prepare_pymol_visualization(structure_merged, complex_value):
    """
    Prepares the data for PyMOL visualization
    selects only the relevant complex based on the input
    gets the PDB ID, for fetching the structure
    """

    if not (structure_merged['complex'] == complex_value).any():
        raise ValueError(f"Complex '{complex_value}' not found in structure metadata")

    # Select only the specified complex
    structure_merged_complex = structure_merged[structure_merged["complex"] == complex_value]

    # Get the PBD ID of the protein structure from the data frame
    structure_id = structure_merged_complex["pdb"].unique()[0]

    # Print the information about the structure
    print(f"Visualizing complex {complex_value}\nUsing PDB ID: {structure_id}\n\n")

    return structure_merged_complex, structure_id

# ----- Color chains in PyMOL based on the Hex colors in the merged data frame --------------------

def color_pymol_chains(structure_merged_complex):
    """
    Color each chain based on the hex colors listed in the input data frame
    """
    
    for index, row in structure_merged_complex.iterrows():
        chain = row['chain']  # Get the chain name
        color = row['color']  # Get the hex color code
        
        # Generate the color name
        color_name = f"custom_color_{index}"
        
        # Color the specific chain using the RGB-converted hex code from the DataFrame
        rgb = hex_to_rgb(color)
        
        # Create the specific color in PyMol
        cmd.set_color(color_name, rgb)
        
        # Color the specific chain with the new color 
        cmd.color(color_name, chain)  

# ----- Main --------------------

def main():

    args = parse_arguments()

    file_path = Path(args.file_path)
    structure_info_path = Path(args.structure_info_path)
    logfc_limit = args.logfc_limit
    complex_value = args.complex

    if complex_value not in ["CI", "CII", "CIII", "CIV", "CV"]:
        raise ValueError("error: invalid complex name")

    # Print information about the input data
    print(f"\n\nStructure metadata:\n{structure_info_path}\n\nFile path:\n{file_path}\n\nLogFC limit:\n{logfc_limit}\n\n")

    # ----- Load data --------------------

    fc_data, structure_info = load_data(file_path, structure_info_path)

    # ----- Filter: select only oxphos genes --------------------

    fc_data_filtered = filter_and_limit_log2fc(fc_data, structure_info, logfc_limit)

    # ----- Convert logFCs to rgb values --------------------

    fc_data_colors = convert_logfc_to_hex(fc_data_filtered)

    # ----- Merge the datasets --------------------

    structure_merged = merge_and_color_nas(structure_info, fc_data_colors)

    # ----- Save csv file --------------------
    
    write_processed_data_to_csv(structure_merged, file_path)

    # ----- Prepare data for PyMOL --------------------

    structure_merged_complex, structure_id = prepare_pymol_visualization(structure_merged, complex_value)

    # ----- PyMOL commands --------------------

    # Clear the display
    cmd.delete("all")

    # Get the structure based on the PDB ID
    cmd.fetch(structure_id)

    # Remove all H2O molecules from the visualization
    cmd.remove("resn HOH")
    
    # Color chains based on the generated hex colors
    color_pymol_chains(structure_merged_complex)

main()